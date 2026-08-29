use super::{ImportOptions, call, locus};
use crate::{Error, Result, model::*};
use std::{io::BufRead, str::FromStr};

fn below_min<T: FromStr + PartialOrd>(val: Option<&str>, min: Option<T>) -> bool {
    min.is_some_and(|min| {
        val.and_then(|txt| txt.parse::<T>().ok())
            .is_none_or(|val| val < min)
    })
}

pub(super) fn read(
    mut reader: Box<dyn BufRead>,
    options: &ImportOptions,
    builder: &mut DatasetBuilder,
) -> Result<()> {
    let mut line: String = String::new();
    let (sample, expected_cols) = loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            return Err(Error::InvalidInput("VCF #CHROM header is missing".into()));
        }
        builder.source.observe(&line)?;
        if line.starts_with("#CHROM\t") {
            let cols: Vec<_> = line.trim_end_matches(['\r', '\n']).split('\t').collect();
            if cols.len() < 10 || cols.get(8) != Some(&"FORMAT") {
                return Err(Error::InvalidInput("VCF must contain FORMAT/GT and at least one individual; site lists are not genotypes".into()));
            }
            let samples: Vec<SampleId> = cols[9..]
                .iter()
                .map(|name| name.parse())
                .collect::<Result<_>>()?;
            let sample: usize = options.sample.choose(&samples)?;
            builder.source.sample = Some(samples[sample].clone());
            break (sample + 9, cols.len());
        }
        if !line.starts_with('#') && !line.trim().is_empty() {
            return Err(Error::InvalidInput(
                "VCF data precedes the #CHROM header".into(),
            ));
        }
    };
    builder.source.notes.push("Only explicit diploid SNP GT calls are used; gVCF blocks are not expanded. Variant-only or imputed data can bias ROH estimates.".into());
    while {
        line.clear();
        reader.read_line(&mut line)? != 0
    } {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let row: Vec<_> = line.trim_end_matches(['\r', '\n']).split('\t').collect();
        let snp = if row.len() != expected_cols {
            Err(Rejection::Malformed)
        } else {
            locus(row[2], row[0], row[1], genotype(&row, sample, options))
        };
        builder.push(snp)?;
    }
    Ok(())
}

fn genotype(
    row: &[&str],
    sample: usize,
    options: &ImportOptions,
) -> std::result::Result<Genotype, Rejection> {
    if !options.include_filtered && !matches!(row[6], "PASS" | ".") {
        return Err(Rejection::Filtered);
    }
    if row[7]
        .split(';')
        .filter_map(|val| val.strip_prefix("END="))
        .any(|end| {
            end.parse::<u64>()
                .ok()
                .zip(row[1].parse::<u64>().ok())
                .is_some_and(|(end, start)| end > start)
        })
    {
        return Err(Rejection::ReferenceBlock);
    }
    let format: Vec<_> = row[8].split(':').collect();
    let vals: Vec<_> = row[sample].split(':').collect();
    let field = |name| {
        format
            .iter()
            .position(|key| *key == name)
            .and_then(|idx| vals.get(idx).copied())
    };
    if !options.include_filtered && field("FT").is_some_and(|val| !matches!(val, "PASS" | ".")) {
        return Err(Rejection::Filtered);
    }
    if field("LEN")
        .and_then(|val| val.parse::<u64>().ok())
        .is_some_and(|len| len > 1)
    {
        return Err(Rejection::ReferenceBlock);
    }
    if below_min(field("DP"), options.min_depth)
        || below_min(field("GQ"), options.min_genotype_quality)
    {
        return Err(Rejection::Filtered);
    }
    let genotype: &str = field("GT").ok_or(Rejection::Missing)?;
    let genotype: &str = genotype.strip_prefix(['/', '|']).unwrap_or(genotype);
    let indices: Vec<_> = genotype.split(['/', '|']).collect();
    if indices
        .iter()
        .any(|allele| *allele == "." || allele.is_empty())
    {
        return Err(Rejection::Missing);
    }
    if indices.len() != 2 {
        return Err(Rejection::Haploid);
    }
    if Base::parse(row[3]).is_none() {
        return Err(Rejection::NonSnp);
    }
    let alternates: Vec<_> = row[4].split(',').collect();
    if !alternates
        .iter()
        .any(|allele| Base::parse(allele).is_some() || matches!(*allele, "." | "<NON_REF>" | "<*>"))
    {
        return Err(Rejection::NonSnp);
    }
    let allele = |idx: &str| -> std::result::Result<&str, Rejection> {
        let idx = idx.parse::<usize>().map_err(|_| Rejection::Malformed)?;
        if idx == 0 {
            Ok(row[3])
        } else {
            alternates.get(idx - 1).copied().ok_or(Rejection::Malformed)
        }
    };
    call(&[allele(indices[0])?, allele(indices[1])?], CallCoding::Dna)
}
