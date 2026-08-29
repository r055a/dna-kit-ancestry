use super::{ImportOptions, call, containers, locus};
use crate::{Error, Result, model::*};
use std::{
    io::{BufRead, Read},
    path::Path,
};

struct Variant {
    chromosome: String,
    id: String,
    pos: String,
    alleles: Option<[String; 2]>,
}

fn variants(path: &Path, bim: bool, options: &ImportOptions) -> Result<Vec<Variant>> {
    let reader = containers::text_reader(path, options.max_bytes)?;
    let mut variants: Vec<Variant> = Vec::new();
    for (line, text) in reader.lines().enumerate() {
        let text: String = text?;
        if text.trim().is_empty() || text.starts_with('#') {
            continue;
        }
        let row: Vec<_> = text.split_whitespace().collect();
        if row.len() != if bim { 6 } else { 4 } {
            return Err(Error::InvalidInput(format!(
                "invalid PLINK metadata row {}",
                line + 1
            )));
        }
        variants.push(Variant {
            chromosome: row[0].into(),
            id: row[1].into(),
            pos: row[3].into(),
            alleles: if bim {
                Some([row[4].into(), row[5].into()])
            } else {
                None
            },
        });
    }
    Ok(variants)
}

fn sample_id(row: &[&str]) -> Result<SampleId> {
    format!("{}:{}", row[0], row[1]).parse()
}

fn read_samples(path: &Path, options: &ImportOptions) -> Result<Vec<SampleId>> {
    let mut samples: Vec<SampleId> = Vec::new();
    for txt in containers::text_reader(path, options.max_bytes)?.lines() {
        let txt: String = txt?;
        let row: Vec<_> = txt.split_whitespace().collect();
        if row.is_empty() {
            continue;
        }
        if row.len() != 6 {
            return Err(Error::InvalidInput("FAM/TFAM requires six columns".into()));
        }
        samples.push(sample_id(&row)?);
    }
    Ok(samples)
}

pub(super) fn read_tped(
    path: &Path,
    options: &ImportOptions,
    builder: &mut DatasetBuilder,
) -> Result<()> {
    let samples: Vec<SampleId> = read_samples(&containers::companion(path, "tfam")?, options)?;
    let chosen: usize = options.sample.choose(&samples)?;
    builder.source.sample = Some(samples[chosen].clone());
    for txt in containers::text_reader(path, options.max_bytes)?.lines() {
        let txt: String = txt?;
        if txt.trim().is_empty() {
            continue;
        }
        let row: Vec<_> = txt.split_whitespace().collect();
        if row.len() != 4 + 2 * samples.len() {
            return Err(Error::InvalidInput(
                "TPED allele count does not match TFAM individuals".into(),
            ));
        }
        builder.push(locus(
            row[1],
            row[0],
            row[3],
            call(&row[4 + 2 * chosen..6 + 2 * chosen], CallCoding::Auto),
        ))?;
    }
    Ok(())
}

fn bed_variant(variant: &Variant, code: u8, builder: &mut DatasetBuilder) -> Result<()> {
    let alleles: &[String; 2] = variant
        .alleles
        .as_ref()
        .ok_or_else(|| Error::InvalidInput("BIM alleles missing".into()))?;
    let genotype = match code {
        0 => call(&[&alleles[0], &alleles[0]], CallCoding::Auto),
        2 => call(&[&alleles[0], &alleles[1]], CallCoding::Auto),
        3 => call(&[&alleles[1], &alleles[1]], CallCoding::Auto),
        _ => Err(Rejection::Missing),
    };
    builder.push(locus(
        &variant.id,
        &variant.chromosome,
        &variant.pos,
        genotype,
    ))
}

fn skip(reader: &mut impl Read, bytes: usize) -> Result<()> {
    if std::io::copy(&mut reader.take(bytes as u64), &mut std::io::sink())? != bytes as u64 {
        return Err(Error::InvalidInput("truncated BED matrix".into()));
    }
    Ok(())
}

pub(super) fn read_ped(
    path: &Path,
    options: &ImportOptions,
    builder: &mut DatasetBuilder,
) -> Result<()> {
    let map: Vec<Variant> = variants(&containers::companion(path, "map")?, false, options)?;
    let mut reader = containers::text_reader(path, options.max_bytes)?;
    let mut samples: Vec<SampleId> = Vec::new();
    let mut line: String = String::new();
    while {
        line.clear();
        reader.read_line(&mut line)? != 0
    } {
        let prefix: Vec<_> = line.split_whitespace().take(6).collect();
        if prefix.is_empty() {
            continue;
        }
        if prefix.len() != 6 {
            return Err(Error::InvalidInput(
                "PED requires six pedigree columns".into(),
            ));
        }
        samples.push(sample_id(&prefix)?);
    }
    let chosen: usize = options.sample.choose(&samples)?;
    builder.source.sample = Some(samples[chosen].clone());
    let mut reader = containers::text_reader(path, options.max_bytes)?;
    let mut idx: usize = 0;
    while {
        line.clear();
        reader.read_line(&mut line)? != 0
    } {
        if line.trim().is_empty() {
            continue;
        }
        if idx == chosen {
            let row: Vec<_> = line.split_whitespace().collect();
            if row.len() != 6 + 2 * map.len() {
                return Err(Error::InvalidInput(
                    "PED allele count does not match MAP variants".into(),
                ));
            }
            for (variant, pair) in map.iter().zip(row[6..].as_chunks::<2>().0) {
                builder.push(locus(
                    &variant.id,
                    &variant.chromosome,
                    &variant.pos,
                    call(pair, CallCoding::Auto),
                ))?;
            }
            break;
        }
        idx += 1;
    }
    Ok(())
}

pub(super) fn read_bed(
    path: &Path,
    options: &ImportOptions,
    builder: &mut DatasetBuilder,
) -> Result<()> {
    let bim: Vec<Variant> = variants(&containers::companion(path, "bim")?, true, options)?;
    let samples: Vec<SampleId> = read_samples(&containers::companion(path, "fam")?, options)?;
    let chosen: usize = options.sample.choose(&samples)?;
    builder.source.sample = Some(samples[chosen].clone());
    let mut reader = containers::binary_reader(path, options.max_bytes)?;
    let mut header: [u8; 3] = [0; 3];
    reader.read_exact(&mut header)?;
    if header[..2] != [0x6c, 0x1b] || header[2] > 1 {
        return Err(Error::InvalidInput(
            "BED must be a PLINK binary file, not UCSC BED intervals".into(),
        ));
    }
    if header[2] == 1 {
        let mut block: Vec<u8> = vec![0; samples.len().div_ceil(4)];
        for variant in &bim {
            reader.read_exact(&mut block)?;
            bed_variant(
                variant,
                (block[chosen / 4] >> (2 * (chosen % 4))) & 3,
                builder,
            )?;
        }
    } else {
        let bytes: usize = bim.len().div_ceil(4);
        skip(
            &mut reader,
            chosen.checked_mul(bytes).ok_or(Error::Overflow)?,
        )?;
        let mut block: Vec<u8> = vec![0; bytes];
        reader.read_exact(&mut block)?;
        for (idx, variant) in bim.iter().enumerate() {
            bed_variant(variant, (block[idx / 4] >> (2 * (idx % 4))) & 3, builder)?;
        }
        skip(
            &mut reader,
            (samples.len() - chosen - 1)
                .checked_mul(bytes)
                .ok_or(Error::Overflow)?,
        )?;
    }
    if reader.read(&mut [0])? != 0 {
        return Err(Error::InvalidInput(
            "BED has trailing bytes inconsistent with BIM/FAM".into(),
        ));
    }
    Ok(())
}
