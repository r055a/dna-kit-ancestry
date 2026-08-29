pub(crate) mod containers;
mod plink;
mod table;
mod vcf;

use crate::{
    Error, Result,
    model::*,
    units::{InputByteLimit, PhredScore, ReadDepth},
};
use std::{
    io::BufRead,
    num::NonZeroUsize,
    path::{Path, PathBuf},
};

pub use table::{ColumnMapping, Delimiter};

#[derive(Clone, Debug, Default)]
pub enum SampleSelection {
    #[default]
    Single,
    Name(SampleId),
    Index(NonZeroUsize),
}

impl SampleSelection {
    pub(crate) fn choose(&self, samples: &[SampleId]) -> Result<usize> {
        match self {
            Self::Single if samples.len() == 1 => Ok(0),
            Self::Name(name) => {
                let matches: Vec<_> = samples
                    .iter()
                    .enumerate()
                    .filter(|(_, id)| {
                        *id == name
                            || id
                                .as_str()
                                .split_once(':')
                                .is_some_and(|(_, iid)| iid == name.as_str())
                    })
                    .map(|(i, _)| i)
                    .collect();
                if matches.len() == 1 {
                    Ok(matches[0])
                } else {
                    Err(sample_error(samples))
                }
            }
            Self::Index(idx) if idx.get() <= samples.len() => Ok(idx.get() - 1),
            _ => Err(sample_error(samples)),
        }
    }
}

fn sample_error(samples: &[SampleId]) -> Error {
    if samples.is_empty() {
        return Error::InvalidInput("no individuals are declared in the sample metadata".into());
    }
    Error::Ambiguous {
        kind: "sample",
        option: "--sample or --sample-index",
        choices: samples
            .iter()
            .take(20)
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(", "),
    }
}

#[derive(Clone, Debug)]
pub struct ImportOptions {
    pub format: InputFormat,
    pub assembly: Assembly,
    pub sample: SampleSelection,
    pub member: Option<String>,
    pub cols: Option<ColumnMapping>,
    pub delimiter: Delimiter,
    pub coding: Option<CallCoding>,
    pub strict: bool,
    pub include_filtered: bool,
    pub min_depth: Option<ReadDepth>,
    pub min_genotype_quality: Option<PhredScore>,
    pub max_bytes: InputByteLimit,
}

impl Default for ImportOptions {
    fn default() -> Self {
        Self {
            format: InputFormat::Auto,
            assembly: Assembly::Unknown,
            sample: SampleSelection::Single,
            member: None,
            cols: None,
            delimiter: Delimiter::Auto,
            coding: None,
            strict: false,
            include_filtered: false,
            min_depth: None,
            min_genotype_quality: None,
            max_bytes: InputByteLimit::try_from(16 * 1024 * 1024 * 1024).expect("nonzero constant"),
        }
    }
}

#[derive(Clone, Debug)]
pub struct ToolOptions {
    pub plink: PathBuf,
    pub plink2: PathBuf,
    pub bcftools: PathBuf,
    pub samtools: PathBuf,
    pub bwa: PathBuf,
    pub reference: Option<PathBuf>,
    pub sites: Option<PathBuf>,
    pub mate: Option<PathBuf>,
    pub threads: NonZeroUsize,
}

impl Default for ToolOptions {
    fn default() -> Self {
        Self {
            plink: crate::install::default_command(),
            plink2: "plink2".into(),
            bcftools: "bcftools".into(),
            samtools: "samtools".into(),
            bwa: "bwa".into(),
            reference: None,
            sites: None,
            mate: None,
            threads: NonZeroUsize::new(1).expect("nonzero constant"),
        }
    }
}

pub fn load(path: &Path, options: &ImportOptions, tools: &ToolOptions) -> Result<Dataset> {
    let prepared = containers::prepare(path, options)?;
    let mut reader = containers::text_reader(&prepared.path, options.max_bytes)?;
    let format: InputFormat = if options.format == InputFormat::Auto {
        detect(&prepared.path, reader.fill_buf()?)?
    } else {
        options.format
    };
    if format != InputFormat::Table
        && (options.cols.is_some()
            || options.coding.is_some()
            || options.delimiter != Delimiter::Auto)
    {
        return Err(Error::InvalidInput(
            "--columns, --coding, and --delimiter apply only to genotype tables".into(),
        ));
    }
    if matches!(
        format,
        InputFormat::Table | InputFormat::Ped | InputFormat::Tped | InputFormat::Bed
    ) && (options.min_depth.is_some() || options.min_genotype_quality.is_some())
    {
        return Err(Error::InvalidInput(
            "--min-depth and --min-gq require VCF FORMAT/DP and FORMAT/GQ fields".into(),
        ));
    }
    let src = Source {
        format,
        assembly: options.assembly,
        archive_member: prepared.member.clone(),
        ..Source::default()
    };
    let mut builder: DatasetBuilder = DatasetBuilder::new(src, options.strict);
    match format {
        InputFormat::Table => table::read(reader, options, &mut builder)?,
        InputFormat::Vcf => vcf::read(reader, options, &mut builder)?,
        InputFormat::Ped => plink::read_ped(&prepared.path, options, &mut builder)?,
        InputFormat::Tped => plink::read_tped(&prepared.path, options, &mut builder)?,
        InputFormat::Bed => plink::read_bed(&prepared.path, options, &mut builder)?,
        InputFormat::Bcf
        | InputFormat::Pgen
        | InputFormat::Bgen
        | InputFormat::Sam
        | InputFormat::Bam
        | InputFormat::Cram
        | InputFormat::Fastq => {
            let converted = crate::process::to_vcf(&prepared.path, format, options, tools)?;
            builder.source.conversions = converted.steps.clone();
            vcf::read(
                containers::text_reader(&converted.path, options.max_bytes)?,
                options,
                &mut builder,
            )?;
        }
        InputFormat::Fasta => {
            return Err(Error::NeedsConversion {
                format: "FASTA",
                requirement: "diploid variant calls; a consensus or haploid sequence cannot establish autosomal heterozygosity",
            });
        }
        InputFormat::Auto => {
            return Err(Error::InvalidInput("format could not be determined".into()));
        }
    }
    Ok(builder.finish())
}

fn detect(path: &Path, prefix: &[u8]) -> Result<InputFormat> {
    if prefix.starts_with(&[0x6c, 0x1b]) {
        return Ok(if matches!(prefix.get(2), Some(0 | 1)) {
            InputFormat::Bed
        } else {
            InputFormat::Pgen
        });
    }
    if prefix.starts_with(b"BCF") {
        return Ok(InputFormat::Bcf);
    }
    if prefix.starts_with(b"BAM\x01") {
        return Ok(InputFormat::Bam);
    }
    if prefix.starts_with(b"CRAM") {
        return Ok(InputFormat::Cram);
    }
    let name: String = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_lowercase();
    let name: &str = containers::uncompressed_name(&name);
    let extension: &str = name.rsplit('.').next().unwrap_or_default();
    let explicit: Option<InputFormat> = match extension {
        "ped" => Some(InputFormat::Ped),
        "tped" => Some(InputFormat::Tped),
        "bed" => Some(InputFormat::Bed),
        "pgen" => Some(InputFormat::Pgen),
        "bgen" => Some(InputFormat::Bgen),
        "bcf" => Some(InputFormat::Bcf),
        "bam" => Some(InputFormat::Bam),
        "cram" => Some(InputFormat::Cram),
        "sam" => Some(InputFormat::Sam),
        "fastq" | "fq" => Some(InputFormat::Fastq),
        "fasta" | "fa" | "fna" => Some(InputFormat::Fasta),
        "vcf" | "gvcf" => Some(InputFormat::Vcf),
        _ => None,
    };
    if let Some(format) = explicit {
        return Ok(format);
    }
    if matches!(
        extension,
        "tbi"
            | "csi"
            | "bai"
            | "crai"
            | "gzi"
            | "fai"
            | "map"
            | "bim"
            | "fam"
            | "tfam"
            | "pvar"
            | "psam"
            | "sample"
    ) {
        return Err(Error::NeedsConversion {
            format: "index/metadata sidecar",
            requirement: "the associated primary VCF, BAM/CRAM, PED/TPED, BED or PGEN file",
        });
    }
    if prefix.starts_with(b"7z\xbc\xaf\x27\x1c") || prefix.starts_with(b"Rar!") {
        return Err(Error::NeedsConversion {
            format: "7z/RAR archive",
            requirement: "local extraction to a supported raw data file",
        });
    }
    let text = String::from_utf8_lossy(prefix);
    if text.starts_with("##fileformat=VCF") || text.starts_with("#CHROM\tPOS\t") {
        return Ok(InputFormat::Vcf);
    }
    if text.starts_with("@HD\t") || text.starts_with("@SQ\t") {
        return Ok(InputFormat::Sam);
    }
    if text.starts_with('>') {
        return Ok(InputFormat::Fasta);
    }
    if text.starts_with('@')
        && text
            .lines()
            .nth(2)
            .is_some_and(|line| line.starts_with('+'))
    {
        return Ok(InputFormat::Fastq);
    }
    if prefix.starts_with(b"%PDF")
        || prefix.starts_with(b"\x89PNG")
        || prefix.starts_with(b"SQLite")
    {
        return Err(Error::NeedsConversion {
            format: "report/image/database",
            requirement: "a raw genotype table or VCF export with chromosome, position, and diploid calls",
        });
    }
    Ok(InputFormat::Table)
}

pub(crate) fn locus(
    id: &str,
    chromosome: &str,
    position: &str,
    genotype: std::result::Result<Genotype, Rejection>,
) -> std::result::Result<Snp, Rejection> {
    let chromosome: Autosome = chromosome
        .parse::<Autosome>()
        .map_err(|_| Rejection::NonAutosomal)?;
    let position: Position = position
        .parse::<Position>()
        .map_err(|_| Rejection::Malformed)?;
    let genotype: Genotype = genotype?;
    let id: MarkerId = id.parse().unwrap_or_else(|_| {
        format!("chr{chromosome}:{position}")
            .parse()
            .expect("generated valid identifier")
    });
    Ok(Snp {
        id,
        chromosome,
        position,
        genotype,
    })
}

pub(crate) fn call(
    alleles: &[&str],
    coding: CallCoding,
) -> std::result::Result<Genotype, Rejection> {
    if alleles.iter().any(|a| {
        matches!(
            a.trim().to_ascii_uppercase().as_str(),
            "" | "." | "0" | "-" | "--" | "00" | "NN" | "N" | "?" | "NC" | "NOCALL" | "NO_CALL"
        )
    }) {
        return Err(Rejection::Missing);
    }
    if alleles.len() != 2 {
        return Err(Rejection::Haploid);
    }
    if coding == CallCoding::Assay {
        let parse = |s: &str| match s.trim().to_ascii_uppercase().as_str() {
            "A" => Ok(AssayAllele::A),
            "B" => Ok(AssayAllele::B),
            _ => Err(Rejection::Malformed),
        };
        Ok(Genotype::Assay([parse(alleles[0])?, parse(alleles[1])?]))
    } else {
        if matches!(coding, CallCoding::Auto | CallCoding::Numeric)
            && let (Some(a), Some(b)) = (
                NumericAllele::parse(alleles[0]),
                NumericAllele::parse(alleles[1]),
            )
        {
            return Ok(Genotype::Numeric([a, b]));
        }
        if coding == CallCoding::Numeric {
            return Err(Rejection::Malformed);
        }
        let a: Base = Base::parse(alleles[0]).ok_or(Rejection::NonSnp)?;
        let b: Base = Base::parse(alleles[1]).ok_or(Rejection::NonSnp)?;
        Ok(Genotype::Dna([a, b]))
    }
}

pub(crate) fn combined_call(
    val: &str,
    coding: CallCoding,
) -> std::result::Result<Genotype, Rejection> {
    let val: String = val.trim().to_ascii_uppercase();
    if val.is_empty() {
        return Err(Rejection::Missing);
    }
    if val.contains(['/', '|', ' ']) {
        return call(
            &val.split(['/', '|', ' '])
                .filter(|s| !s.is_empty())
                .collect::<Vec<_>>(),
            coding,
        );
    }
    if val.len() == 1 && coding == CallCoding::Iupac {
        let iupac = match val.as_str() {
            "A" => Some([Base::A; 2]),
            "C" => Some([Base::C; 2]),
            "G" => Some([Base::G; 2]),
            "T" => Some([Base::T; 2]),
            "R" => Some([Base::A, Base::G]),
            "Y" => Some([Base::C, Base::T]),
            "S" => Some([Base::C, Base::G]),
            "W" => Some([Base::A, Base::T]),
            "K" => Some([Base::G, Base::T]),
            "M" => Some([Base::A, Base::C]),
            _ => None,
        };
        if let Some(pair) = iupac {
            return Ok(Genotype::Dna(pair));
        }
    }
    let chars: Vec<String> = val.chars().map(|c| c.to_string()).collect();
    call(
        &chars.iter().map(String::as_str).collect::<Vec<_>>(),
        coding,
    )
}
