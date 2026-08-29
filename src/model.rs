use crate::{Error, Result};
use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, fmt, num::NonZeroU32, str::FromStr};

macro_rules! identifier {
    ($name:ident, $kind:literal, $valid:expr) => {
        #[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
        #[serde(try_from = "String", into = "String")]
        pub struct $name(Box<str>);
        impl $name {
            pub fn as_str(&self) -> &str {
                &self.0
            }
        }
        impl TryFrom<String> for $name {
            type Error = Error;
            fn try_from(value: String) -> Result<Self> {
                let text = value.trim();
                if text.is_empty() || !($valid)(text) {
                    return Err(Error::InvalidValue { kind: $kind, value });
                }
                Ok(Self(text.into()))
            }
        }
        impl FromStr for $name {
            type Err = Error;
            fn from_str(value: &str) -> Result<Self> {
                value.to_owned().try_into()
            }
        }
        impl From<$name> for String {
            fn from(value: $name) -> Self {
                value.0.into()
            }
        }
        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(self.as_str())
            }
        }
    };
}

identifier!(MarkerId, "marker identifier", |s: &str| !s
    .chars()
    .any(|c| c.is_whitespace() || c.is_control()));
identifier!(SampleId, "sample identifier", |s: &str| !s
    .chars()
    .any(char::is_control));

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(try_from = "u8", into = "u8")]
pub struct Autosome(u8);

impl Autosome {
    pub const fn get(self) -> u8 {
        self.0
    }
}
impl TryFrom<u8> for Autosome {
    type Error = Error;
    fn try_from(val: u8) -> Result<Self> {
        if (1..=22).contains(&val) {
            Ok(Self(val))
        } else {
            Err(Error::InvalidValue {
                kind: "autosome",
                value: val.to_string(),
            })
        }
    }
}
impl From<Autosome> for u8 {
    fn from(val: Autosome) -> Self {
        val.0
    }
}
impl FromStr for Autosome {
    type Err = Error;
    fn from_str(val: &str) -> Result<Self> {
        let txt = val.trim().trim_matches('"').to_ascii_uppercase();
        let txt = txt.strip_prefix("CHR").unwrap_or(&txt);
        let txt = txt
            .strip_prefix("NC_")
            .unwrap_or(txt)
            .split('.')
            .next()
            .unwrap_or(txt);
        txt.parse::<u8>()
            .ok()
            .and_then(|v| Self::try_from(v).ok())
            .ok_or_else(|| Error::InvalidValue {
                kind: "autosome",
                value: val.into(),
            })
    }
}
impl fmt::Display for Autosome {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(try_from = "u32", into = "u32")]
pub struct Position(NonZeroU32);

impl Position {
    pub const fn get(self) -> u32 {
        self.0.get()
    }
}
impl TryFrom<u32> for Position {
    type Error = Error;
    fn try_from(val: u32) -> Result<Self> {
        NonZeroU32::new(val)
            .map(Self)
            .ok_or_else(|| Error::InvalidValue {
                kind: "position",
                value: val.to_string(),
            })
    }
}
impl From<Position> for u32 {
    fn from(val: Position) -> Self {
        val.get()
    }
}
impl FromStr for Position {
    type Err = Error;
    fn from_str(val: &str) -> Result<Self> {
        val.trim()
            .parse::<u32>()
            .ok()
            .and_then(|v| v.try_into().ok())
            .ok_or_else(|| Error::InvalidValue {
                kind: "position",
                value: val.into(),
            })
    }
}
impl fmt::Display for Position {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.get().fmt(f)
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum Base {
    A,
    C,
    G,
    T,
}

impl Base {
    pub fn parse(val: &str) -> Option<Self> {
        match val.trim().to_ascii_uppercase().as_str() {
            "A" => Some(Self::A),
            "C" => Some(Self::C),
            "G" => Some(Self::G),
            "T" => Some(Self::T),
            _ => None,
        }
    }
}
impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum AssayAllele {
    A,
    B,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum NumericAllele {
    #[serde(rename = "1")]
    One,
    #[serde(rename = "2")]
    Two,
    #[serde(rename = "3")]
    Three,
    #[serde(rename = "4")]
    Four,
}
impl NumericAllele {
    pub fn parse(val: &str) -> Option<Self> {
        match val {
            "1" => Some(Self::One),
            "2" => Some(Self::Two),
            "3" => Some(Self::Three),
            "4" => Some(Self::Four),
            _ => None,
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Self::One => "1",
            Self::Two => "2",
            Self::Three => "3",
            Self::Four => "4",
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, ValueEnum)]
#[serde(rename_all = "snake_case")]
pub enum CallCoding {
    #[default]
    Auto,
    Dna,
    Assay,
    Numeric,
    Iupac,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(tag = "coding", content = "alleles", rename_all = "snake_case")]
pub enum Genotype {
    Dna([Base; 2]),
    Assay([AssayAllele; 2]),
    Numeric([NumericAllele; 2]),
}

impl Genotype {
    pub fn canonical(self) -> Self {
        match self {
            Self::Dna(mut pair) => {
                pair.sort_unstable();
                Self::Dna(pair)
            }
            Self::Assay(mut pair) => {
                pair.sort_unstable();
                Self::Assay(pair)
            }
            Self::Numeric(mut pair) => {
                pair.sort_unstable();
                Self::Numeric(pair)
            }
        }
    }
    pub fn is_homozygous(self) -> bool {
        match self {
            Self::Dna(pair) => pair[0] == pair[1],
            Self::Assay(pair) => pair[0] == pair[1],
            Self::Numeric(pair) => pair[0] == pair[1],
        }
    }
    pub fn plink_bases(self) -> [Base; 2] {
        match self {
            Self::Dna(pair) => pair,
            Self::Assay(pair) => pair.map(|a| {
                if a == AssayAllele::A {
                    Base::A
                } else {
                    Base::C
                }
            }),
            Self::Numeric(pair) => {
                if pair[0] == pair[1] {
                    [Base::A, Base::A]
                } else {
                    [Base::A, Base::C]
                }
            }
        }
    }
}
impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Dna(pair) => write!(f, "{}{}", pair[0], pair[1]),
            Self::Assay(pair) => write!(f, "{:?}{:?}", pair[0], pair[1]),
            Self::Numeric(pair) => write!(f, "{}/{}", pair[0].label(), pair[1].label()),
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct Snp {
    pub id: MarkerId,
    pub chromosome: Autosome,
    pub position: Position,
    pub genotype: Genotype,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, ValueEnum)]
#[serde(rename_all = "kebab-case")]
pub enum InputFormat {
    #[default]
    Auto,
    Table,
    Vcf,
    Ped,
    Tped,
    Bed,
    Bcf,
    Pgen,
    Bgen,
    Sam,
    Bam,
    Cram,
    Fastq,
    Fasta,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub enum Assembly {
    #[default]
    Unknown,
    GRCh36,
    GRCh37,
    GRCh38,
}

impl FromStr for Assembly {
    type Err = Error;
    fn from_str(val: &str) -> Result<Self> {
        match val.to_ascii_lowercase().as_str() {
            "unknown" | "auto" => Ok(Self::Unknown),
            "36" | "grch36" | "hg18" | "b36" => Ok(Self::GRCh36),
            "37" | "grch37" | "hg19" | "b37" => Ok(Self::GRCh37),
            "38" | "grch38" | "hg38" | "b38" => Ok(Self::GRCh38),
            _ => Err(Error::InvalidValue {
                kind: "genome assembly",
                value: val.into(),
            }),
        }
    }
}

impl Assembly {
    pub(crate) fn observe(&mut self, line: &str) -> Result<()> {
        let lbl: String = line.trim().to_ascii_lowercase();
        if lbl.starts_with("##")
            && !lbl.starts_with("##reference=")
            && !lbl.starts_with("##assembly=")
            && !(lbl.starts_with("##contig=") && lbl.contains("assembly="))
        {
            return Ok(());
        }
        if lbl.starts_with("#chrom\t") {
            return Ok(());
        }
        let txt: String = line
            .to_ascii_lowercase()
            .chars()
            .filter(|c| c.is_ascii_alphanumeric())
            .collect();
        let obs: Assembly = if txt.contains("grch38")
            || txt.contains("hg38")
            || txt.contains("build38")
            || txt.contains("assemblyb38")
        {
            Self::GRCh38
        } else if txt.contains("grch37")
            || txt.contains("hg19")
            || txt.contains("build37")
            || txt.contains("hs37d5")
            || txt.contains("assemblyb37")
        {
            Self::GRCh37
        } else if txt.contains("grch36")
            || txt.contains("hg18")
            || txt.contains("build36")
            || txt.contains("assemblyb36")
        {
            Self::GRCh36
        } else {
            return Ok(());
        };
        if *self != Self::Unknown && *self != obs {
            return Err(Error::InvalidInput(
                "conflicting genome assemblies; coordinates are not lifted over automatically"
                    .into(),
            ));
        }
        *self = obs;
        Ok(())
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct Source {
    pub format: InputFormat,
    pub assembly: Assembly,
    pub sample: Option<SampleId>,
    pub archive_member: Option<String>,
    pub conversions: Vec<String>,
    pub notes: Vec<String>,
}

impl Source {
    pub(crate) fn observe(&mut self, line: &str) -> Result<()> {
        if let Some(name) = line.trim().strip_prefix("# sample:") {
            let sample: SampleId = name.parse()?;
            if self.sample.as_ref().is_some_and(|cur| cur != &sample) {
                return Err(Error::InvalidInput("conflicting sample metadata".into()));
            }
            self.sample = Some(sample);
        } else {
            self.assembly.observe(line)?;
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Rejection {
    Missing,
    NonAutosomal,
    NonSnp,
    Haploid,
    Filtered,
    ReferenceBlock,
    Malformed,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct ImportStats {
    pub records_seen: u64,
    pub accepted_snps: u64,
    pub missing_calls: u64,
    pub non_autosomal: u64,
    pub non_snps: u64,
    pub haploid_or_polyploid: u64,
    pub filtered: u64,
    pub reference_blocks: u64,
    pub malformed: u64,
    pub duplicate_rows: u64,
    pub conflicting_loci: u64,
    pub encoded_calls: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct Dataset {
    source: Source,
    stats: ImportStats,
    snps: Vec<Snp>,
}

impl Dataset {
    pub fn snps(&self) -> &[Snp] {
        &self.snps
    }
    pub fn src(&self) -> &Source {
        &self.source
    }
    pub fn stats(&self) -> &ImportStats {
        &self.stats
    }
    pub fn from_snps(snps: impl IntoIterator<Item = Snp>) -> Result<Self> {
        let mut builder = DatasetBuilder::new(Source::default(), false);
        for snp in snps {
            builder.push(Ok(snp))?;
        }
        Ok(builder.finish())
    }
}

pub(crate) struct DatasetBuilder {
    pub source: Source,
    stats: ImportStats,
    loci: BTreeMap<(Autosome, Position), Option<Snp>>,
    strict: bool,
}

impl DatasetBuilder {
    pub fn new(source: Source, strict: bool) -> Self {
        Self {
            source,
            stats: ImportStats::default(),
            loci: BTreeMap::new(),
            strict,
        }
    }
    pub fn push(&mut self, snp: std::result::Result<Snp, Rejection>) -> Result<()> {
        self.stats.records_seen += 1;
        let mut snp: Snp = match snp {
            Ok(snp) => snp,
            Err(reason) => {
                let counter = match reason {
                    Rejection::Missing => &mut self.stats.missing_calls,
                    Rejection::NonAutosomal => &mut self.stats.non_autosomal,
                    Rejection::NonSnp => &mut self.stats.non_snps,
                    Rejection::Haploid => &mut self.stats.haploid_or_polyploid,
                    Rejection::Filtered => &mut self.stats.filtered,
                    Rejection::ReferenceBlock => &mut self.stats.reference_blocks,
                    Rejection::Malformed => &mut self.stats.malformed,
                };
                *counter += 1;
                if self.strict && reason == Rejection::Malformed {
                    return Err(Error::InvalidInput(format!(
                        "malformed genotype record {}",
                        self.stats.records_seen
                    )));
                }
                return Ok(());
            }
        };
        snp.genotype = snp.genotype.canonical();
        let key = (snp.chromosome, snp.position);
        match self.loci.entry(key) {
            std::collections::btree_map::Entry::Vacant(entry) => {
                entry.insert(Some(snp));
            }
            std::collections::btree_map::Entry::Occupied(mut entry) => {
                self.stats.duplicate_rows += 1;
                if entry
                    .get()
                    .as_ref()
                    .is_some_and(|old| old.genotype != snp.genotype)
                {
                    entry.insert(None);
                    self.stats.conflicting_loci += 1;
                }
            }
        }
        Ok(())
    }
    pub fn finish(mut self) -> Dataset {
        let snps: Vec<_> = self.loci.into_values().flatten().collect();
        self.stats.accepted_snps = snps.len() as u64;
        self.stats.encoded_calls = snps
            .iter()
            .filter(|s| !matches!(s.genotype, Genotype::Dna(_)))
            .count() as u64;
        if self.stats.encoded_calls > 0 {
            self.source.notes.push("Assay/numeric labels are preserved as labels; only their zygosity is used, and nucleotide identities are not inferred.".into());
        }
        Dataset {
            source: self.source,
            stats: self.stats,
            snps,
        }
    }
}
