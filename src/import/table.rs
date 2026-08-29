use super::{ImportOptions, SampleSelection, combined_call, locus};
use crate::{Error, Result, model::*};
use clap::ValueEnum;
use std::{io::BufRead, str::FromStr};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
pub enum Delimiter {
    #[default]
    Auto,
    Tab,
    Comma,
    Semicolon,
    Whitespace,
}

impl Delimiter {
    fn detect(line: &str) -> Self {
        if line.contains('\t') {
            Self::Tab
        } else if line.contains(',') {
            Self::Comma
        } else if line.contains(';') {
            Self::Semicolon
        } else {
            Self::Whitespace
        }
    }
    fn byte(self) -> u8 {
        match self {
            Self::Tab => b'\t',
            Self::Comma => b',',
            Self::Semicolon => b';',
            _ => b' ',
        }
    }
}

#[derive(Clone, Debug)]
enum Calls {
    Combined {
        idx: usize,
        coding: CallCoding,
    },
    Split {
        first: usize,
        second: usize,
        coding: CallCoding,
    },
}

#[derive(Clone, Debug)]
pub struct ColumnMapping {
    id: Option<usize>,
    chromosome: usize,
    pos: usize,
    calls: Calls,
    sample: Option<usize>,
    coding: Option<usize>,
}

fn norm(val: &str) -> String {
    val.chars()
        .filter(|c| c.is_alphanumeric())
        .flat_map(char::to_lowercase)
        .collect()
}

impl ColumnMapping {
    fn header(fields: &[String]) -> Option<Self> {
        let fields: Vec<_> = fields.iter().map(|s| norm(s)).collect();
        let find = |aliases: &[&str]| {
            aliases
                .iter()
                .find_map(|alias| fields.iter().position(|name| name == alias))
        };
        let chromosome: usize = find(&[
            "chromosome",
            "chrom",
            "chr",
            "chromossome",
            "cromossomo",
            "cromosoma",
            "chromosom",
        ])?;
        let pos: usize = find(&[
            "position",
            "pos",
            "bp",
            "basepairposition",
            "physicalposition",
            "positionbp",
            "positiongrch37",
            "positiongrch38",
            "positionhg19",
            "positionhg38",
            "mapinfo",
            "posição",
            "posicao",
            "posicion",
            "posición",
        ])?;
        let first: Option<usize> = find(&[
            "allele1forward",
            "allele1plus",
            "allele1forwardstrand",
            "allele1",
            "alleleone",
            "allelea",
            "a1",
        ]);
        let second: Option<usize> = find(&[
            "allele2forward",
            "allele2plus",
            "allele2forwardstrand",
            "allele2",
            "alleletwo",
            "alleleb",
            "a2",
        ]);
        let calls: Calls = if let (Some(first), Some(second)) = (first, second) {
            Calls::Split {
                first,
                second,
                coding: CallCoding::Auto,
            }
        } else if let (Some(first), Some(second)) = (find(&["allele1ab"]), find(&["allele2ab"])) {
            Calls::Split {
                first,
                second,
                coding: CallCoding::Assay,
            }
        } else if let Some(idx) = find(&[
            "genotype",
            "result",
            "gt",
            "call",
            "alleles",
            "genotypeforward",
            "genotypeplus",
            "genotipo",
            "genótipo",
        ]) {
            Calls::Combined {
                idx,
                coding: CallCoding::Auto,
            }
        } else if let Some(idx) = find(&["iupac", "genotypeiupac"]) {
            Calls::Combined {
                idx,
                coding: CallCoding::Iupac,
            }
        } else {
            Calls::Combined {
                idx: find(&["abgenotype", "genotypeab"])?,
                coding: CallCoding::Assay,
            }
        };
        Some(Self {
            id: find(&[
                "rsid",
                "snp",
                "snpid",
                "snpname",
                "marker",
                "markername",
                "id",
                "probesetid",
            ]),
            chromosome,
            pos,
            calls,
            sample: find(&["sampleid", "sample", "individual", "individualid"]),
            coding: find(&["coding", "encoding", "allelecoding"]),
        })
    }

    fn infer(fields: &[String]) -> Option<Self> {
        if !matches!(fields.len(), 4 | 5)
            || !fields[0].chars().any(char::is_alphabetic)
            || fields[1].parse::<Autosome>().is_err()
            || fields[2].parse::<Position>().is_err()
        {
            return None;
        }
        Some(Self {
            id: Some(0),
            chromosome: 1,
            pos: 2,
            sample: None,
            coding: None,
            calls: if fields.len() == 4 {
                Calls::Combined {
                    idx: 3,
                    coding: CallCoding::Auto,
                }
            } else {
                Calls::Split {
                    first: 3,
                    second: 4,
                    coding: CallCoding::Auto,
                }
            },
        })
    }

    fn max_index(&self) -> usize {
        let allele: usize = match self.calls {
            Calls::Combined { idx, .. } => idx,
            Calls::Split { first, second, .. } => first.max(second),
        };
        allele
            .max(self.chromosome)
            .max(self.pos)
            .max(self.id.unwrap_or(0))
            .max(self.sample.unwrap_or(0))
            .max(self.coding.unwrap_or(0))
    }

    fn snp(
        &self,
        fields: &[String],
        override_coding: Option<CallCoding>,
    ) -> std::result::Result<Snp, Rejection> {
        if fields.len() <= self.max_index() {
            return Err(Rejection::Malformed);
        }
        let coding: Option<CallCoding> = match self
            .coding
            .map(|idx| fields[idx].trim().to_ascii_lowercase())
        {
            None => None,
            Some(val) if matches!(val.as_str(), "dna" | "nucleotide" | "bases") => {
                Some(CallCoding::Dna)
            }
            Some(val) if val == "numeric" => Some(CallCoding::Numeric),
            Some(val) if val == "iupac" => Some(CallCoding::Iupac),
            Some(val) if matches!(val.as_str(), "assay" | "ab") => Some(CallCoding::Assay),
            Some(_) => return Err(Rejection::Malformed),
        };
        let genotype = match self.calls {
            Calls::Combined {
                idx,
                coding: default,
            } => combined_call(&fields[idx], override_coding.or(coding).unwrap_or(default)),
            Calls::Split {
                first,
                second,
                coding: default,
            } => super::call(
                &[&fields[first], &fields[second]],
                override_coding.or(coding).unwrap_or(default),
            ),
        };
        locus(
            self.id.map(|idx| fields[idx].as_str()).unwrap_or(""),
            &fields[self.chromosome],
            &fields[self.pos],
            genotype,
        )
    }
}

impl FromStr for ColumnMapping {
    type Err = Error;
    fn from_str(val: &str) -> Result<Self> {
        let mut cols = std::collections::BTreeMap::new();
        for assignment in val.split(',') {
            let (name, idx) = assignment
                .split_once('=')
                .ok_or_else(|| Error::InvalidInput("columns use name=one-based-index".into()))?;
            let name = name.trim();
            if ![
                "id", "chr", "pos", "gt", "a1", "a2", "sample", "ab", "coding",
            ]
            .contains(&name)
            {
                return Err(Error::InvalidInput(format!("unknown column role {name:?}")));
            }
            let idx = idx
                .trim()
                .parse::<usize>()
                .ok()
                .and_then(|i| i.checked_sub(1))
                .ok_or_else(|| {
                    Error::InvalidInput("column indices must be positive integers".into())
                })?;
            if cols.insert(name, idx).is_some() {
                return Err(Error::InvalidInput("duplicate column role".into()));
            }
        }
        let req = |name| {
            cols.get(name)
                .copied()
                .ok_or_else(|| Error::InvalidInput(format!("missing {name} column")))
        };
        let calls: Calls = match (
            cols.get("gt"),
            cols.get("ab"),
            cols.get("a1"),
            cols.get("a2"),
        ) {
            (Some(idx), None, None, None) => Calls::Combined {
                idx: *idx,
                coding: CallCoding::Auto,
            },
            (None, Some(idx), None, None) => Calls::Combined {
                idx: *idx,
                coding: CallCoding::Assay,
            },
            (None, None, Some(first), Some(second)) => Calls::Split {
                first: *first,
                second: *second,
                coding: CallCoding::Auto,
            },
            _ => {
                return Err(Error::InvalidInput(
                    "choose gt, ab, or both a1 and a2 columns".into(),
                ));
            }
        };
        let unique: std::collections::BTreeSet<_> = cols.values().collect();
        if unique.len() != cols.len() {
            return Err(Error::InvalidInput(
                "column roles must use distinct indices".into(),
            ));
        }
        Ok(Self {
            id: cols.get("id").copied(),
            chromosome: req("chr")?,
            pos: req("pos")?,
            calls,
            sample: cols.get("sample").copied(),
            coding: cols.get("coding").copied(),
        })
    }
}

fn fields(line: &str, delimiter: Delimiter) -> Result<Vec<String>> {
    if delimiter == Delimiter::Whitespace {
        return Ok(line.split_whitespace().map(str::to_owned).collect());
    }
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .flexible(true)
        .delimiter(delimiter.byte())
        .from_reader(line.as_bytes());
    Ok(reader
        .records()
        .next()
        .transpose()?
        .map(|row| row.iter().map(str::to_owned).collect())
        .unwrap_or_default())
}

struct TableSamples<'a> {
    selection: &'a SampleSelection,
    names: Vec<SampleId>,
}
impl TableSamples<'_> {
    fn accept(&mut self, name: &str) -> Result<bool> {
        let id: SampleId = name.parse()?;
        let idx = if let Some(idx) = self.names.iter().position(|existing| existing == &id) {
            idx
        } else {
            self.names.push(id.clone());
            self.names.len() - 1
        };
        match self.selection {
            SampleSelection::Single if self.names.len() > 1 => {
                Err(super::sample_error(&self.names))
            }
            SampleSelection::Single => Ok(true),
            SampleSelection::Name(req) => Ok(req == &id),
            SampleSelection::Index(req) => Ok(req.get() - 1 == idx),
        }
    }
}

pub(super) fn read(
    mut reader: Box<dyn BufRead>,
    options: &ImportOptions,
    builder: &mut DatasetBuilder,
) -> Result<()> {
    let mut line: String = String::new();
    let mut first = None;
    let mut mapping = None;
    let mut delimiter: Delimiter = options.delimiter;
    while reader.read_line(&mut line)? != 0 {
        let txt: &str = line.trim().trim_start_matches('\u{feff}');
        if txt.is_empty() {
            line.clear();
            continue;
        }
        let unmarked: &str = txt.trim_start_matches('#').trim();
        let selected: Delimiter = if options.delimiter == Delimiter::Auto {
            Delimiter::detect(unmarked)
        } else {
            options.delimiter
        };
        let row: Vec<String> = fields(unmarked, selected)?;
        if let Some(header) = ColumnMapping::header(&row) {
            builder.source.observe(&line)?;
            mapping = Some(options.cols.clone().unwrap_or(header));
            delimiter = selected;
            break;
        }
        if !txt.starts_with('#')
            && !txt.starts_with('[')
            && let Some(schema) = options.cols.clone().or_else(|| ColumnMapping::infer(&row))
            && row.len() > schema.max_index()
        {
            mapping = Some(schema);
            delimiter = selected;
            first = Some(row);
            break;
        }
        builder.source.observe(&line)?;
        line.clear();
    }
    let mapping: ColumnMapping = mapping.ok_or_else(|| Error::InvalidInput("no genotype schema found; use --columns id=1,chr=2,pos=3,gt=4 and --delimiter for an unlabelled layout".into()))?;
    let mut samples = TableSamples {
        selection: &options.sample,
        names: Vec::new(),
    };
    let mut ingest = |row: Vec<String>| -> Result<()> {
        if row.is_empty() {
            return Ok(());
        }
        if let Some(col) = mapping.sample {
            if row.len() <= col {
                return builder.push(Err(Rejection::Malformed));
            }
            if !samples.accept(&row[col])? {
                return Ok(());
            }
        }
        builder.push(mapping.snp(&row, options.coding))
    };
    if let Some(row) = first {
        ingest(row)?;
    }
    if delimiter == Delimiter::Whitespace {
        while {
            line.clear();
            reader.read_line(&mut line)? != 0
        } {
            if !line.trim().is_empty() && !line.trim_start().starts_with('#') {
                ingest(fields(&line, delimiter)?)?;
            }
        }
    } else {
        let mut records = csv::ReaderBuilder::new()
            .has_headers(false)
            .flexible(true)
            .comment(Some(b'#'))
            .delimiter(delimiter.byte())
            .from_reader(reader);
        for record in records.records() {
            ingest(record?.iter().map(str::to_owned).collect())?;
        }
    }
    if mapping.sample.is_some() {
        let selected: usize = options.sample.choose(&samples.names)?;
        builder.source.sample = Some(samples.names[selected].clone());
    } else {
        match options.sample {
            SampleSelection::Single => {}
            SampleSelection::Index(idx) if idx.get() == 1 => {}
            SampleSelection::Name(ref name) if builder.source.sample.as_ref() == Some(name) => {}
            _ => {
                return Err(Error::InvalidInput(
                    "this table has no sample-ID column; omit --sample".into(),
                ));
            }
        }
    }
    Ok(())
}
