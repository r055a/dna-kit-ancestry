use crate::{
    Error, Result,
    import::ToolOptions,
    model::*,
    units::{BasePairCount, BasePairLength, SnpDensity},
};
use clap::ValueEnum;
use serde::Serialize;
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufRead, BufReader, Write},
    num::NonZeroUsize,
    process::Command,
};

#[derive(Clone, Copy, Debug, Default, Serialize, ValueEnum)]
#[serde(rename_all = "snake_case")]
pub enum Method {
    #[default]
    Plink,
    Native,
}

#[derive(Clone, Debug, Serialize)]
pub struct RohConfig {
    pub min_snps: NonZeroUsize,
    pub min_length_bp: BasePairLength,
    pub max_gap_bp: BasePairLength,
    pub max_heterozygous: usize,
    pub max_density_bp_per_snp: SnpDensity,
    pub window_snps: NonZeroUsize,
}

impl Default for RohConfig {
    fn default() -> Self {
        Self {
            min_snps: NonZeroUsize::new(100).expect("nonzero constant"),
            min_length_bp: BasePairLength::try_from(1_000_000).expect("nonzero constant"),
            max_gap_bp: BasePairLength::try_from(1_000_000).expect("nonzero constant"),
            max_heterozygous: 1,
            max_density_bp_per_snp: SnpDensity::try_from(50_000).expect("nonzero constant"),
            window_snps: NonZeroUsize::new(50).expect("nonzero constant"),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Serialize)]
#[serde(tag = "kind", content = "length_bp", rename_all = "snake_case")]
pub enum Denominator {
    #[default]
    ObservedSpan,
    Provided(BasePairLength),
}

#[derive(Clone, Debug, Serialize)]
pub struct RohSegment {
    pub chromosome: Autosome,
    start: Position,
    end: Position,
    pub marker_count: Option<NonZeroUsize>,
}

impl RohSegment {
    pub fn new(
        chromosome: Autosome,
        start: Position,
        end: Position,
        marker_count: Option<NonZeroUsize>,
    ) -> Result<Self> {
        if end < start {
            return Err(Error::InvalidInput(
                "ROH endpoint precedes its start".into(),
            ));
        }
        Ok(Self {
            chromosome,
            start,
            end,
            marker_count,
        })
    }
    pub fn start(&self) -> Position {
        self.start
    }
    pub fn end(&self) -> Position {
        self.end
    }
    pub fn length_bp(&self) -> BasePairLength {
        BasePairLength::try_from(u64::from(self.end.get() - self.start.get()) + 1)
            .expect("validated inclusive interval")
    }
}

#[derive(Clone, Copy, Debug, Serialize)]
#[serde(transparent)]
pub struct Fraction(f64);
impl Fraction {
    pub fn value(self) -> f64 {
        self.0
    }
    pub fn percent(self) -> f64 {
        self.0 * 100.0
    }
}

#[derive(Debug, Serialize)]
pub struct Analysis {
    src: Source,
    import: ImportStats,
    method: Method,
    config: RohConfig,
    homozygous_snps: usize,
    heterozygous_snps: usize,
    denominator: Denominator,
    denominator_bp: BasePairLength,
    total_roh_bp: BasePairCount,
    f_roh: Fraction,
    segments: Vec<RohSegment>,
}

impl Analysis {
    pub fn src(&self) -> &Source {
        &self.src
    }
    pub fn import(&self) -> &ImportStats {
        &self.import
    }
    pub fn method(&self) -> Method {
        self.method
    }
    pub fn config(&self) -> &RohConfig {
        &self.config
    }
    pub fn homozygous_snps(&self) -> usize {
        self.homozygous_snps
    }
    pub fn heterozygous_snps(&self) -> usize {
        self.heterozygous_snps
    }
    pub fn denominator(&self) -> Denominator {
        self.denominator
    }
    pub fn denominator_bp(&self) -> BasePairLength {
        self.denominator_bp
    }
    pub fn total_roh_bp(&self) -> BasePairCount {
        self.total_roh_bp
    }
    pub fn f_roh(&self) -> Fraction {
        self.f_roh
    }
    pub fn segments(&self) -> &[RohSegment] {
        &self.segments
    }
}

pub fn analyze(
    dataset: &Dataset,
    config: &RohConfig,
    method: Method,
    denominator: Denominator,
    tools: &ToolOptions,
) -> Result<Analysis> {
    if dataset.snps().is_empty() {
        return Err(Error::NoUsableCalls);
    }
    let segments: Vec<RohSegment> = match method {
        Method::Native => native(dataset, config)?,
        Method::Plink => plink(dataset, config, tools)?,
    };
    let segments: Vec<RohSegment> = union(segments);
    let total_roh_bp: u64 = segments.iter().try_fold(0_u64, |sum, segment| {
        sum.checked_add(segment.length_bp().get())
            .ok_or(Error::Overflow)
    })?;
    let denominator_bp: BasePairLength = match denominator {
        Denominator::ObservedSpan => observed_span(dataset)?,
        Denominator::Provided(val) => val,
    };
    if total_roh_bp > denominator_bp.get() {
        return Err(Error::InvalidInput(
            "ROH length exceeds the denominator; check coordinates and denominator".into(),
        ));
    }
    let homozygous_snps: usize = dataset
        .snps()
        .iter()
        .filter(|snp| snp.genotype.is_homozygous())
        .count();
    Ok(Analysis {
        src: dataset.src().clone(),
        import: dataset.stats().clone(),
        method,
        config: config.clone(),
        homozygous_snps,
        heterozygous_snps: dataset.snps().len() - homozygous_snps,
        denominator,
        denominator_bp,
        total_roh_bp: total_roh_bp.into(),
        f_roh: Fraction(total_roh_bp as f64 / denominator_bp.get() as f64),
        segments,
    })
}

fn observed_span(dataset: &Dataset) -> Result<BasePairLength> {
    let mut ranges = BTreeMap::<Autosome, (Position, Position)>::new();
    for snp in dataset.snps() {
        ranges
            .entry(snp.chromosome)
            .and_modify(|range| range.1 = snp.position)
            .or_insert((snp.position, snp.position));
    }
    let total: u64 = ranges.values().try_fold(0_u64, |sum, (start, end)| {
        sum.checked_add(u64::from(end.get() - start.get()) + 1)
            .ok_or(Error::Overflow)
    })?;
    total.try_into()
}

struct Candidate {
    chromosome: Autosome,
    start: Position,
    end: Position,
    count: usize,
    heterozygous: usize,
}
impl Candidate {
    fn finish(self, config: &RohConfig, segments: &mut Vec<RohSegment>) -> Result<()> {
        let segment: RohSegment = RohSegment::new(
            self.chromosome,
            self.start,
            self.end,
            NonZeroUsize::new(self.count),
        )?;
        if self.count >= config.min_snps.get()
            && segment.length_bp() >= config.min_length_bp
            && u128::from(segment.length_bp().get())
                <= u128::from(config.max_density_bp_per_snp.get()) * self.count as u128
        {
            segments.push(segment);
        }
        Ok(())
    }
}

fn native(dataset: &Dataset, config: &RohConfig) -> Result<Vec<RohSegment>> {
    let mut candidate: Option<Candidate> = None;
    let mut segments: Vec<RohSegment> = Vec::new();
    for snp in dataset.snps() {
        let is_boundary: bool = candidate.as_ref().is_some_and(|run| {
            run.chromosome != snp.chromosome
                || u64::from(snp.position.get().saturating_sub(run.end.get()))
                    > config.max_gap_bp.get()
                || (!snp.genotype.is_homozygous() && run.heterozygous >= config.max_heterozygous)
        });
        if is_boundary && let Some(run) = candidate.take() {
            run.finish(config, &mut segments)?;
        }
        match &mut candidate {
            Some(run) => {
                run.end = snp.position;
                run.count += 1;
                run.heterozygous += usize::from(!snp.genotype.is_homozygous());
            }
            None if snp.genotype.is_homozygous() => {
                candidate = Some(Candidate {
                    chromosome: snp.chromosome,
                    start: snp.position,
                    end: snp.position,
                    count: 1,
                    heterozygous: 0,
                });
            }
            None => {}
        }
    }
    if let Some(run) = candidate {
        run.finish(config, &mut segments)?;
    }
    Ok(segments)
}

fn union(mut segments: Vec<RohSegment>) -> Vec<RohSegment> {
    segments.sort_by_key(|segment| (segment.chromosome, segment.start));
    let mut merged: Vec<RohSegment> = Vec::new();
    for segment in segments {
        if let Some(last) = merged
            .last_mut()
            .filter(|last| last.chromosome == segment.chromosome && segment.start <= last.end)
        {
            last.end = last.end.max(segment.end);
            last.marker_count = None;
        } else {
            merged.push(segment);
        }
    }
    merged
}

fn kb(bp: u64) -> String {
    format!("{}.{:03}", bp / 1000, bp % 1000)
}

fn plink(dataset: &Dataset, config: &RohConfig, tools: &ToolOptions) -> Result<Vec<RohSegment>> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("genotypes");
    let output = dir.path().join("roh");
    let mut bed: File = File::create(input.with_extension("bed"))?;
    let mut bim = std::io::BufWriter::new(File::create(input.with_extension("bim"))?);
    bed.write_all(&[0x6c, 0x1b, 1])?;
    for snp in dataset.snps() {
        let [a, b] = snp.genotype.plink_bases();
        writeln!(
            bim,
            "{}\tchr{}:{}\t0\t{}\t{}\t{}",
            snp.chromosome,
            snp.chromosome,
            snp.position,
            snp.position,
            a,
            if a == b {
                "0".to_owned()
            } else {
                b.to_string()
            }
        )?;
        bed.write_all(&[if a == b { 0 } else { 2 }])?;
    }
    bim.flush()?;
    drop(bim);
    drop(bed);
    std::fs::write(input.with_extension("fam"), "DNA SAMPLE 0 0 0 -9\n")?;
    let mut command = Command::new(&tools.plink);
    command
        .arg("--bfile")
        .arg(&input)
        .arg("--threads")
        .arg(tools.threads.to_string())
        .args(["--homozyg", "--homozyg-snp"])
        .arg(config.min_snps.to_string())
        .arg("--homozyg-kb")
        .arg(kb(config.min_length_bp.get()))
        .arg("--homozyg-density")
        .arg(kb(config.max_density_bp_per_snp.get()))
        .arg("--homozyg-gap")
        .arg(kb(config.max_gap_bp.get()))
        .arg("--homozyg-het")
        .arg(config.max_heterozygous.to_string())
        .arg("--homozyg-window-het")
        .arg(config.max_heterozygous.to_string())
        .arg("--homozyg-window-snp")
        .arg(config.window_snps.to_string())
        .args(["--allow-no-sex", "--out"])
        .arg(&output);
    crate::process::run(&mut command).map_err(|error| match error {
        Error::ToolNotFound { .. } => Error::PlinkNotFound {
            executable: tools.plink.clone(),
        },
        other => other,
    })?;
    let path = output.with_extension("hom");
    let mut lines =
        BufReader::new(File::open(&path).map_err(|source| Error::File { path, source })?).lines();
    let header: Vec<_> = lines
        .next()
        .transpose()?
        .unwrap_or_default()
        .split_whitespace()
        .map(str::to_owned)
        .collect();
    let idx = |name| {
        header
            .iter()
            .position(|column| column == name)
            .ok_or_else(|| Error::InvalidInput(format!("PLINK .hom is missing {name}")))
    };
    let chr: usize = idx("CHR")?;
    let start: usize = idx("POS1")?;
    let end: usize = idx("POS2")?;
    let markers: Option<usize> = header.iter().position(|column| column == "NSNP");
    let mut segments: Vec<RohSegment> = Vec::new();
    for line in lines {
        let line: String = line?;
        if line.trim().is_empty() {
            continue;
        }
        let row: Vec<_> = line.split_whitespace().collect();
        let field = |column| {
            row.get(column)
                .copied()
                .ok_or_else(|| Error::InvalidInput("truncated PLINK .hom row".into()))
        };
        let chromosome: Autosome = field(chr)?.parse()?;
        segments.push(RohSegment::new(
            chromosome,
            field(start)?.parse()?,
            field(end)?.parse()?,
            markers
                .map(|idx| {
                    field(idx)?.parse::<NonZeroUsize>().map_err(|_| {
                        Error::InvalidInput("PLINK NSNP must be a positive count".into())
                    })
                })
                .transpose()?,
        )?);
    }
    Ok(segments)
}
