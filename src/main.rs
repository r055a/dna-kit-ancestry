use clap::Parser;
use dna_kit_ancestry::{
    Error, Result,
    analysis::{self, Denominator, Method, RohConfig},
    import::{self, ColumnMapping, Delimiter, ImportOptions, SampleSelection, ToolOptions},
    install,
    model::*,
    units::{BasePairLength, InputByteLimit, PhredScore, ReadDepth, SnpDensity},
};
use serde::Serialize;
use std::{
    ffi::OsString,
    fs::File,
    io::{self, Write},
    num::{NonZeroU64, NonZeroUsize},
    path::{Path, PathBuf},
    process::ExitCode,
};

#[derive(Parser)]
#[command(
    version,
    about = "Import raw DNA genotype data and estimate runs of homozygosity"
)]
struct Args {
    /// DNA download, archive, or primary PLINK file
    input: Option<PathBuf>,
    /// Legacy-compatible alternative to the positional input
    #[arg(long, conflicts_with = "input")]
    file: Option<PathBuf>,
    #[arg(long)]
    formats: bool,
    /// Download and verify PLINK 1.9 into the current user's data directory
    #[arg(long, conflicts_with_all = ["input", "file", "formats", "inspect", "export"])]
    install_plink: bool,
    /// Import and report QC without running ROH analysis
    #[arg(long)]
    inspect: bool,
    #[arg(long)]
    json: bool,
    #[arg(long, value_enum, default_value = "auto")]
    format: InputFormat,
    #[arg(long, value_enum, default_value = "plink")]
    method: Method,
    #[arg(long)]
    assembly: Option<Assembly>,
    #[arg(long, conflicts_with = "sample_index")]
    sample: Option<SampleId>,
    #[arg(long, conflicts_with = "sample")]
    sample_index: Option<NonZeroUsize>,
    #[arg(long)]
    member: Option<String>,
    /// One-based mapping, e.g. id=1,chr=2,pos=3,gt=4
    #[arg(long)]
    cols: Option<ColumnMapping>,
    #[arg(long, value_enum, default_value = "auto")]
    delimiter: Delimiter,
    /// Explicit allele interpretation; IUPAC ambiguity codes require iupac
    #[arg(long, value_enum)]
    coding: Option<CallCoding>,
    #[arg(long)]
    strict: bool,
    #[arg(long)]
    include_filtered: bool,
    #[arg(long)]
    min_depth: Option<ReadDepth>,
    #[arg(long)]
    min_gq: Option<PhredScore>,
    #[arg(long, default_value = "17179869184")]
    max_input_bytes: InputByteLimit,
    /// Write sorted, deduplicated autosomal calls; destination must not exist
    #[arg(long)]
    export: Option<PathBuf>,
    #[arg(long, default_value = "100")]
    min_snps: NonZeroUsize,
    #[arg(long, default_value = "1000000")]
    min_length_bp: BasePairLength,
    #[arg(long, default_value = "1000000")]
    max_gap_bp: BasePairLength,
    #[arg(long, default_value = "1")]
    max_heterozygous: usize,
    #[arg(long, default_value = "50")]
    max_density_kb: NonZeroU64,
    #[arg(long, default_value = "50")]
    window_snps: NonZeroUsize,
    /// Override the observed-marker-span denominator
    #[arg(long)]
    denominator_bp: Option<BasePairLength>,
    #[arg(long)]
    plink: Option<PathBuf>,
    #[arg(long, default_value = "plink2")]
    plink2: PathBuf,
    #[arg(long, default_value = "bcftools")]
    bcftools: PathBuf,
    #[arg(long, default_value = "samtools")]
    samtools: PathBuf,
    #[arg(long, default_value = "bwa")]
    bwa: PathBuf,
    /// Matching, idx reference FASTA for read/alignment inputs
    #[arg(long)]
    reference: Option<PathBuf>,
    /// Matching SNP-panel positions for read/alignment inputs
    #[arg(long)]
    sites: Option<PathBuf>,
    #[arg(long)]
    mate: Option<PathBuf>,
    #[arg(long, default_value = "1")]
    threads: NonZeroUsize,
}

fn legacy_args() -> impl Iterator<Item = OsString> {
    std::env::args_os().map(|arg: OsString| match arg.to_str() {
        Some("-file") => "--file".into(),
        Some("-min_snp") => "--min-snps".into(),
        Some("-min_bp_len") => "--min-length-bp".into(),
        Some("-max_bp_range") => "--max-gap-bp".into(),
        Some("-max_snp_het") => "--max-heterozygous".into(),
        Some("-max_snp_avg_range") => "--max-density-kb".into(),
        _ => arg,
    })
}

fn export(path: &Path, dataset: &Dataset) -> Result<()> {
    let mut output: File = File::options().write(true).create_new(true).open(path)?;
    writeln!(output, "# assembly: {:?}", dataset.src().assembly)?;
    if let Some(sample) = &dataset.src().sample {
        writeln!(output, "# sample: {sample}")?;
    }
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(output);
    writer.write_record(["rsid", "chromosome", "position", "genotype", "coding"])?;
    for snp in dataset.snps() {
        let coding = match snp.genotype {
            Genotype::Dna(_) => "dna",
            Genotype::Assay(_) => "assay",
            Genotype::Numeric(_) => "numeric",
        };
        writer.write_record([
            snp.id.to_string(),
            snp.chromosome.to_string(),
            snp.position.to_string(),
            snp.genotype.to_string(),
            coding.into(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

fn write_json(val: &impl Serialize, output: &mut impl Write) -> Result<()> {
    serde_json::to_writer_pretty(&mut *output, val)?;
    writeln!(output)?;
    Ok(())
}

fn import_summary(dataset: &Dataset, output: &mut impl Write) -> Result<()> {
    writeln!(
        output,
        "Format: {:?}; assembly: {:?}",
        dataset.src().format,
        dataset.src().assembly
    )?;
    writeln!(
        output,
        "Usable autosomal SNPs: {}",
        dataset.stats().accepted_snps
    )?;
    writeln!(
        output,
        "Excluded: {} missing, {} filtered, {} non-autosomal, {} non-SNP, {} haploid/polyploid, {} gVCF blocks, {} malformed",
        dataset.stats().missing_calls,
        dataset.stats().filtered,
        dataset.stats().non_autosomal,
        dataset.stats().non_snps,
        dataset.stats().haploid_or_polyploid,
        dataset.stats().reference_blocks,
        dataset.stats().malformed
    )?;
    writeln!(
        output,
        "Duplicates: {}; conflicting loci removed: {}",
        dataset.stats().duplicate_rows,
        dataset.stats().conflicting_loci
    )?;
    Ok(())
}

fn run(args: Args, output: &mut impl Write) -> Result<()> {
    if args.install_plink {
        let path = install::install_plink()?;
        if args.json {
            write_json(&path, output)?;
        } else {
            writeln!(
                output,
                "PLINK {} ready: {}",
                install::PLINK_VERSION,
                path.display()
            )?;
        }
        return Ok(());
    }
    if args.formats {
        return Ok(());
    }
    let path: PathBuf = args
        .input
        .or(args.file)
        .ok_or_else(|| Error::InvalidInput("provide a DNA file or use --formats".into()))?;
    let options = ImportOptions {
        format: args.format,
        assembly: args.assembly.unwrap_or_default(),
        sample: match (args.sample, args.sample_index) {
            (Some(id), _) => SampleSelection::Name(id),
            (_, Some(idx)) => SampleSelection::Index(idx),
            _ => SampleSelection::Single,
        },
        member: args.member,
        cols: args.cols,
        delimiter: args.delimiter,
        coding: args.coding,
        strict: args.strict,
        include_filtered: args.include_filtered,
        min_depth: args.min_depth,
        min_genotype_quality: args.min_gq,
        max_bytes: args.max_input_bytes,
    };
    let tools = ToolOptions {
        plink: args.plink.unwrap_or_else(install::default_command),
        plink2: args.plink2,
        bcftools: args.bcftools,
        samtools: args.samtools,
        bwa: args.bwa,
        reference: args.reference,
        sites: args.sites,
        mate: args.mate,
        threads: args.threads,
    };
    let dataset: Dataset = import::load(&path, &options, &tools)?;
    if let Some(path) = &args.export {
        export(path, &dataset)?;
    }
    if args.inspect {
        if args.json {
            #[derive(Serialize)]
            struct Inspection<'a> {
                src: &'a Source,
                import: &'a ImportStats,
            }
            write_json(
                &Inspection {
                    src: dataset.src(),
                    import: dataset.stats(),
                },
                output,
            )?;
        } else {
            import_summary(&dataset, output)?;
        }
        return Ok(());
    }
    let density: SnpDensity = SnpDensity::from_kilobases(args.max_density_kb)?;
    let config = RohConfig {
        min_snps: args.min_snps,
        min_length_bp: args.min_length_bp,
        max_gap_bp: args.max_gap_bp,
        max_heterozygous: args.max_heterozygous,
        max_density_bp_per_snp: density,
        window_snps: args.window_snps,
    };
    let report = analysis::analyze(
        &dataset,
        &config,
        args.method,
        args.denominator_bp
            .map(Denominator::Provided)
            .unwrap_or_default(),
        &tools,
    )?;
    if args.json {
        write_json(&report, output)?;
    } else {
        import_summary(&dataset, output)?;
        writeln!(
            output,
            "Method: {:?}; homozygous: {}; heterozygous: {}",
            report.method(),
            report.homozygous_snps(),
            report.heterozygous_snps()
        )?;
        writeln!(
            output,
            "ROH: {} segments, {} bp",
            report.segments().len(),
            report.total_roh_bp()
        )?;
        writeln!(
            output,
            "Denominator: {:?}, {} bp",
            report.denominator(),
            report.denominator_bp()
        )?;
        writeln!(
            output,
            "F_ROH: {:.6} ({:.2}%)",
            report.f_roh().value(),
            report.f_roh().percent()
        )?;
        for note in &dataset.src().notes {
            writeln!(output, "Note: {note}")?;
        }
    }
    Ok(())
}

fn main() -> ExitCode {
    match run(Args::parse_from(legacy_args()), &mut io::stdout().lock()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(Error::Io(error)) if error.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("dna-kit: {error}");
            ExitCode::FAILURE
        }
    }
}
