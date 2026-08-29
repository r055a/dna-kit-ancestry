use crate::{
    Error, Result,
    import::{ImportOptions, ToolOptions, containers},
    model::InputFormat,
};
use std::{
    fs::{self, File},
    io::{Read, Seek, SeekFrom},
    path::{Path, PathBuf},
    process::{Child, Command, Stdio},
};
use tempfile::TempDir;

fn label(command: &Command) -> String {
    command.get_program().to_string_lossy().into_owned()
}

fn spawn_error(command: &Command, error: std::io::Error) -> Error {
    if error.kind() == std::io::ErrorKind::NotFound {
        Error::ToolNotFound {
            tool: label(command),
        }
    } else {
        Error::Tool {
            tool: label(command),
            detail: error.to_string(),
        }
    }
}

fn tool_error(tool: String, mut log: File, code: Option<i32>) -> Error {
    let mut bytes: Vec<u8> = Vec::new();
    let _ = log.seek(SeekFrom::Start(0));
    let _ = log.take(16 * 1024).read_to_end(&mut bytes);
    let detail = String::from_utf8_lossy(&bytes);
    Error::Tool {
        tool,
        detail: format!("exit {code:?}: {}", detail.trim()),
    }
}

pub(crate) fn run(command: &mut Command) -> Result<()> {
    let log: File = tempfile::tempfile()?;
    let status = command
        .stdout(log.try_clone()?)
        .stderr(log.try_clone()?)
        .status()
        .map_err(|error| spawn_error(command, error))?;
    if !status.success() {
        return Err(tool_error(label(command), log, status.code()));
    }
    Ok(())
}

fn stop(children: &mut [Child]) {
    for child in children.iter_mut() {
        let _ = child.kill();
    }
    for child in children {
        let _ = child.wait();
    }
}

fn pipeline(mut commands: Vec<Command>, output: &Path) -> Result<()> {
    let log: File = tempfile::tempfile()?;
    let mut children: Vec<Child> = Vec::new();
    let mut prev = None;
    let ttl: usize = commands.len();
    let result = (|| {
        for (idx, command) in commands.iter_mut().enumerate() {
            if let Some(stdout) = prev.take() {
                command.stdin(Stdio::from(stdout));
            }
            command.stderr(log.try_clone()?);
            if idx + 1 == ttl {
                command.stdout(File::create(output)?);
            } else {
                command.stdout(Stdio::piped());
            }
            match command.spawn() {
                Ok(mut child) => {
                    prev = child.stdout.take();
                    children.push(child);
                }
                Err(error) => {
                    return Err(spawn_error(command, error));
                }
            }
        }
        let mut failed = None;
        for idx in (0..children.len()).rev() {
            let status = match children[idx].wait() {
                Ok(status) => status,
                Err(error) => {
                    return Err(Error::Tool {
                        tool: label(&commands[idx]),
                        detail: error.to_string(),
                    });
                }
            };
            if !status.success() {
                failed = Some((label(&commands[idx]), status.code()));
                break;
            }
        }
        if let Some((tool, code)) = failed {
            return Err(tool_error(tool, log, code));
        }
        Ok(())
    })();
    stop(&mut children);
    result
}

pub(crate) struct Converted {
    pub path: PathBuf,
    pub steps: Vec<String>,
    _directory: TempDir,
}

fn req(path: &Option<PathBuf>, format: &'static str, requirement: &'static str) -> Result<PathBuf> {
    let path: &PathBuf = path.as_ref().ok_or(Error::NeedsConversion {
        format,
        requirement,
    })?;
    fs::canonicalize(path).map_err(|source| Error::File {
        path: path.clone(),
        source,
    })
}

pub(crate) fn to_vcf(
    path: &Path,
    format: InputFormat,
    options: &ImportOptions,
    tools: &ToolOptions,
) -> Result<Converted> {
    let dir = tempfile::tempdir()?;
    let primary: PathBuf = containers::materialize(
        path,
        dir.path(),
        options.max_bytes,
        matches!(format, InputFormat::Bam | InputFormat::Bcf),
    )?;
    let mut output: PathBuf = dir.path().join("calls.vcf");
    let mut steps: Vec<String> = Vec::new();
    match format {
        InputFormat::Bcf => {
            let mut command: Command = Command::new(&tools.bcftools);
            command.args(["view", "-Ov", "--"]).arg(&primary);
            pipeline(vec![command], &output)?;
            steps.push("bcftools view (BCF to VCF)".into());
        }
        InputFormat::Pgen | InputFormat::Bgen => {
            let mut command: Command = Command::new(&tools.plink2);
            if format == InputFormat::Pgen {
                let pvar: PathBuf = containers::materialize(
                    &containers::companion(path, "pvar")?,
                    dir.path(),
                    options.max_bytes,
                    false,
                )?;
                let psam: PathBuf = containers::materialize(
                    &containers::companion(path, "psam")?,
                    dir.path(),
                    options.max_bytes,
                    false,
                )?;
                command
                    .arg("--pgen")
                    .arg(&primary)
                    .arg("--pvar")
                    .arg(pvar)
                    .arg("--psam")
                    .arg(psam);
            } else {
                command.arg("--bgen").arg(&primary).arg("ref-first");
                if let Some(sample) = containers::optional_companion(path, "sample")? {
                    let sample =
                        containers::materialize(&sample, dir.path(), options.max_bytes, false)?;
                    command.arg("--sample").arg(sample);
                }
            }
            let prefix: PathBuf = dir.path().join("converted");
            command
                .args([
                    "--chr",
                    "1-22",
                    "--export",
                    "vcf",
                    "id-paste=fid,iid",
                    "id-delim=:",
                    "--hard-call-threshold",
                    "0.1",
                    "--out",
                ])
                .arg(&prefix);
            command.arg("--threads").arg(tools.threads.to_string());
            run(&mut command)?;
            output = prefix.with_extension("vcf");
            if !output.is_file() {
                return Err(Error::InvalidInput("PLINK 2 did not produce a VCF".into()));
            }
            steps.push("PLINK 2 hard-call VCF export (threshold 0.1)".into());
        }
        InputFormat::Sam | InputFormat::Bam | InputFormat::Cram | InputFormat::Fastq => {
            let reference: PathBuf = req(
                &tools.reference,
                "reads/alignments",
                "--reference FASTA and --sites SNP-panel positions; reads do not contain genotype calls",
            )?;
            let sites: PathBuf = req(
                &tools.sites,
                "reads/alignments",
                "--sites CHROM/POS or BED SNP-panel file, matching the reference assembly",
            )?;
            let aligned: PathBuf = dir.path().join("sorted.bam");
            let mut sort: Command = Command::new(&tools.samtools);
            sort.args(["sort", "-@"])
                .arg(tools.threads.to_string())
                .arg("-o")
                .arg(&aligned);
            if format == InputFormat::Fastq {
                let mut align: Command = Command::new(&tools.bwa);
                align
                    .args(["mem", "-t"])
                    .arg(tools.threads.to_string())
                    .arg("-R")
                    .arg(r"@RG\tID:DNA\tSM:SAMPLE")
                    .arg(&reference)
                    .arg(&primary);
                if let Some(mate) = &tools.mate {
                    let mate: PathBuf = fs::canonicalize(mate).map_err(|source| Error::File {
                        path: mate.clone(),
                        source,
                    })?;
                    let mate_dir: PathBuf = dir.path().join("mate");
                    fs::create_dir(&mate_dir)?;
                    align.arg(containers::materialize(
                        &mate,
                        &mate_dir,
                        options.max_bytes,
                        false,
                    )?);
                }
                sort.arg("-");
                pipeline(vec![align, sort], &dir.path().join("sort.stdout"))?;
                steps.push("bwa mem and samtools sort".into());
            } else {
                if format == InputFormat::Cram {
                    sort.arg("--reference").arg(&reference);
                }
                sort.arg(&primary);
                run(&mut sort)?;
                steps.push("samtools sort".into());
            }
            let mut pileup: Command = Command::new(&tools.bcftools);
            pileup
                .args(["mpileup", "-Ou", "-f"])
                .arg(&reference)
                .arg("-T")
                .arg(&sites)
                .args(["-q", "20", "-Q", "20", "-a", "FORMAT/DP,FORMAT/AD"])
                .arg(&aligned);
            let mut call: Command = Command::new(&tools.bcftools);
            call.args(["call", "-m", "-f", "GQ", "-Ov"]);
            pipeline(vec![pileup, call], &output)?;
            steps.push(
                "bcftools mpileup/call at supplied SNP sites, including reference calls".into(),
            );
        }
        _ => {
            return Err(Error::InvalidInput(
                "no external bridge for this format".into(),
            ));
        }
    }
    Ok(Converted {
        path: output,
        steps,
        _directory: dir,
    })
}
