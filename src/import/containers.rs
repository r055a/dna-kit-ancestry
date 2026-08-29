use super::ImportOptions;
use crate::{Error, Result, units::InputByteLimit};
use std::{
    fs::{self, File},
    io::{self, BufRead, BufReader, Read, Write},
    path::{Component, Path, PathBuf},
};
use tempfile::TempDir;

const COMPRESS_SFX: [&str; 7] = [".gz", ".bgz", ".bgzf", ".bz2", ".xz", ".zst", ".zstd"];

fn wrapped_names(stem: String) -> impl Iterator<Item = String> {
    std::iter::once(stem.clone()).chain(
        COMPRESS_SFX
            .into_iter()
            .map(move |suffix| format!("{stem}{suffix}")),
    )
}

pub(super) struct PreparedInput {
    pub path: PathBuf,
    pub member: Option<String>,
    _directories: Vec<TempDir>,
}

struct Limited<R> {
    inner: R,
    remaining: u64,
}
impl<R: Read> Read for Limited<R> {
    fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
        if bytes.is_empty() {
            return Ok(0);
        }
        if self.remaining == 0 {
            return if self.inner.read(&mut [0])? == 0 {
                Ok(0)
            } else {
                Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "decompressed input exceeds --max-input-bytes",
                ))
            };
        }
        let len: usize = bytes
            .len()
            .min(self.remaining.min(usize::MAX as u64) as usize);
        let count: usize = self.inner.read(&mut bytes[..len])?;
        self.remaining -= count as u64;
        Ok(count)
    }
}

pub(crate) fn binary_reader(path: &Path, maximum: InputByteLimit) -> Result<Box<dyn BufRead>> {
    let file: File = File::open(path).map_err(|source| Error::File {
        path: path.into(),
        source,
    })?;
    let mut reader: Box<dyn BufRead> = Box::new(BufReader::new(file));
    for depth in 0..=4 {
        let magic: &[u8] = reader.fill_buf()?;
        let decoded: Box<dyn Read> = if magic.starts_with(&[0x1f, 0x8b]) {
            Box::new(flate2::read::MultiGzDecoder::new(reader))
        } else if magic.starts_with(b"BZh") {
            Box::new(bzip2::read::MultiBzDecoder::new(reader))
        } else if magic.starts_with(&[0xfd, b'7', b'z', b'X', b'Z', 0]) {
            Box::new(lzma_rust2::XzReader::new(reader, true))
        } else if magic.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]) {
            Box::new(zstd::stream::read::Decoder::new(reader)?)
        } else {
            return Ok(Box::new(BufReader::new(Limited {
                inner: reader,
                remaining: maximum.get(),
            })));
        };
        if depth == 4 {
            break;
        }
        reader = Box::new(BufReader::new(decoded));
    }
    Err(Error::InvalidInput(
        "compression nesting exceeds four layers".into(),
    ))
}

pub(crate) fn text_reader(path: &Path, maximum: InputByteLimit) -> Result<Box<dyn BufRead>> {
    let reader = binary_reader(path, maximum)?;
    Ok(Box::new(BufReader::new(
        encoding_rs_io::DecodeReaderBytesBuilder::new().build(reader),
    )))
}

pub(crate) fn uncompressed_name(name: &str) -> &str {
    let mut cur: &str = name;
    for _ in 0..4 {
        let nxt: Option<&str> = COMPRESS_SFX.iter().find_map(|sfx| {
            let split: usize = cur.len().checked_sub(sfx.len())?;
            cur.get(split..)?
                .eq_ignore_ascii_case(sfx)
                .then(|| &cur[..split])
        });
        match nxt {
            Some(nxt) => cur = nxt,
            None => break,
        }
    }
    cur
}

fn candidate(name: &str) -> bool {
    let lower: String = name.to_lowercase();
    if lower
        .split('/')
        .any(|part| part.starts_with('.') || part == "__macosx")
    {
        return false;
    }
    let stem: &str = lower.rsplit('/').next().unwrap_or(&lower);
    if stem.starts_with("readme") || stem.starts_with("license") {
        return false;
    }
    let ext: &str = uncompressed_name(&lower)
        .rsplit('.')
        .next()
        .unwrap_or_default();
    matches!(
        ext,
        "txt"
            | "csv"
            | "tsv"
            | "vcf"
            | "gvcf"
            | "bcf"
            | "bed"
            | "ped"
            | "tped"
            | "pgen"
            | "bgen"
            | "bam"
            | "sam"
            | "cram"
            | "fastq"
            | "fq"
            | "fasta"
            | "fa"
            | "zip"
            | "tar"
    )
}

fn choose(names: &[String], req: Option<&str>) -> Result<String> {
    if let Some(req) = req {
        if names.iter().filter(|name| name.as_str() == req).count() == 1 {
            return Ok(req.into());
        }
        return Err(Error::InvalidInput(format!(
            "archive member {req:?} is missing or duplicated"
        )));
    }
    let candidates: Vec<_> = names.iter().filter(|name| candidate(name)).collect();
    if candidates.len() == 1 {
        return Ok(candidates[0].to_owned());
    }
    Err(Error::Ambiguous {
        kind: "archive member",
        option: "--member",
        choices: candidates
            .into_iter()
            .take(30)
            .cloned()
            .collect::<Vec<_>>()
            .join(", "),
    })
}

fn sidecars(name: &str) -> Vec<String> {
    let lower: String = name.to_lowercase();
    let plain: &str = uncompressed_name(&lower);
    let (stem, ext) = plain.rsplit_once('.').unwrap_or((plain, ""));
    let suffixes: &[&str] = match ext {
        "ped" => &["map"],
        "tped" => &["tfam"],
        "bed" => &["bim", "fam"],
        "pgen" => &["pvar", "psam"],
        "bgen" => &["sample"],
        _ => &[],
    };
    suffixes
        .iter()
        .flat_map(|suffix| wrapped_names(format!("{stem}.{suffix}")))
        .collect()
}

fn dest(dir: &Path, name: &str) -> Result<PathBuf> {
    let path: &Path = Path::new(name);
    if path.components().any(|part| {
        matches!(
            part,
            Component::ParentDir | Component::RootDir | Component::Prefix(_)
        )
    }) || name.contains('\\')
    {
        return Err(Error::InvalidInput("unsafe archive member path".into()));
    }
    let output: PathBuf = dir.join(path);
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }
    Ok(output)
}

fn copy_limited(reader: impl Read, output: &Path, remaining: &mut u64) -> Result<()> {
    let mut input = reader.take(remaining.saturating_add(1));
    let mut writer: File = File::options().write(true).create_new(true).open(output)?;
    let count: u64 = io::copy(&mut input, &mut writer)?;
    writer.flush()?;
    if count > *remaining {
        return Err(Error::InvalidInput(
            "archive contents exceed --max-input-bytes".into(),
        ));
    }
    *remaining -= count;
    Ok(())
}

fn unpack_zip(
    path: &Path,
    dir: &Path,
    options: &ImportOptions,
    req: Option<&str>,
) -> Result<(PathBuf, String)> {
    let mut archive = zip::ZipArchive::new(File::open(path)?)?;
    if archive.has_duplicate_names() {
        return Err(Error::InvalidInput(
            "duplicate ZIP member names are ambiguous".into(),
        ));
    }
    let names: Vec<String> = (0..archive.len())
        .map(|idx: usize| {
            archive
                .name_for_index(idx)
                .transpose()?
                .map(|name| name.into_owned())
                .ok_or(zip::result::ZipError::FileNotFound)
        })
        .collect::<std::result::Result<_, _>>()?;
    let selected: String = choose(&names, req)?;
    let companions: Vec<String> = sidecars(&selected);
    let password: Option<String> = std::env::var("DNA_KIT_ANCESTRY_ARCHIVE_PASSWORD").ok();
    let mut remaining: u64 = options.max_bytes.get();
    for idx in 0..archive.len() {
        let name: &String = names
            .get(idx)
            .ok_or_else(|| Error::InvalidInput("invalid ZIP index".into()))?;
        if name != &selected && !companions.contains(&name.to_lowercase()) {
            continue;
        }
        let mut entry = match password.as_deref() {
            Some(password) => archive.by_index_decrypt(idx, password.as_bytes())?,
            None => archive.by_index(idx)?,
        };
        if entry.is_dir()
            || entry
                .unix_mode()
                .is_some_and(|mode: u32| mode & 0o170000 == 0o120000)
        {
            return Err(Error::InvalidInput(
                "archive genotype members must be regular files".into(),
            ));
        }
        let output: PathBuf = dest(dir, name)?;
        copy_limited(&mut entry, &output, &mut remaining)?;
    }
    Ok((dir.join(&selected), selected))
}

fn unpack_tar(
    path: &Path,
    dir: &Path,
    options: &ImportOptions,
    req: Option<&str>,
) -> Result<(PathBuf, String)> {
    let mut archive = tar::Archive::new(binary_reader(path, options.max_bytes)?);
    let mut names: Vec<String> = Vec::new();
    for entry in archive.entries()? {
        let entry = entry?;
        if entry.header().entry_type().is_file() {
            names.push(entry.path()?.to_string_lossy().into_owned());
        }
    }
    let selected: String = choose(&names, req)?;
    let companions: Vec<String> = sidecars(&selected);
    let mut remaining: u64 = options.max_bytes.get();
    let mut archive = tar::Archive::new(binary_reader(path, options.max_bytes)?);
    for entry in archive.entries()? {
        let mut entry = entry?;
        let name: String = entry.path()?.to_string_lossy().into_owned();
        if name != selected && !companions.contains(&name.to_lowercase()) {
            continue;
        }
        if !entry.header().entry_type().is_file() {
            return Err(Error::InvalidInput(
                "archive genotype members must be regular files".into(),
            ));
        }
        let output: PathBuf = dest(dir, &name)?;
        copy_limited(&mut entry, &output, &mut remaining)?;
    }
    Ok((dir.join(&selected), selected))
}

pub(super) fn prepare(path: &Path, options: &ImportOptions) -> Result<PreparedInput> {
    let mut prepared = PreparedInput {
        path: fs::canonicalize(path).map_err(|source| Error::File {
            path: path.into(),
            source,
        })?,
        member: None,
        _directories: Vec::new(),
    };
    let req: Vec<_> = options
        .member
        .as_deref()
        .map(|name| name.split('!').collect())
        .unwrap_or_default();
    for depth in 0..=4 {
        let mut raw = BufReader::new(File::open(&prepared.path)?);
        let is_raw_zip: bool = raw.fill_buf()?.starts_with(b"PK\x03\x04");
        let mut decoded = binary_reader(&prepared.path, options.max_bytes)?;
        let prefix: &[u8] = decoded.fill_buf()?;
        let is_zip: bool = prefix.starts_with(b"PK\x03\x04");
        let is_tar: bool = prefix.get(257..262) == Some(b"ustar")
            || uncompressed_name(&prepared.path.to_string_lossy().to_lowercase()).ends_with(".tar");
        if !is_zip && !is_tar {
            if depth < req.len() {
                return Err(Error::InvalidInput(
                    "--member has an unused archive selection".into(),
                ));
            }
            return Ok(prepared);
        }
        if depth == 4 {
            break;
        }
        if is_zip && !is_raw_zip {
            let dir: TempDir = tempfile::tempdir()?;
            let path: PathBuf = dir.path().join("decoded.zip");
            let mut remaining: u64 = options.max_bytes.get();
            copy_limited(
                binary_reader(&prepared.path, options.max_bytes)?,
                &path,
                &mut remaining,
            )?;
            prepared.path = path;
            prepared._directories.push(dir);
        }
        let dir: TempDir = tempfile::tempdir()?;
        let selection: Option<&str> = req.get(depth).copied();
        let (path, member) = if is_zip {
            unpack_zip(&prepared.path, dir.path(), options, selection)?
        } else {
            unpack_tar(&prepared.path, dir.path(), options, selection)?
        };
        prepared.path = path;
        prepared.member = Some(match prepared.member {
            Some(outer) => format!("{outer}!{member}"),
            None => member,
        });
        prepared._directories.push(dir);
    }
    Err(Error::InvalidInput(
        "archive nesting exceeds four layers".into(),
    ))
}

fn companion_path(path: &Path, extension: &str) -> PathBuf {
    let stem: PathBuf = PathBuf::from(uncompressed_name(&path.to_string_lossy()));
    stem.with_extension(extension)
}

pub(crate) fn companion(path: &Path, extension: &str) -> Result<PathBuf> {
    optional_companion(path, extension)?.ok_or_else(|| {
        Error::InvalidInput(format!(
            "missing companion {}",
            companion_path(path, extension).display()
        ))
    })
}

pub(crate) fn optional_companion(path: &Path, extension: &str) -> Result<Option<PathBuf>> {
    let plain: PathBuf = companion_path(path, extension);
    for candidate in wrapped_names(plain.to_string_lossy().into_owned()).map(PathBuf::from) {
        if candidate.is_file() {
            return Ok(Some(candidate));
        }
    }
    if let Some(dir) = plain.parent() {
        let names: Vec<_> = wrapped_names(
            plain
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .to_lowercase(),
        )
        .collect();
        let matches: Vec<_> = fs::read_dir(dir)?
            .filter_map(std::result::Result::ok)
            .filter(|entry| names.contains(&entry.file_name().to_string_lossy().to_lowercase()))
            .map(|entry| entry.path())
            .collect();
        if matches.len() == 1 {
            return Ok(Some(matches[0].clone()));
        }
        if matches.len() > 1 {
            return Err(Error::InvalidInput(format!(
                "ambiguous companions for {}",
                path.display()
            )));
        }
    }
    Ok(None)
}

pub(crate) fn materialize(
    path: &Path,
    dir: &Path,
    maximum: InputByteLimit,
    keep_bgzf: bool,
) -> Result<PathBuf> {
    let mut file: BufReader<File> = BufReader::new(File::open(path)?);
    let magic: &[u8] = file.fill_buf()?;
    let is_compressed: bool = magic.starts_with(&[0x1f, 0x8b])
        || magic.starts_with(b"BZh")
        || magic.starts_with(&[0xfd, b'7', b'z', b'X', b'Z', 0])
        || magic.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]);
    if !is_compressed || (keep_bgzf && magic.get(12..14) == Some(b"BC")) {
        return Ok(path.into());
    }
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    let output: PathBuf = dir.join(uncompressed_name(&name));
    let mut remaining: u64 = maximum.get();
    copy_limited(binary_reader(path, maximum)?, &output, &mut remaining)?;
    Ok(output)
}
