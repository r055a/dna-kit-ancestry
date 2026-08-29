use sha2::{Digest, Sha256};
use std::{
    env,
    fs::{self, File},
    io::{self, Cursor, Read, Write},
    path::{Path, PathBuf},
    time::Duration,
};
use tempfile::NamedTempFile;
use thiserror::Error;
use ureq::{
    Agent,
    tls::{RootCerts, TlsConfig},
};

pub const PLINK_VERSION: &str = "1.9.0";
const MAX_BYTES: u64 = 64 * 1024 * 1024;
type Result<T> = std::result::Result<T, InstallError>;

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum InstallError {
    #[error(
        "automatic PLINK installation is unsupported on {os}/{arch}; install PLINK 1.9 for this platform and use --plink PATH"
    )]
    UnsupportedPlatform {
        os: &'static str,
        arch: &'static str,
    },
    #[error("cannot determine a user data directory; set DNA_KIT_ANCESTRY_TOOLS_DIR to an absolute path")]
    MissingDirectory,
    #[error("the tools directory must be absolute: {0}")]
    RelativeDirectory(PathBuf),
    #[error("PLINK {target} failed SHA-256 verification")]
    Integrity { target: &'static str },
    #[error("refusing to replace an unexpected existing file: {0}; remove it before reinstalling")]
    ExistingFile(PathBuf),
    #[error("PLINK download or executable exceeds the 64 MiB limit")]
    TooLarge,
    #[error("PLINK archive must contain one regular executable with the expected name")]
    InvalidArchive,
    #[error("PLINK download failed: {0}")]
    Download(#[from] ureq::Error),
    #[error(transparent)]
    Io(#[from] io::Error),
    #[error(transparent)]
    Zip(#[from] zip::result::ZipError),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Checksum([u8; 32]);

impl Checksum {
    const fn hex(value: &str) -> Self {
        const fn digit(value: u8) -> u8 {
            match value {
                b'0'..=b'9' => value - b'0',
                b'a'..=b'f' => value - b'a' + 10,
                _ => panic!("invalid checksum literal"),
            }
        }
        let bytes = value.as_bytes();
        assert!(bytes.len() == 64, "invalid checksum length");
        let mut output = [0; 32];
        let mut i = 0;
        while i < output.len() {
            output[i] = digit(bytes[2 * i]) * 16 + digit(bytes[2 * i + 1]);
            i += 1;
        }
        Self(output)
    }
    fn of(bytes: &[u8]) -> Self {
        Self(Sha256::digest(bytes).into())
    }
    fn verify(self, bytes: &[u8], target: &'static str) -> Result<()> {
        if self == Self::of(bytes) {
            Ok(())
        } else {
            Err(InstallError::Integrity { target })
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Asset {
    url: &'static str,
    executable: &'static str,
    archive_checksum: Checksum,
    binary_checksum: Checksum,
}

impl Asset {
    fn current() -> Result<Self> {
        Self::for_platform(env::consts::OS, env::consts::ARCH)
    }

    fn for_platform(os: &'static str, arch: &'static str) -> Result<Self> {
        // Stable binaries linked from https://www.cog-genomics.org/plink/1.9/.
        let (url, executable, archive_checksum, binary_checksum) = match (os, arch) {
            ("linux", "x86_64") => (
                "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20260927.zip",
                "plink",
                Checksum::hex("ee4d8616dc64c04ed395e81a7bb04f9b3d92e4bf23d12a8b6233b19f7227eccd"),
                Checksum::hex("47ce3069a1349128bdc5193aa8af33ec20dc99ece920843c80453c34eea1a062"),
            ),
            ("macos", "x86_64" | "aarch64") => (
                "https://s3.amazonaws.com/plink1-assets/plink_mac_20260927.zip",
                "plink",
                Checksum::hex("927ecadbefb3e0ceefbe4978d71bf63566e9ea422d6595c2ee23a83ad76eeb62"),
                Checksum::hex("13a7fc9b908d0feb968d0835108b7fa5416fe27f44cad9bb7a0887676a644c40"),
            ),
            ("windows", "x86_64") => (
                "https://s3.amazonaws.com/plink1-assets/plink_win64_20260927.zip",
                "plink.exe",
                Checksum::hex("625538e0062eb9fad6c26a47c9f36b19d3ec0afe5053ec68ff5ff877a56676c6"),
                Checksum::hex("efc55270aba7216e4fcebdf781ee04e9fa556c2fe5be1f71c836ba73930a6b6b"),
            ),
            _ => return Err(InstallError::UnsupportedPlatform { os, arch }),
        };
        Ok(Self {
            url,
            executable,
            archive_checksum,
            binary_checksum,
        })
    }
}

fn absolute(path: PathBuf) -> Result<PathBuf> {
    if path.is_absolute() {
        Ok(path)
    } else {
        Err(InstallError::RelativeDirectory(path))
    }
}

fn tools_directory() -> Result<PathBuf> {
    if let Some(path) = env::var_os("DNA_KIT_ANCESTRY_TOOLS_DIR") {
        return absolute(path.into());
    }
    let home = || env::var_os("HOME").map(PathBuf::from);
    let data = match env::consts::OS {
        "windows" => env::var_os("LOCALAPPDATA").map(PathBuf::from),
        "macos" => home().map(|home| home.join("Library/Application Support")),
        _ => env::var_os("XDG_DATA_HOME")
            .filter(|path| !path.is_empty())
            .map(PathBuf::from)
            .or_else(|| home().map(|home| home.join(".local/share"))),
    }
    .ok_or(InstallError::MissingDirectory)?;
    absolute(data.join("dna-kit/tools"))
}

fn managed_executable(asset: Asset) -> Result<PathBuf> {
    Ok(tools_directory()?
        .join(format!("plink-{PLINK_VERSION}"))
        .join(format!("{}-{}", env::consts::OS, env::consts::ARCH))
        .join(asset.executable))
}

pub fn default_command() -> PathBuf {
    if let Ok(asset) = Asset::current()
        && let Ok(path) = managed_executable(asset)
        && path.is_file()
    {
        return path;
    }
    if let Some(paths) = env::var_os("PATH") {
        for name in ["plink1.9", "plink"] {
            let name = format!("{name}{}", env::consts::EXE_SUFFIX);
            for directory in env::split_paths(&paths) {
                let path = directory.join(&name);
                if path.is_file() {
                    return path;
                }
            }
        }
    }
    "plink".into()
}

fn read_bounded(reader: impl Read) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    reader.take(MAX_BYTES + 1).read_to_end(&mut bytes)?;
    if bytes.len() as u64 > MAX_BYTES {
        return Err(InstallError::TooLarge);
    }
    Ok(bytes)
}

fn existing(path: &Path, checksum: Checksum) -> Result<bool> {
    let metadata = match fs::symlink_metadata(path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(false),
        Err(error) => return Err(error.into()),
    };
    if !metadata.is_file()
        || metadata.len() > MAX_BYTES
        || Checksum::of(&read_bounded(File::open(path)?)?) != checksum
    {
        return Err(InstallError::ExistingFile(path.into()));
    }
    Ok(true)
}

fn publish(bytes: &[u8], asset: Asset, destination: &Path) -> Result<()> {
    asset.archive_checksum.verify(bytes, "archive")?;
    let mut archive = zip::ZipArchive::new(Cursor::new(bytes))?;
    if archive.has_duplicate_names() {
        return Err(InstallError::InvalidArchive);
    }
    let mut entry = archive.by_name(asset.executable)?;
    if entry.is_dir()
        || entry
            .unix_mode()
            .is_some_and(|mode| mode & 0o170000 == 0o120000)
    {
        return Err(InstallError::InvalidArchive);
    }
    let binary = read_bounded(&mut entry)?;
    asset.binary_checksum.verify(&binary, "executable")?;
    let parent = destination.parent().ok_or(InstallError::MissingDirectory)?;
    fs::create_dir_all(parent)?;
    let mut pending = NamedTempFile::new_in(parent)?;
    pending.write_all(&binary)?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        pending
            .as_file()
            .set_permissions(fs::Permissions::from_mode(0o755))?;
    }
    pending.as_file().sync_all()?;
    if let Err(error) = pending.persist_noclobber(destination) {
        if error.error.kind() != io::ErrorKind::AlreadyExists {
            return Err(error.error.into());
        }
        if !existing(destination, asset.binary_checksum)? {
            return Err(error.error.into());
        }
    }
    Ok(())
}

pub fn install_plink() -> Result<PathBuf> {
    let asset = Asset::current()?;
    let destination = managed_executable(asset)?;
    if existing(&destination, asset.binary_checksum)? {
        return Ok(destination);
    }
    let agent = Agent::config_builder()
        .https_only(true)
        .timeout_global(Some(Duration::from_secs(60)))
        .tls_config(
            TlsConfig::builder()
                .root_certs(RootCerts::PlatformVerifier)
                .build(),
        )
        .build()
        .new_agent();
    let mut response = agent.get(asset.url).call()?;
    let bytes = read_bounded(response.body_mut().as_reader())?;
    publish(&bytes, asset, &destination)?;
    Ok(destination)
}
