pub mod analysis;
pub mod import;
pub mod install;
pub mod model;
mod process;
pub mod units;

use std::path::PathBuf;
use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    #[error("{0}")]
    InvalidInput(String),
    #[error("invalid {kind}: {value:?}")]
    InvalidValue { kind: &'static str, value: String },
    #[error("{path}: {source}")]
    File {
        path: PathBuf,
        source: std::io::Error,
    },
    #[error("{format} is recognized but needs {requirement}")]
    NeedsConversion {
        format: &'static str,
        requirement: &'static str,
    },
    #[error("choose one {kind} with {option}; available: {choices}")]
    Ambiguous {
        kind: &'static str,
        option: &'static str,
        choices: String,
    },
    #[error("{tool} failed: {detail}")]
    Tool { tool: String, detail: String },
    #[error("executable not found: {tool}")]
    ToolNotFound { tool: String },
    #[error(
        "PLINK executable not found: {executable}. Run `dna-kit --install-plink`, or select an installed PLINK 1.9 with `--plink PATH`"
    )]
    PlinkNotFound { executable: PathBuf },
    #[error("no usable autosomal diploid SNP calls; inspect the import statistics")]
    NoUsableCalls,
    #[error("numeric range exceeded")]
    Overflow,
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Csv(#[from] csv::Error),
    #[error(transparent)]
    Json(#[from] serde_json::Error),
    #[error(transparent)]
    Zip(#[from] zip::result::ZipError),
    #[error(transparent)]
    Install(#[from] install::InstallError),
}
