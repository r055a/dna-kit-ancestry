use crate::{Error, Result};
use serde::{Deserialize, Serialize};
use std::{fmt, num::NonZeroU64, str::FromStr};

macro_rules! quantity {
    (@conversion checked, $name:ident, $inner:ty, $raw:ty, $kind:literal) => {
        impl TryFrom<$raw> for $name {
            type Error = Error;
            fn try_from(val: $raw) -> Result<Self> {
                <$inner>::try_from(val).map(Self).map_err(|_| Error::InvalidValue {
                    kind: $kind,
                    value: val.to_string(),
                })
            }
        }
    };
    (@conversion total, $name:ident, $inner:ty, $raw:ty, $kind:literal) => {
        impl From<$raw> for $name {
            fn from(val: $raw) -> Self { Self(val) }
        }
    };
    ($(#[$doc:meta])* $name:ident, $inner:ty, $raw:ty, $kind:literal, $conversion:ident) => {
        $(#[$doc])*
        #[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
        #[serde(transparent)]
        pub struct $name($inner);

        impl $name {
            pub fn get(self) -> $raw {
                self.0.into()
            }
        }
        quantity!(@conversion $conversion, $name, $inner, $raw, $kind);
        impl From<$name> for $raw {
            fn from(val: $name) -> Self {
                val.get()
            }
        }
        impl FromStr for $name {
            type Err = Error;
            fn from_str(val: &str) -> Result<Self> {
                val.trim().parse::<$inner>().map(Self)
                    .map_err(|_| Error::InvalidValue { kind: $kind, value: val.into() })
            }
        }
        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                self.get().fmt(f)
            }
        }
    };
}

quantity!(
    BasePairLength,
    NonZeroU64,
    u64,
    "positive base-pair length",
    checked
);
quantity!(BasePairCount, u64, u64, "base-pair count", total);
quantity!(
    SnpDensity,
    NonZeroU64,
    u64,
    "positive base pairs per SNP",
    checked
);
quantity!(
    InputByteLimit,
    NonZeroU64,
    u64,
    "positive input byte limit",
    checked
);
quantity!(ReadDepth, u32, u32, "read depth", total);
quantity!(PhredScore, u32, u32, "Phred score", total);

impl SnpDensity {
    pub fn from_kilobases(val: NonZeroU64) -> Result<Self> {
        val.get()
            .checked_mul(1000)
            .ok_or(Error::Overflow)?
            .try_into()
    }
}
