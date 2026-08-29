# DNA Ancestry Estimate

Compute F_ROH from raw genetic DNA data kit file using [PLINK](https://github.com/chrchang/plink-ng).

F_ROH estimates provide an approximate representation of genomic inbreeding or shared parental ancestry:
 
* ~25.000% - Parent-child or siblings (example: King Charles II of Spain)
* ~12.500% - Half siblings, grandparent-grandchild, uncle–niece or aunt–nephew
* ~6.2500% - First cousins, half uncle–niece or half aunt–nephew
* ~3.1250% - First cousins once removed
* ~1.5625% - Second cousins

Refer [here](https://en.wikipedia.org/wiki/Coefficient_of_inbreeding#Iterated_sibling_mating) for coefficients of multigenerational inbreeding.

## Requirements

> A raw genetic DNA data kit file used for genealogical research.

## Build

```sh
cargo build
```

## Format

```sh
cargo fmt --all
cargo clippy --all-targets
```

## Install

```sh
cargo run -- --install-plink
```

## Run

```sh
cargo run -- [OPTIONS]
```
