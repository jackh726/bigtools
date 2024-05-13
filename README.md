# Bigtools <a href="https://github.com/jackh726/bigtools"><img align="right" src="https://github.com/jackh726/bigtools/raw/master/assets/bigtools-logo.svg" height="38"></img></a>

[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/jackh726/bigtools/blob/master/LICENSE)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.10606493.svg)](https://doi.org/10.5281/zenodo.10606493)

<table>
    <tr>
        <td>Rust, CLI</td>
        <td>
            <a href="https://crates.io/crates/bigtools">
                <img src="https://img.shields.io/crates/v/bigtools.svg" alt="crates.io">
            </a>
            <a href="https://bioconda.github.io/recipes/bigtools/README.html">
                <img src="https://img.shields.io/conda/vn/bioconda/bigtools.svg?color=green" alt="bigtools on Bioconda">
            </a>
            <a href="https://docs.rs/bigtools">
                <img src="https://img.shields.io/docsrs/bigtools/latest?label=docs.rs" alt="Rust Docs">
            </a>
        </td>
    </tr>
    <tr>
        <td>Python</td>
        <td>
            <a href="https://pypi.org/project/pybigtools/">
                <img src="https://img.shields.io/pypi/v/pybigtools?color=orange" alt="PyPI">
            </a>
            <a href="https://bioconda.github.io/recipes/pybigtools/README.html">
                <img src="https://img.shields.io/conda/vn/bioconda/pybigtools?color=green" alt="pybigtools on Bioconda">
            </a>
            <a href="https://bigtools.readthedocs.io/">
                <img src="https://img.shields.io/readthedocs/bigtools/latest?label=docs" alt="Python Docs">
            </a>
        </td>
    </tr>
</table>

Bigtools is a library and associated tools for reading and writing bigwig and bigbed files.

The primary goals of the project are to be
- Performant
- Extensible
- Modern

### Performant

Bigtools uses `async/await` internally to allow for efficient, multi-core computation when possible. In addition, tools are optimized for minimal memory usage. See [Benchmarks] for more details.

### Extensible

Bigtools is designed to be as modular as possible. This, in addition to the safety and reliability of Rust, allows both flexibility and correctness as a library. In addition, its extremely easy to quickly create new tools or binaries. A number of binaries are available that parallel related existing binaries from [UCSC](https://hgdownload.soe.ucsc.edu/admin/exe/), with drop-in compatibility for the most common flags.

### Modern

Bigtools is written in Rust and published to `crates.io`, so binaries can be installed with `cargo install bigtools` or it can be used as a library by simply including it in your `cargo.toml`.

## Library

See the `bigtools` 🦀 [Documentation](https://docs.rs/bigtools).

### Example

```rust,norun
use bigtools::bbiread::BigWigRead;

let mut reader = BigWigRead::open("test.bigWig").unwrap();
let chr1 = reader.get_interval("chr1", 0, 10000).unwrap();
for interval in chr1 {
    println!("{:?}", interval);
}
```

## Binaries

The following binaries are available:

|binary|description|
| ---- | ----- |
|bigtools|Provides access to multiple subcommands, including all below|
|bedgraphtobigwig|Writes a bigWig from a given bedGraph file|
|bedtobigbed|Writes a bigBed from a given bed file|
|bigbedinfo|Shows info about a provided bigBed|
|bigbedtobed|Writes a bed from the data in a bigBed|
|bigwigaverageoverbed|Calculate statistics over the regions of a bed file using values from a bigWig|
|bigwiginfo|Shows info about a provided bigWig|
|bigwigmerge|Merges multiple bigWigs, outputting to either a new bigWig or a bedGraph|
|bigwigtobedgraph|Writes a bedGraph from the data in a bigWig|
|bigwigvaluesoverbed|Get the per-base values from a bigWig over the regions of a bed file using values|

Renaming the `bigtools` binary to any of the subcommands (case-insensitive) allows you to run that subcommand directly.

The `bigtools` CLI binaries can be installed through `crates.io` (`cargo install bigtools`) or [`conda`](https://anaconda.org/bioconda/bigtools/). Additionally, pre-built binaries can be downloaded through [Github releases](https://github.com/jackh726/bigtools/releases).

## Python wrapper

Also included in this repo is a Python wrapper, `pybigtools` written using [`PyO3`](https://pyo3.rs/).

See the `pybigtools` 🐍 [API Documentation](https://bigtools.readthedocs.io/en/latest).

The `pybigtools` package can be used as a dpendency either through [pypi](https://pypi.org/project/pybigtools/) or [conda](https://anaconda.org/bioconda/pybigtools/).

## How to build from source

In order to build the bigtools binaries, you can run

```
cargo build --release
```

and the binaries can be found in `target/release/`.

Otherwise, you can install the binaries from source by running

```
cargo install --path bigtools/
```

Building the python wheels for pybigtools requires [maturin](https://pypi.org/project/maturin/). To build the pybigtools wheel for installation (and install), you can run

```
maturin build --release -m pybigtools/Cargo.toml
pip install target/wheels/pybigtools*.whl
```

or

```
maturin develop --release -m pybigtools/Cargo.toml
```

## Benchmarks
[Benchmarks]: #Benchmarks

Benchmarks are included in the `./bench` directory. They require `python` to run.

Multiple tools are compared against the comparable UCSC tools. For completeness, both single-threaded and multi-threaded (when available) benchmarks are included. Multiple different configuration options are benchmarked across multiple replicates, but a summar is available in the table below:

<img src="https://github.com/jackh726/bigtools/raw/master/assets/bigtools-bench.png"></img></a>

## How to cite

This repository contains contains a `CITATION.cff` file with citation information. Github allows you to get a citation in either APA or BibTeX format; this is available in "Cite this repository" under About.
