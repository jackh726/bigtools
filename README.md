# Bigtools <a href="https://github.com/jackh726/bigtools"><img align="right" src="https://github.com/jackh726/bigtools/raw/master/assets/bigtools-logo.svg" height="38"></img></a>

[![crates.io](https://img.shields.io/crates/v/bigtools.svg)](https://crates.io/crates/bigtools)
[![PyPI](https://img.shields.io/pypi/v/pybigtools?color=green)](https://pypi.org/project/pybigtools/)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/jackh726/bigtools/blob/master/LICENSE)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.10606493.svg)](https://doi.org/10.5281/zenodo.10606493)

[![Rust Docs](https://img.shields.io/docsrs/bigtools/latest?label=docs.rs)](https://docs.rs/bigtools)
[![Python Docs](https://img.shields.io/readthedocs/bigtools/latest?label=docs|rtd)](https://bigtools.readthedocs.io/)

Bigtools is a library and associated tools for reading and writing bigwig and bigbed files.

The primary goals of the project are to be
- Performant
- Extensible
- Modern

### Performant

Bigtools uses `async/await` internally to allow for efficient, multi-core computation when possible. In addition, tools are optimized for minimal memory usage. See [Benchmarks] for more details.

### Extensible

Bigtools is designed to be as modular as possible. This, in addition to the safety and reliability of Rust, allows both flexibility and correctness as a library. In addition, its extremely easy to quickly create new tools or binaries.. (TODO: mention python wrapper)

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

TODO

## Python wrapper

See the `pybigtools` 🐍 [API Documentation](https://bigtools.readthedocs.io/en/latest).

## Benchmarks
[Benchmarks]: #Benchmarks

Benchmarks are included in the `./bench` directory. The require `python` to run.

Multiple tools are compared against the comparable UCSC tools. For completeness, both single-threaded and multi-threaded (when available) benchmarks are included.

TODO: include image of benchmarks
