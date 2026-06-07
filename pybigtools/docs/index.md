# Bigtools

```{toctree}
:maxdepth: 2
:hidden:

api
Rust API Reference <https://docs.rs/bigtools>
```

Bigtools is a modern, high-performance library and associated tools for reading and writing **BigWig** and **BigBed** files. Bigtools is written in Rust and bindings for Python are available.

## CLI binaries

Drop-in replacements for the UCSC binaries. Install with `cargo` or via bioconda (conda, mamba, pixi, etc.).

```sh
cargo install bigtools
```

```sh
conda install -c bioconda bigtools
```

## Python package

Python bindings to the bigtools Rust library is provided by the [pybigtools](https://pypi.org/project/pybigtools/) package. Install with `pip`, `uv` or via bioconda (conda, mamba, pixi, etc.).


```sh
pip install pybigtools
```

```sh
conda install -c bioconda pybigtools
```

To open a file for reading, {func}`~pybigtools.open` auto-detects whether the file is a BigWig or BigBed and returns a {class}`~pybigtools.BBIReader`:

```python
import pybigtools

b = pybigtools.open("path/to/file.bigWig")  # also accepts an http(s) URL
print(b.chroms())          # {'chr1': 248956422, ...}
print(b.info())            # version, summary stats, zoom levels, ...
```

To rasterize values over a region as a NumPy array:

```python
values = b.values("chr1", 0, 1000)        # one value per base
binned = b.values("chr1", 0, 1_000_000, bins=1000, summary="mean")
```

To iterate over raw records (intervals for BigWig, BED entries for BigBed):

```python
for start, end, value in b.records("chr1"):
    ...
```

Files can be used as context managers, and file-like objects are accepted in
place of a path:

```python
with pybigtools.open(open("path/to/file.bigBed", "rb")) as b:
    schema = b.sql(parse=True)
```

For writing:

```python
import pybigtools

w = pybigtools.open("out.bigWig", "w")
w.write(
    {"chr1": 248956422},
    [("chr1", 0, 100, 1.5), ("chr1", 100, 200, 2.0)],
)
# the file is closed automatically once write() completes
```

See the [API reference](api.rst) for the full set of methods and options.


## Rust crate

To use [bigtools](https://crates.io/crates/bigtools) in your Rust project, add bigtools to your Cargo.toml or run:

```sh
cargo add bigtools
```

See the Rust documentation [here](https://docs.rs/bigtools).
