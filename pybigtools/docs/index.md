# pybigtools

Python bindings to the [Bigtools](https://github.com/jackh726/bigtools) Rust
library for high-performance reading and writing of **BigWig** and **BigBed**
files.

```{toctree}
:maxdepth: 2
:hidden:

api
```

## Installation

```sh
pip install pybigtools
```

## Quickstart

Open a file for reading. {func}`~pybigtools.open` auto-detects whether the file
is a BigWig or BigBed and returns a {class}`~pybigtools.BBIReader`:

```python
import pybigtools

b = pybigtools.open("path/to/file.bigWig")  # also accepts an http(s) URL
print(b.chroms())          # {'chr1': 248956422, ...}
print(b.info())            # version, summary stats, zoom levels, ...
```

Rasterize values over a region as a NumPy array:

```python
values = b.values("chr1", 0, 1000)        # one value per base
binned = b.values("chr1", 0, 1_000_000, bins=1000, summary="mean")
```

Iterate over raw records (intervals for BigWig, BED entries for BigBed):

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

### Writing

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
