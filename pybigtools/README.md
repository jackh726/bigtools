The `pybigtools` python package wraps the `bigtools` Rust library and allows effecient reading and writing of bigWig and bigBed files.

## Development

Pybigtools is a mixed Python-Rust project using the pyo3 `maturin` build system to create a Python package. We recommend using `uv` for manage your Python development environment.

To set up your environment, from the top-level pybigtools directory:
```sh
uv sync
```

Building:
```sh
maturin develop --release --uv
```

Testing:
```sh
cargo test  # Rust code
uv run pytest  # Python code
```

Linting:
```sh
cargo fmt --all -- --check  # Rust code
uv run ruff check  # Python code
```

Formatting:
```sh
cargo fmt  # Rust code
uv run ruff format  # Python code
```

## Documenation

Documentation is available on [readthedocs](https://bigtools.readthedocs.io/en/latest/pybigtools.html).

## Examples

(Replace `<path>` with the path to a bigWig file or url)

### Iterator of intervals

```python
import pybigtools
b = pybigtools.open(<path>)
i = b.intervals("chr1")
print(next(i))
```

### Numpy array of values

```python
import pybigtools
b = pybigtools.open(<path>)
a = b.values("chr1")
print(a.shape)
```

### Open a file-like object

```python
import pybigtools
b = pybigtools.open(open(<path>, 'rb'))
a = b.values("chr1")
print(a.shape)
```
