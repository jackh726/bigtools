The `pybigtools` python package wraps the `bigtools` Rust library and allows effecient reading and writing of bigWig and bigBed files.

See the [documentation](https://bigtools.readthedocs.io) for more information.

## Development

Pybigtools is a mixed Python-Rust project using the pyo3 `maturin` build system to create a Python package. We recommend using `uv` to manage your Python development environment.

To set up your environment, from the top-level pybigtools directory:
```sh
uv sync
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

To build the documentation with hot reloading in the browser:

```sh
cd pybigtools
uv run --group docs sphinx-autobuild docs docs/_build/html/
```