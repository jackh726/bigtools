To build and install in a local virtualenv run
```
maturin build --release
pip install -I target/wheels/pybigtools-0.1.0-cp311-cp311-manylinux_2_28_x86_64.whl
```


## Documenation

To generate documenation install pdoc3 with `pip install pdoc3`.
Then make sure you have pybigtools install as above.
Then run `pdoc3 --html --force pybigtools`.