To build and install in a local virtualenv run
```
maturin build --release
pip install -I target/wheels/pybigtools-0.1.0-cp311-cp311-manylinux_2_28_x86_64.whl
```

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

## Open a file-like object

```python
import pybigtools
b = pybigtools.open(open(<path>, 'rb'))
a = b.values("chr1")
print(a.shape)
```

## Documenation
To generate documenation install pdoc3 with `pip install pdoc3`.
Then make sure you have pybigtools install as above.
Then run `pdoc3 --html --force pybigtools`.