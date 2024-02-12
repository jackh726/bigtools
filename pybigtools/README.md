The `pybigtools` python package wraps the `bigtools` Rust library and allows effecient reading and writing of bigWig and bigBed files.

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
