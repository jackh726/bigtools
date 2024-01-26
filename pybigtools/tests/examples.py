
# While the built-in `bigWigAverageOverBed` uses compiled code for speed,
# you could also implement as so
def bigWigAverageOverBed(bigWig, bed, names = None):
    """
    Gets the average values from a bigWig over the entries of a bed file.

    Parameters:
        bigWig (str): The path to the bigWig.
        bed (str): The path to the bed.
        names (None, bool, or int):  
            If `None`, then each return value will be a single `float`,
                the average value over an interval in the bed file.  
            If `True`, then each return value will be a tuple of the value of column 4
                and the average value over the interval with that name in the bed file.  
            If `False`, then each return value will be a tuple of the interval in the format
                `{chrom}:{start}-{end}` and the average value over that interval.  
            If `0`, then each return value will match as if `False` was passed.  
            If a `1+`, then each return value will be a tuple of the value of column of this
                parameter (1-based) and the average value over the interval.  

    Returns:
        This returns a generator of values. (Therefore, to save to a list, do `list(bigWigAverageOverBed(...))`)  
        If no name is specified (see the `names` parameter above), then returns a generator of `float`s.  
        If a name column is specified (see above), then returns a generator of tuples `({name}, {average})`  
    """
    if names is None:
        namecol = None
    elif names is False:
        namecol = 0
    elif names is True:
        namecol = 4
    elif isinstance(names, int):
        if names < 0:
            raise Exception("Invalid names argument. Must be >= 0.")
        namecol = names
    else:
        raise Exception("Invalid names argument. Should be either `None`, a `bool`, or an `int`")

    from pybigtools import openbig
    from math import isnan

    bigWigFile = openbig(bigWig)
    with open(bed) as bed_file:
        for line in bed_file:
            split = line.strip().split()
            chrom = split[0]
            start = int(split[1])
            end = int(split[2])
            if namecol is None:
                name = None
            else:
                if len(split) <= namecol:
                    raise Exception(f"Specified name column ({namecol + 1}) exceeds number of columns {len(split)}")
                if namecol == 0:
                    name = f"{chrom}:{start}-{end}"
                else:
                    name = split[namecol - 1]
            length = end - start
            vals_iter = bigWigFile.values(chrom, start, end)
            validonly = True
            if validonly:
                vals = [v for v in vals_iter if not isnan(v)]
            else:
                vals = [v if not isnan(v) else 0 for v in vals_iter]
            avg = (sum(vals) / len(vals))
            if namecol is not None:
                yield (name, avg)
            else:
                yield avg
