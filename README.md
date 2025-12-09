# ATLabData

Julia package for loading and handling data from the 
[ATLab](https://github.com/turbulencia/atlab) solver.


## Loading data
Loading is done by _load_, which consists of methods for loading data from visuals.x and raw dns.x output:
```
visuals_data = load("path/to/dir/Buoyancy000000")
raw_data = load("path/to/dir/scal.000000.1")
```

Data of 3 different files can also be loaded into a single vector-like data structure:
```
vector_data = load(
    "path/to/dir/flow.000000.1",
    "path/to/dir/flow.000000.2",
    "path/to/dir/flow.000000.3",
)
```

The _grid_ file (inigrid.x) has to be present in the same directory.


## Data structure

### Grid
Composite type containing the grid information.  
Attributes are:  
- nx (grid points)
- nx (grid points)
- nz (grid points)
- x (x-axis vector)
- y (y-axis vector)
- z (z-axis vector)
- scalex (x-axis length)
- scaley (y-axis length)
- scalez (z-axis length)


### ScalarData
Composite type returned by _load_ with in single String argument.

### VectorData
Composite type returned by _load_ with three String arguments.


## Handling data
Standard operations are overloaded for _ScalarData_ and _VectorData_:

```
data1 = load("file1")
data2 = load("file2")

data = data1 + data2
data = data1 - data2
data = data1 / data2
data = data1 * data2
```

## Statistics
Basic statistical operations are implemented returning _ScalarData_:
```
data = load("/path/to/file")
meandata = mean(data)
flucsdata = flucs(data)
```