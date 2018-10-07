# calcimp.jl
Calcimp julia version

# Description
Calcimp.jl is a Julia version of input impedance calculation program for wind instruments.

It reads .xmen file and return (frq,imp) dataframe.

# Requirements

[DataFrames.jl](https://github.com/JuliaData/DataFrames.jl), [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl), [Struve.jl](https://github.com/gwater/Struve.jl)

# Basic usage

```julia
include("calcimp.jl")
using .Calcimp

calcimp("file.xmen"; params...)

# params can be created using ;
prm = initcalcparams(fmin=0.0, fmax=2000, fstep=2.5, tmpr=24.0, "PIPE")
# arguments can be omitted when using default values.
# radiation types are "PIPE/BAFFLE/NONE"

```

# Documents

1. [calcimp basics](./doc/calcimp_basics.md)

# Examples

See [example](./example) folder.

# ChangeLog

version | date | description
----|------|------------
0.3.0 | 08.10.2018 | Added magnitude column to output.
0.2.0 | 07.10.2018 | Tune up. Runs almost as fast as C.
0.1.0 | 06.10.2018 | First completion 


