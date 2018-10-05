__precompile__()

"""
    calcimp.jl

Input impedance calculation program by Julia.
Returns [freq::Float, impedance::Complex] Dataframe.

Available parameters :
- minfreq : default = 0.0
- maxfreq : default = 2000.0
- stepfreq : default = 2.5
- temperature : default = 24.0 (celusius)
- radiation : default = "PIPE". Candidates are PIPE, BAFFLE, NONE.

# Example:

```julia
calcimp(file.xmen) # with default params
```

```julia
calcimp(file.xmen, params) # with custom params
```
"""
module Calcimp

include("mensur.jl")

# export structs, functions
export Men
export calcimp, initcalcparam

end 