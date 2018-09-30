"""
calcimp.jl
Julia version of calcimp.
"""

version = v"0.1.0"

include("mensur.jl")
using .Mensur

# command line parser
fpath = ARGS[1]
# @show fpath
# execute 
mentable = men_readfile(fpath)
# print output

men_printtable(mentable)

# end 