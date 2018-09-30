"""
calcimp.jl
Julia version of calcimp.
"""

include("mensur.jl")
using .Mensur

function showhelp()
    txt = """
        calcimp.jl
        Input impedance calculation program by Julia.

        Available parameters :
            minfreq : default = 0.0
            maxfreq : default = 2000.0
            stepfreq : default = 2.5
            temperature : default = 24.0 (celusius)
            radiation : default = "PIPE". Candidates are PIPE, BAFFLE, NONE.
            output : default = "", using same base name from input file. Stdout is used when "-". 
    """
    println(txt)
end

function parseoptions(args)
    params = Dict([
        "minfreq"=>0.0,
        "maxfreq"=>2000.0,
        "stepfreq"=>2.5,
        "temperature"=>24.0,
        "radiation"=>"PIPE",
        "output"=>"stdout",
        "version"=>false,
        "help"=>false
    ])
    l = length(args)
    for i in 2:l
        tt = split(args[i],":")
        if tt[1] in keys(params)
            params[tt[1]] = tt[2]
        end
    end
    return(params)
end

function main()
    version = v"1.0.0"
    author = "Yoshinobu Ishizaki (ysnbiszk@gmail.com)"

    # command line parser
    fpath = ARGS[1]
    params = parseoptions(ARGS)

    if params["version"] != false
        println("calcimp.jl : ",version)
        println(author)
    elseif params["help"] != false
        showhelp()
    else
        # execute 
        mentable = men_readfile(fpath)
        # print output
    end
    
    # men_printtable(mentable) # debug
end
    
main()
# end 