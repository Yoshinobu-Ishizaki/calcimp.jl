"""
calcimp.jl
Julia version of calcimp.
"""

using DataFrames

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
            output : default = "", using same base name from input file. Stdout is used when "stdout". 

        Usage:
            julia calcimp.jl file.xmen
            julia calcimp.jl file.xmen option:arg...
    """
    println(txt)
end

function parse_opt(args)
    params = Dict([
        "minfreq"=>0.0,
        "maxfreq"=>2000.0,
        "stepfreq"=>2.5,
        "temperature"=>24.0,
        "radiation"=>"PIPE",
        "output"=>"",
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

function print_params(param)
    for (k,v) in pairs(param)
        println(k,":",v)
    end
end

function set_output!(fpath, params)
    if length(params["output"]) == 0
        bdy, ext = splitext(fpath)
        fout = bdy * ".imp"
        params["output"] = fout
    end
end

function main()
    version = v"1.0.0"
    author = "Yoshinobu Ishizaki (ysnbiszk@gmail.com)"

    # command line parser
    fpath = ARGS[1]
    params = parse_opt(ARGS)

    set_output!(fpath,params)
    # print_params(params)

    if params["version"] != false
        println("calcimp.jl : ",version)
        println(author)
    elseif params["help"] != false
        showhelp()
    else
        mentable = men_readfile(fpath)
        # execute 
        imped = input_impedance(mentable,params)
        # print output
        print_impedance(imped)
    end
    
    men_printtable(mentable) # debug
end
    
main()
# end 