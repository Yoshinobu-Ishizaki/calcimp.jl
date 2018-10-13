using DataFrames
using SpecialFunctions
using LinearAlgebra
using Struve

version = v"0.4.0"

# constants
const OPEN = 1
const CLOSE = 0
const HALF = 0.5
const HEAD = 0
const LAST = 1

const GMM = 1.4  # specific head ratio
const PR = 0.72  # Prandtl number

const type_keywords = Dict("SPLIT"=>:split, 
    "TONEHOLE"=>:split,
    "|"=>:split,
    "VALVE_OUT"=>:branch,
    "BRANCH"=>:branch,
    "<"=>:branch,
    "MERGE"=>:merge,
    "VALVE_IN"=>:merge,
    ">"=>:merge,
    "OPEN_END"=>:open,
    "CLOSED_END"=>:close,
    "INSERT"=>:insert,
    "@"=>:insert
    )

const group_keywords = Dict("MAIN"=>:main,
    "["=>:main,
    "END_MAIN"=>:endmain,
    "]"=>:endmain,
    "GROUP"=>:group,
    "{"=>:group,
    "END_GROUP"=>:endgroup,
    "}"=>:endgroup
    )

type_keys = keys(type_keywords)
group_keys = keys(group_keywords)

"""Calculate dBSPL of pressure x"""
dBSPL(x) = 20log10(abs(x/2.0e-5))

"""
Men struct which contains basic unit for expressing mensur structure.
"""
mutable struct Men
    # data to store
    df::Float64 # diameter at input 
    db::Float64 # diameter at output
    r::Float64  # length of mensur cell
    cm::AbstractString # comment
    gr::AbstractString # group name
    
    # calculated values
    zi::ComplexF64 # impedance at input 
    zo::ComplexF64 # impedance at output
    vi::Array{ComplexF64,1} # (p,U) vector at input
    vo::Array{ComplexF64,1} # (p,U) vector at output
    tm::Array{ComplexF64,2} # transmission matrix
    
    # child info
    cname::String
    ctype::String
    cratio::Float64

    # list reference
    prev::Union{Men,Nothing}
    next::Union{Men,Nothing}
    child::Union{Men,Nothing}
    parent::Union{Men,Nothing}
    
    # incomplete constructor
    Men() = new()
end

function isbranch(t::String)
    return(type_keywords[t] == :branch)
end

function issplit(t::String)
    return(type_keywords[t] == :split)
end

function isaddon(t::String)
    return(type_keywords[t] == :addon)
end

function ismerge(t::String)
    return(type_keywords[t] == :merge)
end

"""
create new mensur 
"""
function men_create(df::Number=0.0, db::Number=0.0, r::Number=0.0, cm::AbstractString="", gr::AbstractString="")
    men = Men()

    men.df = df
    men.db = db
    men.r = r
    men.cm = cm
    men.gr = gr

    men.zi = complex(0.0)
    men.zo = complex(0.0)
    men.vi=zeros(ComplexF64,2) 
    men.vo=zeros(ComplexF64,2)
    men.tm=zeros(ComplexF64,2,2)

    men.cname = ""
    men.ctype = ""
    men.cratio = 0.0

    men.prev = nothing
    men.next = nothing
    men.parent = nothing
    men.child = nothing

    return(men)
end

function men_append!(prev::Men, curr::Men)
    prev.next = curr
    curr.prev = prev
    return(curr) 
end

"""insert ins after cur"""
function men_insert!(cur::Men, ins::Men)
    nxt = cur.next
    inse = men_end(ins)

    cur.next = ins
    ins.prev = cur

    nxt.prev = inse
    inse.next = nxt
end

function men_setchildinfo!(men::Men, t::String, n::String,r::Float64)
    men.ctype=t
    men.cname=n
    men.cratio=r
end

function men_top(curr::Men)
    c = curr
    if c != nothing
        while c.prev != nothing
            c = c.prev
        end
    end
    return(c)
end

function men_end(curr::Men)
    c = curr
    if c != nothing
        while c.next != nothing
            c = c.next
        end
    end
    return(c)
end

function eat_comment(s)
    m = match(r"#", s)
    if m != nothing
        ss = s[1:m.offset-1]
    else
        ss = s
    end
    return(ss)
end

function getfbr(cur::Men)
    return(cur.df, cur.db,cur.r,cur.cm)
end

function men_bykeyword(cur::Men, wd::Array{String,1})
    w = wd[1]
    if length(wd) > 1
        name = wd[2]
    end
    if length(wd) > 2
        ratio = float(eval(Meta.parse(wd[3])))
    end
    df,db,r,c = getfbr(cur)
    gr = cur.gr

    men = nothing
    if type_keywords[w] == :close
        men = men_create(0.0,0.0,0.0,"",gr )
    elseif type_keywords[w] == :open
        men = men_create(db,0.0,0.0,"",gr )
    else
        if type_keywords[w] == :insert
            ratio = 1.0
        end
        men = men_create(db,db,0.0,"",gr )
        men_setchildinfo!(men,w,name,ratio)
    end
    return(men)
end

function men_setchild!(cur::Men, sid::Men)
    cur.child = sid
    sid.parent = cur
end

function men_resolvechild!(mentable)
    cur = mentable["MAIN"]

    while cur != nothing
        cnm = cur.cname
        typ = cur.ctype
        if cnm != ""
            if type_keywords[typ] == :insert
                sid = mentable[cnm]
                men_insert!(cur,sid)
            else
                if type_keywords[typ] == :merge
                    sid = men_end(mentable[cnm])
                else
                    sid = mentable[cnm]
                end
                men_setchild!(cur,sid)
            end
        end
        cur = cur.next
    end
end

"""Build mensur from text lines""" 
function men_build(lns)
    mentable = Dict{String,Men}()
    grouplist = []
    groupnames = []
    cur = nothing

    for ln in lns
        # @show ln
        ln = eat_comment(ln)
        ln = replace(ln,r"[\t\s]+"=>"") # remove all white spaces
        
        if ln != ""
            if occursin("=",ln)
                # evaluate var definition line
                eval(Meta.parse(ln))
            else
                # normal df,db,r,comment line
                wd = String.(split(ln,","))

                w = wd[1] # df,db,r or group definition
                if w in group_keys
                    # group definitions
                    if group_keywords[w] == :endmain
                        grouplist = []
                        cur = nothing
                    elseif group_keywords[w] == :endgroup
                        pop!(grouplist)
                        if length(grouplist) == 0
                            cur = nothing
                        end
                    elseif group_keywords[w] == :main
                        "MAIN" in groupnames ? throw(AssertionError("Duplicating MAIN.")) : nothing
                        push!(grouplist,"MAIN")
                        push!(groupnames,"MAIN")
                    elseif group_keywords[w] == :group
                        push!(grouplist,wd[2])
                        gname = join(grouplist,":")
                        gname in groupnames ? throw(AssertionError("group name $group_name is not unique.")) : nothing
                        push!(groupnames,gname)
                    end
                elseif w in type_keys
                    # BRANCH, MERGE, SPLIT, etc
                    men = men_bykeyword(cur,wd)
                    cur = men_append!(cur,men)
                elseif length(wd) > 2
                    # normal "df,db,r,cmt" like string
                    df = eval(Meta.parse(wd[1]))
                    db = eval(Meta.parse(wd[2]))
                    r = eval(Meta.parse(wd[3]))
                    cmt = length(wd) > 3 ? wd[4] : ""

                    men = men_create(df*0.001, db*0.001,r*0.001, cmt, groupnames[end])

                    # set mensur group dictionary
                    if cur == nothing
                        cur = men
                        if length(mentable) == 0
                            push!(mentable,"MAIN"=>cur)
                        else
                            push!(mentable,groupnames[end] => cur)
                        end
                    else
                        cur = men_append!(cur,men)
                    end
                else    
                    throw(AssertionError("string $ln is invalid."))
                end
            end
        end
    end
    # men_printtable(mentable) # debug
    men_resolvechild!(mentable)
    return(mentable) 
end

function readxmen(path)
    lns = readlines(path)

    mentable = []
    if length(lns) > 0 
        mentable = men_build(lns)    
    end
    return(mentable)
end

function men_print(men::Men)
    gr = men.group
    if gr != "MAIN"
        println("GROUP,",gr)
    end

    while men != nothing
        df,db,r,c = getfbr(men)
        if men.child != nothing
            println(men.ctype,",",men.cname,",",men.cratio)
        else
            println(df*1000.0,",",db*1000.0,",",r*1000.0,",",c)
        end
        men = men.next
    end
    if gr != "MAIN"
        println("END_GROUP")
    end
end

function men_printtable(mentable)
    for (key, men) in pairs(mentable)
        if key == "MAIN"
            println(key)
            men_print(men)
            println("END_MAIN")
        else
            men_print(men)
        end
    end
end

function joint_mensur(men::Men)
    if isbranch(men.ctype)
        e = men_end(men.child)
        if e.parent != nothing && ismerge(e.parent.ctype)
            return(e.parent)
        else
            return(nothing)
        end
    else
        return(nothing)
    end
end

"""
create calculation parameters

Example:
    initcalcparam(minfreq,maxfreq,stepfreq,temperature,radiation)
    initcalcparam(0.0,2000.0,2.5,"PIPE/BAFFLE/NONE")
"""
function initcalcparam(;minfreq=0.0, maxfreq=2000.0, stepfreq=2.5, temperature=24.0, radiation="PIPE")
    c0 = 331.45 * sqrt(temperature / 273.16 + 1)
    rho = 1.2929 * (273.16 / (273.16 + temperature))
    rhoc0 = rho * c0
    mu = (18.2 + 0.0456*(temperature - 25)) * 1.0e-6  # viscosity constant. Linear approximation from Scientific Dictionary.
    nu = mu/rho  # dynamic viscous constant.
    dmp = (1+(GMM-1)/sqrt(PR)) # wall dumping factor
    
    Dict(
        :minfreq=>minfreq,
        :maxfreq=>maxfreq,
        :stepfreq=>stepfreq,
        :temperature=>temperature,
        :radiation=>radiation,
        :c0=>c0,
        :rho =>rho,
        :rhoc0 => rhoc0,
        :mu => mu,
        :nu => nu,
        :dmp => dmp
    )
end

function radimp(wf::Float64, dia::Float64; c0::Float64, rhoc0::Float64, radiation::String,kwd...)
    if wf > 0
        if dia == 0.0
            zr = complex(Inf) # closed end
        else
            if radiation == "NONE"
                zr = complex(0)  # simple open end impedance
            else
                s = dia^2*pi/4.0
                k = wf/c0
                x = k*dia

                re = rhoc0/s*(1 - besselj1(x)/x*2)  # 1st order bessel function.
                img = rhoc0/s*Struve.H(1,x)/x*2  # 1st order struve function.

                if radiation == "BAFFLE"
                    zr = complex(re,img)
                elseif radiation == "PIPE"
                    # real is about 0.5 times and imaginary is 0.7 times when without frange.
                    zr = complex(0.5*re,0.7*img)
                end
            end
        end
        return(zr) 
    end
end

function calc_transmission(wf::Float64,men::Men; c0::Float64,rhoc0::Float64,nu::Float64,dmp::Float64,kwd...)
    df = men.df
    db = men.db
    r = men.r
    d = (df + db)*0.5
    aa = dmp * sqrt(2*wf*nu)/c0/d  # wall dumping factor
    k = sqrt((wf/c0)*(wf/c0 - 2*complex(-1,1)*aa))  # complex wave number including wall dumping
    x = k * r
    cc = cos(x)
    ss = sin(x)

    tm = zeros(ComplexF64,2,2)

    if df != db
        # taper
        r1 = df*0.5
        r2 = db*0.5
        dr = r2-r1

        tm[1,1] = (r2*x*cc - dr*ss)/(r1*x)
        tm[1,2] = im*rhoc0*ss/(pi*r1*r2)
        tm[2,1] = -im*pi*(dr*dr*x*cc - (dr*dr + x*x*r1*r2)*ss)/(x*x*rhoc0)
        tm[2,2] = (r1*x*cc + dr*ss)/(r2*x)
    else
        # straight
        s1 = pi/4*df*df
        tm[1,1] = tm[2,2] = cc
        tm[1,2] = im*rhoc0*ss/s1
        tm[2,1] = im*s1*ss/rhoc0
    end

    return(tm)
end

function zo2zi(tm::Array{ComplexF64,2},zo::ComplexF64)
    if !isinf(zo)
        zi = (tm[1,1]*zo + tm[1,2])/(tm[2,1]*zo + tm[2,2])
    else
        if tm[2,1] != 0
            zi = tm[1,1]/tm[2,1]
        else
            zi = complex(Inf)
        end
    end
    return(zi)
end

function transmission_matrix(men1::Men, men2::Union{Men,Nothing})
    if men2 == nothing
        men = men_end(men1)
        men = men.prev
    else
        men = men2
    end

    m = Matrix{ComplexF64}(I,2,2) # eye 
    while men != nothing && men != men1
        m *= men.tm
        men = men.prev
    end
    m *= men1.tm
    
    return(m)
end

function child_impedance(wf::Float64,men::Men; params...)
    if issplit(men.ctype)
        # split (tonehole) type.
        # @show men.child
        impedance!(wf,men.child; params...)  # recursive call for input impedance
        if men.cratio == 0
            men.zo = men.next.zi
        else
            z1 = men.child.zi / men.cratio  # adjust blending ratio
            z2 = men.next.zi
            if z1 == complex(0.0) && z2 == complex(0.0)
                z = complex(0.0)
            else
                z = z1*z2/(z1+z2)
            end
            men.zo = z
        end
    elseif isbranch(men.ctype) && men.cratio > 0
        # multiple tube connection
        impedance!(wf, men.child;params...)
        m = transmission_matrix(men.child, nothing)
        jnt = joint_mensur(men)
        n = transmission_matrix(men.next, jnt)

        # section area adjustment
        if men.cratio == 1
            m[1,2] = complex(Inf)
        else
            m[1,2] /= (1 - men.cratio)
        end
        m[2,1] *= (1 - men.cratio)

        if men.cratio == 0
            n[1,2] = complex(Inf)
        else
            n[1,2] /= men.cratio
        end
        n[2,1] *= men.cratio

        z2 = jnt.next.zi
        dv = (m[2,2]*n[1,2] + m[1,2]*n[2,2] + (
            (m[1,2] + n[1,2])*(m[2,1] + n[2,1]) - (m[1,1] - n[1,1])*(m[2,2] - n[2,2]))*z2)
        if dv != complex(0.0)
            z = (m[1,2]*n[1,2] + (m[1,2]*n[1,1] + m[1,1]*n[1,2])*z2)/dv
        else
            z = 0
        end
        men.zo = z
    elseif isaddon(men.ctype) && men.cratio > 0
        # this routine will not called until 'ADDON(LOOP)' type of connection is implemented.
        impedance!(wf, men.child,params)
        m = transmission_matrix(men.child, nothing)
        z1 = m[1,2]/(m[1,2]*m[2,1]-(1-m[1,1])*(1-m[2,2]))
        z2 = men.next.zi
        if men.cratio == 0
            men.zo = men.next.zi
        elseif men.cratio == 1
            men.zo = z1
        else
            z1 /= men.cratio
            z2 /= (1 - men.cratio)
            if z1 == complex(0.0) && z2 == complex(0.0)
                z = complex(0.0)
            else
                z = z1*z2/(z1+z2)
            end
            men.zo = z
        end
    end
end

"""recursively called function to calculate impedance at each mensur node"""
function calc_impedance!(wf::Float64,men::Men; params...)
    if men.child != nothing
        child_impedance(wf, men; params...)
    elseif men.next != nothing
        men.zo = men.next.zi
    end

    if men.r > 0.0
        men.tm = calc_transmission(wf,men;params...)
        men.zi = zo2zi(men.tm, men.zo)
    else
        # length 0
        men.zi = men.zo
    end
end

"""calculate impedance at given frequency """
function impedance!(wf::Float64, men::Men; params...)
    if wf > 0
        cur = men_end(men)
        # end impedance
        cur.zo = radimp(wf,cur.df; params...)

        while cur != men
            calc_impedance!(wf,cur; params...)
            cur = cur.prev
        end
        calc_impedance!(wf,men; params...)
        return(men.zi)
    else
        return(complex(0.0))
    end
end

function input_impedance(mentable;params...)
    minfreq = params[:minfreq]
    maxfreq = params[:maxfreq]
    stepfreq = params[:stepfreq]

    imped = DataFrame()

    frq = range(minfreq, stop=maxfreq, step=stepfreq)
    imped[:frq] = frq
    imped[:imp] = zeros(ComplexF64,length(frq))

    men = mentable["MAIN"]

    ss = men.df^2 * pi/4

    for i in 1:length(frq)
        wf = frq[i]*2*pi
        imped[:imp][i] = ss * impedance!(wf,men;params...)
    end

    return(imped)
end

"Calculate input impedance for given xmensur file."
function calcimp(fpath::String; params...)
    mentable = readxmen(fpath)
    prm = initcalcparam(;params...)
    imped = input_impedance(mentable;prm...)
    imped[:mag] = 20log10.(abs.(imped[:imp]))
    return(imped)
end

"Reuse read data"
function calcimp(mtable::Dict{String,Men};params...)
    prm = initcalcparam(;params...)
    imped = input_impedance(mtable;prm...)
    imped[:mag] = 20log10.(abs.(imped[:imp]))
    return(imped)
end

function actualnext(men::Men)
    if men.ctype == "BRANCH" && men.cratio > HALF
        m = men.child
    elseif men.next != nothing
        m = men.next
    else
        m = men.parent 
    end
    m 
end

function slice_mensur!(men::Men, step::Float64)
    if step > 0
        while men != nothing
            if men.r > step
                if men.db != men.df
                    db2 = (men.db - men.df)/men.r * step + men.df
                else
                    db2 = men.db
                end
                r2 = men.r - step
                nw = men_create(db2, men.db, r2, men.cm, men.gr)
                men_insert!(men, nw)
                # change self data
                men.r = step
                men.db = db2
            end
            men = actualnext(men)
        end
    end
end

function calc_pressure(p::Float64,men::Men)
    v = complex([p,0.0]) # initial pressure, volume speed pair
    x = 0.0 # length from top
    prs = DataFrame(x=[],D=[],p=[],U=[]) # x, dia, pressure, volume velosity
    push!(prs,(x,men.df, v[1],v[2]))
    
    while men != nothing
        # @show men.tm
        if men.child != nothing
            push!(prs, (x,0.0, v[1],v[2]))
        end
        if men.r > 0
            v = inv(men.tm)*v
        end
        x += men.r
        push!(prs,(x,men.df,v[1],v[2]))
        men = actualnext(men)
    end
    return prs
end

function calcprs(mentable, frequency::Float64; slice=1.0, temperature=24.0, pressure = 2.0e-2)
    slice_mensur!(mentable["MAIN"], slice/1000.0) # unit conv here
    
    # calc impedance first
    prm = initcalcparam(temperature=temperature)
    wf = 2pi*frequency
    men = mentable["MAIN"]
    impedance!(wf,men;prm...)
    
    prs = calc_pressure(pressure,men) # return (x,p,U) dataframe
    prs[:x] .*= 1000.0 # x unit conv
    prs[:D] .*= 1000.0 # dia unit conv
    return prs
end

"""
    calcprs(filepath, frequency; kwargs...)
    
    Calculate pressure distribution for given frequency.
    Mensur will be sliced by `slice=` value (default 1mm).
    Pressure at end must be supplied (default = 2.0e-2 Pa = 60dBSPL).

    Returns dataframe of (x, D, p, U),
    where x:length from top, D,p,U are diameter,pressure,volume velosity at x.
"""
function calcprs(fpath::String,frequency::Float64; params...)
    mentable = readxmen(fpath)
    return calcprs(mentable,frequency; params...)
end