__precompile__()

"""
Mensur handling module.
"""
module Mensur
using DataFrames
using SpecialFunctions
using LinearAlgebra
using Struve

# constants
OPEN = 1
CLOSE = 0
HALF = 0.5
HEAD = 0
LAST = 1

type_keywords = Dict("SPLIT"=>:split, 
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

group_keywords = Dict("MAIN"=>:main,
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

"""
Men struct which contains basic unit for expressing mensur structure.
"""
mutable struct Men
    # data to store
    params::Dict{Symbol,Any}
    childinfo::Dict{Symbol,Any}    
    vars::Dict{Symbol,Any}
    
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
function men_create(df::Number, db::Number, r::Number, cm::AbstractString, gr::String)
    men = Men()

    men.params = Dict(:df=>df, :db=>db, :r=>r, :comment=>cm, :group=>gr)
    men.vars = Dict(:zi=>zero(Complex), :zo=>zero(Complex),
        :vi=>zeros(Complex,2), 
        :vo=>zeros(Complex,2),
        :tm=>zeros(Complex,2,2)
    )
    men.childinfo = Dict(:type=>"", :name=>"",:ratio=>0.0)

    men.prev = nothing
    men.next = nothing
    men.parent = nothing
    men.child = nothing

    return(men)
end
men_create() = men_create(0,0,0,"","")

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
    men.childinfo[:type]=t
    men.childinfo[:name]=n
    men.childinfo[:ratio]=r
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
    return(cur.params[:df], cur.params[:db],cur.params[:r],cur.params[:comment])
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
    gr = cur.params[:group]

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
        cnm = cur.childinfo[:name]
        typ = cur.childinfo[:type]
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

function men_readfile(path)
    lns = readlines(path)

    mentable = []
    if length(lns) > 0 
        mentable = men_build(lns)    
    end
    return(mentable)
end

function men_print(men::Men)
    gr = men.params[:group]
    if gr != "MAIN"
        println("GROUP,",gr)
    end

    while men != nothing
        df,db,r,c = getfbr(men)
        if men.child != nothing
            println(men.childinfo[:type],",",men.childinfo[:name],",",men.childinfo[:ratio])
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
    if isbranch(men.childinfo[:type])
        e = men_end(men.child)
        if e.parent != nothing & ismerge(e.parent.childinfo[:type])
            return(e.parent)
        else
            return(nothing)
        end
    else
        return(nothing)
    end
end

function update_calcparams!(params)
    GMM = 1.4  # specific head ratio
    PR = 0.72  # Prandtl number

    tp = params["temperature"]

    params["c0"] = c0 = 331.45 * sqrt(tp / 273.16 + 1)
    params["rho"] = ρ = 1.2929 * (273.16 / (273.16 + tp))
    params["rhoc0"] = ρc0 = ρ * c0
    params["mu"] = μ = (18.2 + 0.0456*(tp - 25)) * 1.0e-6  # viscosity constant. Linear approximation from Scientific Dictionary.
    params["nu"] = ν= μ/ρ  # dynamic viscous constant.
    params["dmp"] = (1+(GMM-1)/sqrt(PR))
end

function radimp(wf::Float64, dia::Float64, params)
    if wf > 0
        if dia == 0.0
            zr = complex(Inf) # closed end
        else
            if params["radiation"] == "NONE"
                zr = complex(0)  # simple open end impedance
            else
                s = dia^2*pi/4.0
                k = wf/params["c0"]
                x = k*dia
                rhoc0 = params["rhoc0"]

                re = rhoc0/s*(1 - besselj1(x)/x*2)  # 1st order bessel function.
                img = rhoc0/s*Struve.H(1,x)/x*2  # 1st order struve function.

                if params["radiation"] == "BAFFLE"
                    zr = complex(re,imag)
                elseif params["radiation"] == "PIPE"
                    # real is about 0.5 times and imaginary is 0.7 times when without frange.
                    zr = complex(0.5*re,0.7*img)
                end
            end
        end
        return(zr) 
    end
end

function calc_transmission(wf::Float64,men::Men, params::Dict{String,Any})
    df = men.params[:df]
    db = men.params[:db]
    r = men.params[:r]
    dmp = params["dmp"]
    nu = params["nu"]
    c0 = params["c0"]
    rhoc0 = params["rhoc0"]

    d = (df + db)*0.5
    aa = dmp * sqrt(2*wf*nu)/c0/d  # wall dumping factor
    k = sqrt((wf/c0)*(wf/c0 - 2*complex(-1,1)*aa))  # complex wave number including wall dumping
    x = k * r
    cc = cos(x)
    ss = sin(x)

    tm = zeros(Complex,2,2)

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

function zo2zi(tm::Array{Complex,2},zo::Complex)
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
    while men != nothing & men != men1
        m *= men.vars[:tm]
        men = men.prev
    end
    m *= men1.vars[:tm]
    
    return(m)
end

function child_impedance(wf::Float64,men::Men, params::Dict{String,Any})
    if issplit(men.childinfo[:type])
        # split (tonehole) type.
        input_impedance(wf, men.child)  # recursive call for input impedance
        if men.childinfo[:ratio] == 0
            men.vars[:zo] = men.next.vars[:zi]
        else
            z1 = men.child.vars[:zi] / men.childinfo[:ratio]  # adjust blending ratio
            z2 = men.next.vars[:zi]
            if z1 == zero(Complex) & z2 == zero(Complex)
                z = zero(Complex)
            else
                z = z1*z2/(z1+z2)
            end
            men.vars[:zo] = z
        end
    elseif isbranch(men.childinfo[:type]) & men.childinfo[:ratio] > 0
        # multiple tube connection
        input_impedance(wf, men.child)
        m = transmission_matrix(men.child, nothing)
        jnt = joint_mensur(men)
        n = transmission_matrix(men.next, jnt)

        # section area adjustment
        if men.childinfo[:ratio] == 1
            m[1,2] = complex(Inf)
        else
            m[1,2] /= (1 - men.childinfo[:ratio])
        end
        m[2,1] *= (1 - men.childinfo[:ratio])

        if men.childinfo[:ratio] == 0
            n[1,2] = complex(Inf)
        else
            n[1,2] /= men.childinfo[:ratio]
        end
        n[2,1] *= men.childinfo[:ratio]

        z2 = jnt.next.vars[:zi]
        dv = (m[2,2]*n[1,2] + m[1,2]*n[2,2] + (
            (m[1,2] + n[1,2])*(m[2,1] + n[2,1]) - (m[1,1] - n[1,1])*(m[2,2] - n[2,2]))*z2)
        if dv != zero(Complex)
            z = (m[1,2]*n[1,2] + (m[1,2]*n[1,1] + m[1,1]*n[1,2])*z2)/dv
        else
            z = 0
        end
        men.vars[:zo] = z
    elseif isaddon(men.childinfo[:type]) & men.childinfo[:ratio] > 0
        # this routine will not called until 'ADDON(LOOP)' type of connection is implemented.
        input_impedance(wf, men.child)
        m = transmission_matrix(men.child, nothing)
        z1 = m[1,2]/(m[1,2]*m[2,1]-(1-m[1,1])*(1-m[2,2]))
        z2 = men.next.vars[:zi]
        if men.childinfo[:ratio] == 0
            men.vars[:zo] = men.next.vars[:zi]
        elseif men.childinfo[:ratio] == 1
            men.vars[:zo] = z1
        else
            z1 /= men.childinfo[:ratio]
            z2 /= (1 - men.childinfo[:ratio])
            if z1 == zero(Complex) & z2 == zero(Complex)
                z = zero(Complex)
            else
                z = z1*z2/(z1+z2)
            end
            men.vars[:zo] = z
        end
    end
end

"""recursively called function to calculate impedance at each mensur node"""
function calc_impedance!(wf::Float64,men::Men, params::Dict{String,Any})
    if men.child != nothing
        child_impedance(wf, men, params)
    elseif men.next != nothing
        men.vars[:zo] = men.next.vars[:zi]
    end

    if men.params[:r] > 0.0
        men.vars[:tm] = calc_transmission(wf,men,params)
        men.vars[:zi] = zo2zi(men.vars[:tm], men.vars[:zo])
    else
        # length 0
        men.vars[:zi] = men.vars[:zo]
    end
end

"""calculate impedance at given frequency """
function impedance!(wf::Float64, men::Men, params::Dict{String,Any})
    if wf > 0
        cur = men_end(men)
        # end impedance
        cur.vars[:zo] = radimp(wf, cur.params[:df],params)

        while cur != men
            calc_impedance!(wf, cur,params)
            cur = cur.prev
        end
        calc_impedance!(wf, men,params)
        return(men.vars[:zi])
    else
        return(zero(Complex))
    end
end

function input_impedance(mentable,params)
    minf = params["minfreq"]
    maxf = params["maxfreq"]
    sf = params["stepfreq"]
    tp = params["temperature"]
    radtype = params["radiation"]

    update_calcparams!(params)

    imped = DataFrame()
    frq = range(minf, stop=maxf, step=sf)
    imped[:frq] = frq
    imped[:imp] = zeros(Complex,length(frq))

    men = mentable["MAIN"]

    ss = men.params[:df]^2 * pi/4

    for i in 1:length(frq)
        wf = frq[i]*2*pi
        imped[:imp][i] = ss * impedance!(wf,men,params)
    end

    return(imped)
end

# export structs, functions
export Men
export men_create, men_append!, men_top, men_end, men_readfile, men_printtable, input_impedance

end # module