__precompile__(false)

"""
Mensur handling module.
"""
module Mensur

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

"""
create new mensur 
"""
function men_create(df::Number, db::Number, r::Number, cm::String, gr::String)
    men = Men()

    men.params = Dict(:df=>df, :db=>db, :r=>r, :comment=>cm, :group=>gr)
    men.vars = Dict(:zi=>zero(Complex), :zo=>zero(Complex),
        :vi=>zeros(Complex,2), 
        :vo=>zeros(Complex,2),
        :tm=>zeros(Complex,2,2)
    )
    men.childinfo = Dict(:type=>0, :name=>"", :ratio=>0.0)

    men.prev = nothing
    men.next = nothing
    men.parent = nothing
    men.child = nothing

    return(men)
end
men_create() = men_create(0,0,0,"","")

# function 
function men_append!(prev::Men, curr::Men)
    prev.next = curr
    curr.prev = prev
    return(curr) 
end

function men_setchildinfo()
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

function men_bykeyword(cur::Men, wd::Array{String})
    # not implemented
end

"""Build mensur from text lines""" 
function men_build(lns)
    mentable = Dict{String,Men}()
    grouplist = []
    groupnames = []
    cur = nothing

    for ln in lns
        ln = eat_comment(ln)
        ln = replace(ln,r"[\t\s]+"=>"") # remove all white spaces
        
        if occursin("=",ln)
            # evaluate var definition line
            eval(Meta.parse(ln))
        else
            # normal df,db,r,comment line
            wd = split(ln,",")

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
                cur = men_bykeyword(cur,wd)
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
                    push!(mentable,groupnames[end] => cur)
                else
                    cur = men_append!(cur,men)
                end
            else    
                throw(AssertionError("string $ln is invalid."))
            end

        end
    end
    men = men_top(cur)
    return(men) 
end

function men_readfile(path)
    lns = readlines(path)

    men = nothing
    if length(lns) > 0 
        men = men_build(lns)    
    end
    return(men)
end

# export structs, functions
export Men
export men_create, men_append!, men_top, men_end

end # module