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

function men_bykeyword!(cur::Men, wd::Array{String,1})
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
        men_setchildinfo!(cur,w,name,ratio)
    end
    cur = men_append!(cur,men)
    return(cur)
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
                    cur = men_bykeyword!(cur,wd)
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
        println(df*1000.0,",",db*1000.0,",",r*1000.0,",",c)
        if men.child != nothing
            println(men.childinfo[:type],",",men.childinfo[:name],",",men.childinfo[:ratio])
            men = men.next # skip one !
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

# export structs, functions
export Men
export men_create, men_append!, men_top, men_end, men_readfile, men_printtable

end # module