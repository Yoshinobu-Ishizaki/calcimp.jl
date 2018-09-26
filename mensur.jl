__precompile__(false)

"""
Mensur handling module.
"""
module Mensur
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

function men_top(curr::Men)
    c = curr
    if c != nothing
        while c.prev != nothing
            c = c.prev
        end
    end
    return(c)
end

# export structs, functions
export Men
export men_create, men_append!, men_top

end # module