static notenames=Dict(0=>"A",1=>"Bb",2=>"B",3=>"C",4=>"Db",5=>"D",6=>"Eb",
7=>"E",8=>"F",9=>"Gb",10=>"G",11=>"Ab")

static convlist = Dict("C#"=>"Db","D#"=>"Eb","F#"=>"Gb","G#"=>"Ab","A#"=>"Bb")

"""Calculate dBSPL of pressure x"""
dBSPL(x) = 20log10(abs(x/2.0e-5))

"""
Hz to Note name
    hztonote(f,pitch=440.0)
"""
function hztonote(f,pitch=440.0)
    c = 1200log2(f/pitch)
    n = Int(mod(c÷100,12)) # note number
    cc = c%100 # cents remainder
    if cc > 50 
        cc = 100-cc
        n+=1
    elseif cc <= -50
        cc = 100+cc
        n-=1
        if n < 0
            n += 12
        end
    end
    octbase = pitch*2^(-95/120) # c'-50cent
    oc =Int(floor(log2(f/octbase)))
    nm = notenames[n]
    if oc > -2
        nm = lowercase(nm)
    end
    if oc > -1
        om = "'"^(oc+1)
    elseif oc < -2
        om = ","^(-2-oc)
    else
        om = ""
    end
    nm *= om
    return nm,cc
end

function notenameconv(name)
    if name in keys(convlist)
        nm = convlist[name]
    else
        nm = name
    end
    nm
end

function notenumber(s)
    ss = notenameconv(uppercasefirst(s))
    for p in notenames
        if p.second == ss
            return p.first
        end
    end
    return -1
end

"""
Convert note and cent to Hz
"""
function notetohz(name, cent, pitch = 440.0)
    if isuppercase(name[1])
        om = replace(name,r"[^,]"=>"")
        nm = replace(name,","=>"")
        oc = -3 - length(om)
    else
        om = replace(name,r"[^']"=>"")
        nm = replace(name,"'"=>"")
        oc = length(om)-2
    end
    nn = notenumber(nm)
    if nn < 3 # lower than c
        oc+=1
    end
    f = pitch * 2^((nn*100+cent)/1200+oc)
#     @show nm,nn,oc
    return f
end

