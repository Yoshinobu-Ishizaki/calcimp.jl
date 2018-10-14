"""Calculate dBSPL of pressure x"""
dBSPL(x) = 20log10(abs(x/2.0e-5))

static notenames=Dict(0=>"A",1=>"Bb",2=>"B",3=>"C",4=>"Db",5=>"D",6=>"Eb",
    7=>"E",8=>"F",9=>"Gb",10=>"G",11=>"Ab")

"""
Hz to Note name
    hztonote(f,pitch=440.0)
"""
function hztonote(f,pitch=440.0)
    # todo : add octave mark
    c = 1200log2(f/pitch)
    n = Int(mod(c÷100,12)) # note number
    cc = c%100 # cents remainder
    if cc > 50 
        cc = 100-cc
        n+=1
    elseif cc < -50
        cc = 100+cc
        n-=1
        if n < 0
            n += 12
        end
    end
    return notenames[n],cc
end

function notenumber(s)
    for p in notenames
        if p.second == s
            return p.first
        end
    end
    return -1
end

function notetohz(name, cent, pitch = 440.0)
end