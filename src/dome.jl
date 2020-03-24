### Saturation Dome ###

#  ℘ V1 V2 vr1 vr2 Tr_sat Pr_sat - are needed

###################################################################################

#Functions

V2(℘) = (2*log(℘)/(℘-1) - 1/℘ - 1) / (2 - (℘+1)*log(℘)/(℘-1))

V1(℘) = V2(℘)*℘

vr1(V1) = (V1 + 1) / 3

vr2(V2) = (V2 + 1) / 3

Tr_sat(V1,V2) = (27*V1*V2*(V1 + V2 + 2)) / (8*((V1 + 1)^2)*((V2 + 1)^2))

Pr_sat(V1,V2) = (27*(-V1 - V2 + (V1 + 1)*(V2 + 1) - 2)) / (((V1 + 1)^2)*((V2 + 1)^2))

#####################################################################################

# Arrays

points = 10000

℘list = range(0.0001, stop = 0.95, length = points)

V2list = []

V1list = []

vr1list = []

vr2list = []

Pr_sat_list = []

Tr_sat_list = [] 

vr1listlog = []

vr2listlog = []

for i in ℘list
    
    append!(V2list,Float64(V2(i)))
    
    append!(V1list,Float64(V1(i)))
    
    append!(vr1list,Float64(vr1(V1(i))))
    
    append!(vr2list,Float64(vr2(V2(i))))
    
    append!(Pr_sat_list,Float64(Pr_sat(V1(i),V2(i))))
    
    append!(Tr_sat_list,Float64(Tr_sat(V1(i),V2(i))))
    
    append!(vr1listlog,Float64(log(vr1(V1(i)))))
    
    append!(vr2listlog,Float64(log(vr2(V2(i)))))
    
end

########################################################################################

#PLOT Prxvr e Prxlogvr

# "using Makie" is needed only to get the plots, but can be used externally, thats why it is commented

#using Makie

function PlotDomoPrVr(name::String)
    
    if name == "Pr x vr"

        x1 = Array{Float64,1}(vr1list)

        x2 = Array{Float64,1}(vr2list)

        y = Array{Float64,1}(Pr_sat_list)

        scene = lines(x1, y, color = :blue)

        lines!(scene, x2, y, color = :green)

        Makie.save("Pr x vr.jpg", scene)
        
    elseif name == "Pr x log(vr)"

        x1 = Array{Float64,1}(vr1listlog)

        x2 = Array{Float64,1}(vr2listlog)

        y = Array{Float64,1}(Pr_sat_list)

        scene = lines(x1, y, color = :blue)

        lines!(scene, x2, y, color = :green)

        Makie.save("Pr x log(vr).jpg", scene)
        
    else
        
        print("name must be Pr x vr or Pr x log(vr)")
        
    end
    
end

########################################################################################

# using the dome properties, it is possible to get the other associated properties
# the properties bellow depends on the constants \phi and integration constants
# it's needed to plot all the domes for each \phi

Zc = 3/8

ar(vr,Tr,ϕ,C1) = C1 + Tr*(C2 - (ϕ/Zc)*log(Tr) + ϕ/Zc) - (8/3)*Tr*log(3vr - 1) - (3/vr)

sr(vr,Tr,ϕ,C2) = (8/3)*((log(3vr - 1)) + ϕ*log(Tr)) - C2

ur(vr,Tr,ϕ,C1) = C1 + (8/3)*Tr*ϕ - (3/vr)

hr(vr,Tr,ϕ,C1) = C1 + (8/3)*Tr*ϕ + ((8*Tr*vr)/(3vr - 1)) - (6/vr)

# C1 = C2 = 0
# case 1 - ϕ = 3/2
# case 2 - ϕ = 5/2
# case 3 - ϕ = 7/2

ϕ1 = 3/2

ϕ2 = 5/2

ϕ3= 7/2

C=0

# Functions to get the plots and properties for ur, hr, and sr

function domoprop(ϕ::Number,C1::Number,C2::Number,plot::String, array::String)
    
    if ϕ == ϕ1 || ϕ == ϕ2 || ϕ == ϕ3

        sr1list = []

        ur1list = []

        hr1list = []

        sr2list = []

        ur2list = []

        hr2list = []

        for i in 1:points

            append!(sr1list,Float64(sr(vr1list[i],Tr_sat_list[i],ϕ,C2)))

            append!(ur1list,Float64(ur(vr1list[i],Tr_sat_list[i],ϕ,C1)))

            append!(hr1list,Float64(hr(vr1list[i],Tr_sat_list[i],ϕ,C1)))

        end

        for i in 1:points

            append!(sr2list,Float64(sr(vr2list[i],Tr_sat_list[i],ϕ,C2)))

            append!(ur2list,Float64(ur(vr2list[i],Tr_sat_list[i],ϕ,C1)))

            append!(hr2list,Float64(hr(vr2list[i],Tr_sat_list[i],ϕ,C1)))

        end
        
        if plot == "yes"

            x1 = Array{Float64,1}(sr1list)

            x2 = Array{Float64,1}(sr2list)

            y = Array{Float64,1}(Pr_sat_list)

            scene = lines(x1, y, color = :blue)

            lines!(scene, x2, y, color = :green)

            Makie.save("Pr x sr.jpg", scene)

            x1 = Array{Float64,1}(ur1list)

            x2 = Array{Float64,1}(ur2list)

            y = Array{Float64,1}(Pr_sat_list)

            scene = lines(x1, y, color = :blue)

            lines!(scene, x2, y, color = :green)

            Makie.save("Pr x ur.jpg", scene)

            x1 = Array{Float64,1}(hr1list)

            x2 = Array{Float64,1}(hr2list)

            y = Array{Float64,1}(Pr_sat_list)

            scene = lines(x1, y, color = :blue)

            lines!(scene, x2, y, color = :green)

            Makie.save("Pr x hr.jpg", scene)
            
        elseif plot == "no" 
            
            if array == "sr1"
            
                return sr1list
            
            elseif array == "ur1"
                
                return ur1list
                
            elseif array == "hr1"
                
                return hr1list
                
            elseif array == "sr2"
                
                return sr2list
                
            elseif array == "ur2"
                
                return ur2list
                
            elseif array == "hr2"
                
                return hr2list
                
            else
                
                print("ERROR")
                
            end
            
        else
            
            print("plot must be yes or no")
            
        end
        
    else 
        
        print("ERROR, ϕ needs to be 3/2, 5/2 or 7/2")
        
    end
        
end

# Function to find the closest number in an array comparing to a specified number

function findclosest(array::Array,x::Number,p::Number)

    for i in 1:points
    
        y = x - array[i]
    
        if y < p
            
            return i
        
            break
            
        end    
        
        if y < 0 
            
            return i - 1
            
            break
            
        end
        
    end

end

# Dome properties for case 1

srlist1c1 = domoprop(ϕ1,C,C,"no","sr1")

urlist1c1 = domoprop(ϕ1,C,C,"no","ur1")

hrlist1c1 = domoprop(ϕ1,C,C,"no","hr1")

srlist2c1 = domoprop(ϕ1,C,C,"no","sr2")

urlist2c1 = domoprop(ϕ1,C,C,"no","ur2")

hrlist2c1 = domoprop(ϕ1,C,C,"no","hr2")

# Dome properties for case 2

srlist1c2 = domoprop(ϕ2,C,C,"no","sr1")

urlist1c2 = domoprop(ϕ2,C,C,"no","ur1")

hrlist1c2 = domoprop(ϕ2,C,C,"no","hr1")

srlist2c2 = domoprop(ϕ2,C,C,"no","sr2")

urlist2c2 = domoprop(ϕ2,C,C,"no","ur2")

hrlist2c2 = domoprop(ϕ2,C,C,"no","hr2")

# Dome properties for case 3

srlist1c3 = domoprop(ϕ3,C,C,"no","sr1")

urlist1c3 = domoprop(ϕ3,C,C,"no","ur1")

hrlist1c3 = domoprop(ϕ3,C,C,"no","hr1")

srlist2c3 = domoprop(ϕ3,C,C,"no","sr2")

urlist2c3 = domoprop(ϕ3,C,C,"no","ur2")

hrlist2c3 = domoprop(ϕ3,C,C,"no","hr2")