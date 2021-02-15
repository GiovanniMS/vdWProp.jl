### Saturation Dome ###

#  ℘ V1 V2 vr1 vr2 Tr_sat Pr_sat - are needed

#Functions

V2(℘) = (2*log(℘)/(℘-1) - 1/℘ - 1) / (2 - (℘+1)*log(℘)/(℘-1))

V1(℘) = V2(℘)*℘

vr1(V1) = (V1 + 1) / 3

vr2(V2) = (V2 + 1) / 3

Tr_sat(V1,V2) = (27*V1*V2*(V1 + V2 + 2)) / (8*((V1 + 1)^2)*((V2 + 1)^2))

Pr_sat(V1,V2) = (27*(-V1 - V2 + (V1 + 1)*(V2 + 1) - 2)) / (((V1 + 1)^2)*((V2 + 1)^2))

# Arrays

#points = 2500

#℘list = range(0.0008, stop = 0.0003, length = points)

#℘list0 = range(0.1, stop = 0.95, length = points)

#℘list1 = range(0.0001, stop = 0.0999, length = points)

#℘list2 = range(0.00000001, stop = 0.0000999, length = points)

#℘list3 = range(0.00000000001, stop = 0.0000000999, length = points)

#℘list4 = range(0.00000000000001, stop = 0.0000000000999, length = points)

#℘list5 = range(0.00000000000000001, stop = 0.0000000000000999, length = points)

#℘list6 = vcat(℘list1, ℘list0)

#℘list = vcat(℘list2, ℘list6)

#℘list8 = vcat(℘list3, ℘list7)

#℘list9 = vcat(℘list4, ℘list8)

#℘list = vcat(℘list5, ℘list9)

#points = 4999 # rev 4

points = 7109

℘list = [0.95]

for n in 1:points
    
    #append!(℘list, ℘list[n]/(1.0000054^n))
    
    #append!(℘list, ℘list[n]/(1.000001793^n))
    
    #append!(℘list, ℘list[n]/(1.00000011202^n))
    
    #append!(℘list, ℘list[n]/(1.0000010694^n))
    
    append!(℘list, ℘list[n]/(1.0000008872^n))
    
    #append!(℘list, ℘list[n]/(1.00000022175^n))
    
    n = n + 1
    
end

℘list2 = [0.99999]

for n in 1:100
    
    append!(℘list2, ℘list2[n]/(1.0005))
    
    n = n + 1
    
end

℘list = vcat(℘list2, ℘list)

℘list2 = []

points = length(℘list)

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

#PLOT Prxvr e Prxlogvr

# needs using Makie

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

# using the dome properties, it is possible to get the other associated properties
# the properties bellow depends on the constants \phi and integration constants
# it's needed to plot all the domes for each \phi

Zc = 3/8

ar(vr,Tr,ϕ,C1,C2) = C1 + Tr*(C2 - (ϕ/Zc)*log(Tr) + ϕ/Zc) - (8/3)*Tr*log(3vr - 1) - (3/vr)

sr(vr,Tr,ϕ,C2) = (8/3)*((log(3vr - 1)) + ϕ*log(Tr)) - C2

ur(vr,Tr,ϕ,C1) = C1 + (8/3)*Tr*ϕ - (3/vr)

hr(vr,Tr,ϕ,C1) = C1 + (8/3)*Tr*ϕ + ((8*Tr*vr)/(3vr - 1)) - (6/vr)

cpr(vr,Tr,ϕ) = (8*(4*Tr*(vr^3) + ϕ*(4*Tr*(vr^3) - (3*vr - 1)^2)))/(3*(4*Tr*(vr^3) - (3*vr - 1)^2))

# C1 = C2 = 0
# case 1 - ϕ = 3/2
# case 2 - ϕ = 5/2
# case 3 - ϕ = 7/2

ϕ1 = 3/2

ϕ2 = 5/2

ϕ3= 7/2

C=0

# Functions to get the plots and properties for ur, hr, and sr

# Plot needs using Makie

function domoprop(ϕ::Number,C1::Number,C2::Number,plot::String, array::String)
    
    if ϕ == ϕ1 || ϕ == ϕ2 || ϕ == ϕ3

        sr1list = []

        ur1list = []

        hr1list = []
        
        ar1list = []
        
        cpr1list = []

        sr2list = []

        ur2list = []

        hr2list = []
        
        ar2list = []
        
        cpr2list = []

        for i in 1:points

            append!(sr1list,Float64(sr(vr1list[i],Tr_sat_list[i],ϕ,C2)))

            append!(ur1list,Float64(ur(vr1list[i],Tr_sat_list[i],ϕ,C1)))

            append!(hr1list,Float64(hr(vr1list[i],Tr_sat_list[i],ϕ,C1)))
            
            append!(ar1list,Float64(ar(vr1list[i],Tr_sat_list[i],ϕ,C1,C2)))
            
            append!(cpr1list,Float64(cpr(vr1list[i],Tr_sat_list[i],ϕ)))

        end

        for i in 1:points

            append!(sr2list,Float64(sr(vr2list[i],Tr_sat_list[i],ϕ,C2)))

            append!(ur2list,Float64(ur(vr2list[i],Tr_sat_list[i],ϕ,C1)))

            append!(hr2list,Float64(hr(vr2list[i],Tr_sat_list[i],ϕ,C1)))
            
            append!(ar2list,Float64(ar(vr2list[i],Tr_sat_list[i],ϕ,C1,C2)))
            
            append!(cpr2list,Float64(cpr(vr2list[i],Tr_sat_list[i],ϕ)))

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
                
            elseif array == "ar1"
                
                return ar1list
                
            elseif array == "cpr1"
                
                return cpr1list
                
            elseif array == "sr2"
                
                return sr2list
                
            elseif array == "ur2"
                
                return ur2list
                
            elseif array == "hr2"
                
                return hr2list
                
            elseif array == "ar2"
                
                return ar2list
                
            elseif array == "cpr2"
                
                return cpr2list
                
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

# Dome properties for case 1

srlist1c1 = domoprop(ϕ1,C,C,"no","sr1")

urlist1c1 = domoprop(ϕ1,C,C,"no","ur1")

hrlist1c1 = domoprop(ϕ1,C,C,"no","hr1")

arlist1c1 = domoprop(ϕ1,C,C,"no","ar1")

cprlist1c1 = domoprop(ϕ1,C,C,"no","cpr1")

srlist2c1 = domoprop(ϕ1,C,C,"no","sr2")

urlist2c1 = domoprop(ϕ1,C,C,"no","ur2")

hrlist2c1 = domoprop(ϕ1,C,C,"no","hr2")

arlist2c1 = domoprop(ϕ1,C,C,"no","ar2")

cprlist2c1 = domoprop(ϕ1,C,C,"no","cpr2")

# Dome properties for case 2

srlist1c2 = domoprop(ϕ2,C,C,"no","sr1")

urlist1c2 = domoprop(ϕ2,C,C,"no","ur1")

hrlist1c2 = domoprop(ϕ2,C,C,"no","hr1")

arlist1c2 = domoprop(ϕ2,C,C,"no","ar1")

cprlist1c2 = domoprop(ϕ2,C,C,"no","cpr1")

srlist2c2 = domoprop(ϕ2,C,C,"no","sr2")

urlist2c2 = domoprop(ϕ2,C,C,"no","ur2")

hrlist2c2 = domoprop(ϕ2,C,C,"no","hr2")

arlist2c2 = domoprop(ϕ2,C,C,"no","ar2")

cprlist2c2 = domoprop(ϕ2,C,C,"no","cpr2")

# Dome properties for case 3

srlist1c3 = domoprop(ϕ3,C,C,"no","sr1")

urlist1c3 = domoprop(ϕ3,C,C,"no","ur1")

hrlist1c3 = domoprop(ϕ3,C,C,"no","hr1")

arlist1c3 = domoprop(ϕ3,C,C,"no","ar1")

cprlist1c3 = domoprop(ϕ3,C,C,"no","cpr1")

srlist2c3 = domoprop(ϕ3,C,C,"no","sr2")

urlist2c3 = domoprop(ϕ3,C,C,"no","ur2")

hrlist2c3 = domoprop(ϕ3,C,C,"no","hr2")

arlist2c3 = domoprop(ϕ3,C,C,"no","ar2")

cprlist2c3 = domoprop(ϕ3,C,C,"no","cpr2")

function Domelist(ϕ::Number, x::String)
    
    if x == "ur1"
        
        if ϕ == ϕ1
            
            return urlist1c1
            
        elseif ϕ == ϕ2
            
            return urlist1c2
            
        elseif ϕ == ϕ3
            
            return urlist1c3
            
        end
        
    elseif x == "ur2"
        
        if ϕ == ϕ1
            
            return urlist2c1
            
        elseif ϕ == ϕ2
            
            return urlist2c2
            
        elseif ϕ == ϕ3
            
            return urlist2c3
            
        end
        
    elseif x == "hr1"
        
        if ϕ == ϕ1
            
            return hrlist1c1
            
        elseif ϕ == ϕ2
            
            return hrlist1c2
            
        elseif ϕ == ϕ3
            
            return hrlist1c3
            
        end
        
    elseif x == "hr2"
        
        if ϕ == ϕ1
            
            return hrlist2c1
            
        elseif ϕ == ϕ2
            
            return hrlist2c2
            
        elseif ϕ == ϕ3
            
            return hrlist2c3
            
        end
        
    elseif x == "sr1"
        
        if ϕ == ϕ1
            
            return srlist1c1
            
        elseif ϕ == ϕ2
            
            return srlist1c2
            
        elseif ϕ == ϕ3
            
            return srlist1c3
            
        end
        
    elseif x == "sr2"
        
        if ϕ == ϕ1
            
            return srlist2c1
            
        elseif ϕ == ϕ2
            
            return srlist2c2
            
        elseif ϕ == ϕ3
            
            return srlist2c3
            
        end
        
    elseif x == "ar1"
        
        if ϕ == ϕ1
            
            return arlist1c1
            
        elseif ϕ == ϕ2
            
            return arlist1c2
            
        elseif ϕ == ϕ3
            
            return arlist1c3
            
        end
        
    elseif x == "ar2"
        
        if ϕ == ϕ1
            
            return arlist2c1
            
        elseif ϕ == ϕ2
            
            return arlist2c2
            
        elseif ϕ == ϕ3
            
            return arlist2c3
            
        end
        
    elseif x == "cpr1"
        
        if ϕ == ϕ1
            
            return cprlist1c1
            
        elseif ϕ == ϕ2
            
            return cprlist1c2
            
        elseif ϕ == ϕ3
            
            return cprlist1c3
            
        end
        
    elseif x == "cpr2"
        
        if ϕ == ϕ1
            
            return cprlist2c1
            
        elseif ϕ == ϕ2
            
            return cprlist2c2
            
        elseif ϕ == ϕ3
            
            return cprlist2c3
            
        end
       
    end
    
end