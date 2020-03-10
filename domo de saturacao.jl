### Domo de Saturacao ###

# para obtencao do domo sao necessarios ℘ V1 V2 vr1 vr2 Tr_sat Pr_sat

###################################################################################

#FUNCOES

#V1(vr1) = 3vr1 - 1

#V2(vr2) = 3vr2 - 1

#℘(V1,V2) = V1 / V2

#V1(℘) = ((℘^2)*(-2℘ + (℘ + 1)*log(℘) + 2)) / ((℘^2) - 2*℘*log(℘) -1)

#V2(℘) = (℘*(-2℘ + (℘ + 1)*log(℘) + 2)) / ((℘^2) - 2*℘*log(℘) -1)

V2(℘) = (2*log(℘)/(℘-1) - 1/℘ - 1) / (2 - (℘+1)*log(℘)/(℘-1))

V1(℘) = V2(℘)*℘

vr1(V1) = (V1 + 1) / 3

vr2(V2) = (V2 + 1) / 3

Tr_sat(V1,V2) = (27*V1*V2*(V1 + V2 + 2)) / (8*((V1 + 1)^2)*((V2 + 1)^2))

Pr_sat(V1,V2) = (27*(-V1 - V2 + (V1 + 1)*(V2 + 1) - 2)) / (((V1 + 1)^2)*((V2 + 1)^2))

#####################################################################################

#ARRAYS

℘list = range(0.0001, stop = 0.95, length = 10000)

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

#######################################################################################

#PLOT Pxv

using Makie

x1 = Array{Float64,1}(vr1list)

x2 = Array{Float64,1}(vr2list)

y = Array{Float64,1}(Pr_sat_list)

scene = lines(x1, y, color = :blue)

lines!(scene, x2, y, color = :green)

display(scene)

##########################################################################################

#PLOT Pxv(log)

using Makie

x1 = Array{Float64,1}(vr1listlog)

x2 = Array{Float64,1}(vr2listlog)

y = Array{Float64,1}(Pr_sat_list)

scene = lines(x1, y, color = :blue)

lines!(scene, x2, y, color = :green)

display(scene)