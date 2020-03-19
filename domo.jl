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

########################################################################################

# utilizando as propriedades do domo podemos encontrar as demais propriedades associadas
# a partir das demais propriedades, havera variacao dependendo do tipo de gas (por causa do ϕ) e dependendo das constantes escolhidas
# portanto deve-se escolher as constantes e fazer um domo para cada tipo de gas (3 domos)

Zc = 3/8

ar(vr,Tr,ϕ,C1) = C1 + Tr*(C2 - (ϕ/Zc)*log(Tr) + ϕ/Zc) - (8/3)*Tr*log(3vr - 1) - (3/vr)

sr(vr,Tr,ϕ,C2) = (8/3)*((log(3vr - 1)) + ϕ*log(Tr)) - C2

ur(vr,Tr,ϕ,C1) = C1 + (8/3)*Tr*ϕ - (3/vr)

hr(vr,Tr,ϕ,C1) = C1 + (8/3)*Tr*ϕ + ((8*Tr*vr)/(3vr - 1)) - (6/vr)

# inicialmente serao considerados C1 = C2 = 0
# caso 1 - ϕ = 3/2
# caso 2 - ϕ = 5/2
# caso 3 - ϕ = 7/2

ϕ1 = 3/2

ϕ2 = 5/2

ϕ3= 7/2

C=0

# Funcao para obter as propriedades e plots em todos os casos

function domoprop(ϕ::Number,C1::Number,C2::Number)
    
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
        
    else 
        
        print("ERROR, ϕ needs to be 3/2, 5/2 or 7/2")
        
    end
        
end