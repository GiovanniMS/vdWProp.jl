using Polynomials

using EngThermBase

# par pressao e temperatura

function p3vr(Pr,Tr) #Polinomio de 3 grau de vr utilizando Pr e Tr
    
    return roots(Poly([(-3/8),(9/8),(-Tr - (Pr/8)),(3*Pr/8)]))
    
end

# par pressao e energia interna

function p3vr(Pr,ur,ϕ,C1) #Polinomio de 3 grau de vr utilizando Pr e ur
    
    return roots(Poly([-3,(9 - (9/ϕ)),(-Pr - (ur - C1)*(3/ϕ)),3*Pr]))
    
end

# par pressao e entalpia

function p3vr2(Pr,hr,ϕ,C1) #Polinomio de 3 grau de vr utilizando Pr e hr
    
    return roots(Poly([(-ϕ),(3*ϕ - 3),(-(ϕ*Pr/3) - (hr - C1)),(ϕ*Pr + Pr)]))
    
end

# par energia interna e entalpia

function p2vr(ur,hr,ϕ,C1) #Polinomio de 3 grau de vr utilizando Pr e hr
    
    return roots(Poly([ϕ,((-ϕ/3)*(ur - C1) + 3*ϕ + 3 + (ϕ/3)*(hr - C1) - 6*ϕ),((ϕ + 1)*(ur - C1) - ϕ*(hr - C1))]))
    
end

# das 7 equacoes mais problematicas, as 4 acima puderam ser tranformadas em polinomios, ja as 3 restantes
# precisarao ser feitas numericamente

using Roots

# par Pressao e entropia 

function ParPs(Pr,sr,ϕ,C2) 
    
    f(vr) = ((8*exp((3/(8*ϕ))*(sr + C2)))/((3*vr - 1)^(1 + (1/ϕ)))) - (3/vr^2) - Pr
    
    return find_zero(f,1,Order1()) #metodo da secante

end

# par energia interna e entropia

function Parus(ur,sr,ϕ,C1,C2) 
    
    f(vr) = (8*ϕ/3)*(exp((3/(8*ϕ))*(sr + C2)))*((3*vr - 1)^(-1/ϕ)) - ur - (3/vr) + C1
    
    return find_zero(f,0.5,Order1()) #metodo da secante

end

# par entalpia e entropia

function Parhs(hr,sr,ϕ,C1,C2) 
    
    f(vr) = ((8*ϕ/3) + (8*vr/(3*vr - 1)))*exp((3/(8*ϕ))*(sr + C2))*((3*vr - 1)^(-1/ϕ)) - hr + C1 - (6/vr)
    
    return find_zero(f,0.5,Order1()) #metodo da secante

end

# resolvendo agora os 3 anteriores por Newton-Raphson

using ForwardDiff

# par Pressao e entropia 

function ParPsN(Pr,sr,ϕ,C2) 
    
    D(f) = vr -> ForwardDiff.derivative(f,float(vr))
    
    x0 = 0.5
    
    f(vr) = ((8*exp((3/(8*ϕ))*(sr + C2)))/((3*vr - 1)^(1 + (1/ϕ)))) - (3/vr^2) - Pr
    
    return find_zero((f,D(f)), x0, Roots.Newton()) #metodo de Newton-Raphson

end

# par energia interna e entropia

function ParusN(ur,sr,ϕ,C1,C2) 
    
    D(f) = vr -> ForwardDiff.derivative(f,float(vr))
    
    x0 = 0.5
    
    f(vr) = (8*ϕ/3)*(exp((3/(8*ϕ))*(sr + C2)))*((3*vr - 1)^(-1/ϕ)) - ur - (3/vr) + C1
    
    return find_zero((f,D(f)), x0, Roots.Newton()) #metodo de Newton-Raphson

end

# par entalpia e entropia

function ParhsN(hr,sr,ϕ,C1,C2) 
    
    D(f) = vr -> ForwardDiff.derivative(f,float(vr))
    
    x0 = 0.5
    
    f(vr) = ((8*ϕ/3) + (8*vr/(3*vr - 1)))*exp((3/(8*ϕ))*(sr + C2))*((3*vr - 1)^(-1/ϕ)) - hr + C1 - (6/vr)
    
    return find_zero((f,D(f)), x0, Roots.Newton()) #metodo de Newton-Raphson
    
end