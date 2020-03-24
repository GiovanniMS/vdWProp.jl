using Polynomials

using EngThermBase

# Pr - Tr

function p3vr(Pr,Tr) # 3 degree polynomial of vr using Pr and Tr
    
    return roots(Poly([(-3/8),(9/8),(-Tr - (Pr/8)),(3*Pr/8)]))
    
end

# Pr - ur

function p3vr(Pr,ur,ϕ,C1) # 3 degree polynomial of vr using Pr and ur
    
    return roots(Poly([-3,(9 - (9/ϕ)),(-Pr - (ur - C1)*(3/ϕ)),3*Pr]))
    
end

# Pr - hr

function p3vr2(Pr,hr,ϕ,C1) # 3 degree polynomial of vr using Pr and hr
    
    return roots(Poly([(-ϕ),(3*ϕ - 3),(-(ϕ*Pr/3) - (hr - C1)),(ϕ*Pr + Pr)]))
    
end

# ur - hr

function p2vr(ur,hr,ϕ,C1) # 3 degree polynomial of vr using Pr and hr
    
    return roots(Poly([ϕ,((-ϕ/3)*(ur - C1) + 3*ϕ + 3 + (ϕ/3)*(hr - C1) - 6*ϕ),((ϕ + 1)*(ur - C1) - ϕ*(hr - C1))]))
    
end

# from the 7 problematic equations, the 4 above could be reduced to polynomials
# in the other 3, numeric solution is needed

using Roots

# Pr - sr

function ParPs(Pr,sr,ϕ,C2) 
    
    f(vr) = ((8*exp((3/(8*ϕ))*(sr + C2)))/((3*vr - 1)^(1 + (1/ϕ)))) - (3/vr^2) - Pr
    
    return find_zero(f,1,Order1()) # secant method

end

# ur - sr

function Parus(ur,sr,ϕ,C1,C2) 
    
    f(vr) = (8*ϕ/3)*(exp((3/(8*ϕ))*(sr + C2)))*((3*vr - 1)^(-1/ϕ)) - ur - (3/vr) + C1
    
    return find_zero(f,0.5,Order1()) # secant method

end

# hr - sr

function Parhs(hr,sr,ϕ,C1,C2) 
    
    f(vr) = ((8*ϕ/3) + (8*vr/(3*vr - 1)))*exp((3/(8*ϕ))*(sr + C2))*((3*vr - 1)^(-1/ϕ)) - hr + C1 - (6/vr)
    
    return find_zero(f,0.5,Order1()) # secant method

end

# same functions using Newton-Raphson method

using ForwardDiff

# Pr - sr

function ParPsN(Pr,sr,ϕ,C2) 
    
    D(f) = vr -> ForwardDiff.derivative(f,float(vr))
    
    x0 = 0.5
    
    f(vr) = ((8*exp((3/(8*ϕ))*(sr + C2)))/((3*vr - 1)^(1 + (1/ϕ)))) - (3/vr^2) - Pr
    
    return find_zero((f,D(f)), x0, Roots.Newton()) # Newton-Raphson method
end

# ur - sr

function ParusN(ur,sr,ϕ,C1,C2) 
    
    D(f) = vr -> ForwardDiff.derivative(f,float(vr))
    
    x0 = 0.5
    
    f(vr) = (8*ϕ/3)*(exp((3/(8*ϕ))*(sr + C2)))*((3*vr - 1)^(-1/ϕ)) - ur - (3/vr) + C1
    
    return find_zero((f,D(f)), x0, Roots.Newton()) # Newton-Raphson method

end

# hr - sr

function ParhsN(hr,sr,ϕ,C1,C2) 
    
    D(f) = vr -> ForwardDiff.derivative(f,float(vr))
    
    x0 = 0.5
    
    f(vr) = ((8*ϕ/3) + (8*vr/(3*vr - 1)))*exp((3/(8*ϕ))*(sr + C2))*((3*vr - 1)^(-1/ϕ)) - hr + C1 - (6/vr)
    
    return find_zero((f,D(f)), x0, Roots.Newton()) # Newton-Raphson method
    
end