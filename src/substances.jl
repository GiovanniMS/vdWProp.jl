# creating a structure to van der Waals gases

struct vdWGas
    
    nome::String
    
    Pc::sysP{Float64,EX}
    
    Tc::sysT{Float64,EX}
        
    M::mAmt{Float64,EX,MO}
    
    R_vdw::RAmt{Float64,EX,MA}
    
    ϕ::Float64
    
end

nome(gas::vdWGas) = gas.nome

Pc(gas::vdWGas) = gas.Pc

Tc(gas::vdWGas) = gas.Tc

M(gas::vdWGas) = gas.M

R_vdw(gas::vdWGas) = gas.R

ϕ(gas::vdWGas) = gas.ϕ

α(gas::vdWGas) = (27/64)*(gas.R_vdw^2)*(gas.Tc^2)/gas.Pc

b(gas::vdWGas) = (gas.R_vdw/8)*gas.Tc/gas.Pc

vc(gas::vdWGas) = 3*b(gas)

# Substances

Ar = vdWGas("Argon", P(4860), T(151), (N(1)^-1)*39.948, R(0.2081), (3/2))

Br2 = vdWGas("Bromine", P(10340), T(584), (N(1)^-1)*159.808, R(0.0520), (5/2))

Cl2 = vdWGas("Chlorine", P(7710), T(417), (N(1)^-1)*70.906, R(0.1173), (5/2))

He = vdWGas("Helium", P(228.9945), T(5.21), (N(1)^-1)*4.003, R(2.0769), (3/2))

H2 = vdWGas("Hydrogen", P(1300), T(33.3), (N(1)^-1)*2.016, R(4.1240), (5/2))

Kr = vdWGas("Krypton", P(5500), T(209.4), (N(1)^-1)*83.80, R(0.09921), (3/2))

Ne = vdWGas("Neon", P(2730), T(44.5), (N(1)^-1)*20.183, R(0.4119), (3/2))

N2 = vdWGas("Nytrogen", P(3390), T(126.2), (N(1)^-1)*28.013, R(0.2968), (5/2))

O2 = vdWGas("Oxygen", P(5080), T(154.8), (N(1)^-1)*31.999, R(0.2598), (5/2))

Xe = vdWGas("Xenon", P(5880), T(289.8), (N(1)^-1)*131.30, R(0.06332), (3/2))

