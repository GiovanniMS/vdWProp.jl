using EngThermBase

include("dome.jl")

# creating a structure to van der Waals gases

struct vdWGas
    nome::String
    Pc::sysP{Float64,EX}
    Tc::sysT{Float64,EX}
    α::_Amt{Float64,EX}
    b::vAmt{Float64,EX,MA}
    M::mAmt{Float64,EX,MO}
    Rvdw::RAmt{Float64,EX,MA}
end

nome(gas::vdWGas) = gas.nome
Pc(gas::vdWGas) = gas.Pc
Tc(gas::vdWGas) = gas.Tc
α(gas::vdWGas) = gas.α
b(gas::vdWGas) = gas.b
M(gas::vdWGas) = gas.M
Rvdw(gas::vdWGas) = gas.R
vc(gas::vdWGas) = 3*gas.b

# Gas example

He = vdWGas("Helio", P(228.9945), T(5.21), P(1)*(v(1)^2)*0.21562534694033178, v(1)*0.005945540844366726, (N(1)^-1)*4.003,R(2.0769))

# Now the functions are implemented using the vdWGas as an argument

Tr(Tc::sysT{Float64,EX}, T::sysT{Float64,EX}) = T/Tc

Pr(Pc::sysT{Float64,EX}, P::sysT{Float64,EX}) = P/Pc

vr(vc::vAmt{Float64,EX,MA}, v::vAmt{Float64,EX,MA}) = v/vc

function Pvdw(gas::vdWGas, T::sysT{Float64,EX}, v::vAmt{Float64,EX,MA}) 
    
    if vr2list[findclosest(Tr_sat_list, T, (10^-3))] < v < vr1list[findclosest(Tr_sat_list, T, (10^-3))]
        
        return 8*Tr(Tc(gas),T)/(3*vr(vc(gas),v) - 1) - 3/(vr(vc(gas),v)^2)
    
    else 
        
        return Pr_sat_list[findclosest(Tr_sat_list, T, (10^-3))]
        
    end
    
end