using EngThermBase

using Polynomials

using Roots

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

# Verification if the pair of properties given is inside or outside of the dome. The function is divided in parts of different pairs because the order changes the logic.

function DomeVerification(gas::vdWGas, v::vAmt{Float64,EX,MA}, u::uAmt{Float64,EX,MA})
    
    if findclosest(Domelist(ϕ(gas), "ur1"), ur(u, Pc(gas), vc(gas)), (10^-3)) == -1 && findclosest(Domelist(ϕ(gas), "ur2"), ur(u, Pc(gas), vc(gas)), (10^-3)) == -1
        
        return "out"
        
    else
        
        if vr1list[findclosest(Domelist(ϕ(gas), "ur1"), ur(u, Pc(gas), vc(gas)), (10^-3))] <= vr(vc(gas), v) <= vr2list[findclosest(Domelist(ϕ(gas), "ur2"), ur(u, Pc(gas), vc(gas)), (10^-3))]
            
            return "in"
            
        else
            
            return "out"
            
        end
        
    end
    
end

function DomeVerification(gas::vdWGas, v::vAmt{Float64,EX,MA}, h::hAmt{Float64,EX,MA})
    
    if findclosest(Domelist(ϕ(gas), "hr1"), hr(h, Pc(gas), vc(gas)), (10^-3)) == -1 && findclosest(Domelist(ϕ(gas), "hr2"), hr(h, Pc(gas), vc(gas)), (10^-3)) == -1
        
        return "out"
        
    else
        
        if vr1list[findclosest(Domelist(ϕ(gas), "hr1"), hr(h, Pc(gas), vc(gas)), (10^-3))] <= vr(vc(gas), v) <= vr2list[findclosest(Domelist(ϕ(gas), "hr2"), hr(h, Pc(gas), vc(gas)), (10^-3))]
            
            return "in"
            
        else
            
            return "out"
            
        end
        
    end
    
end

function DomeVerification(gas::vdWGas, v::vAmt{Float64,EX,MA}, s::sAmt{Float64,EX,MA})
    
    if findclosest(Domelist(ϕ(gas), "sr1"), sr(s, Pc(gas), vc(gas), Tc(gas)), (10^-3)) == -1 && findclosest(Domelist(ϕ(gas), "sr2"), sr(s, Pc(gas), vc(gas), Tc(gas)), (10^-3)) == -1
        
        return "out"
        
    else
        
        if vr1list[findclosest(Domelist(ϕ(gas), "sr1"), sr(s, Pc(gas), vc(gas), T(gas)), (10^-3))] <= vr(vc(gas), v) <= vr2list[findclosest(Domelist(ϕ(gas), "sr2"), sr(s, Pc(gas), vc(gas), Tc(gas)), (10^-3))]
            
            return "in"
            
        else
            
            return "out"
            
        end
        
    end
    
end

function DomeVerification(gas::vdWGas, h::hAmt{Float64,EX,MA}, u::uAmt{Float64,EX,MA})
    
    if findclosest(Domelist(ϕ(gas), "ur1"), ur(u, Pc(gas), vc(gas)), (10^-3)) == -1 && findclosest(Domelist(ϕ(gas), "ur2"), ur(u, Pc(gas), vc(gas)), (10^-3)) == -1
        
        return "out"
        
    else
        
        if Domelist(ϕ(gas), "hr1")[findclosest(Domelist(ϕ(gas), "ur1"), ur(u, Pc(gas), vc(gas)), (10^-3))] <= hr(h, Pc(gas), vc(gas)) <= Domelist(ϕ(gas), "hr2")[findclosest(Domelist(ϕ(gas), "ur2"), ur(u, Pc(gas), vc(gas)), (10^-3))]
            
            return "in"
            
        else
            
            return "out"
            
        end
        
    end
    
end

function DomeVerification(gas::vdWGas, h::hAmt{Float64,EX,MA}, s::sAmt{Float64,EX,MA})
    
    if findclosest(Domelist(ϕ(gas), "sr1"), sr(s, Pc(gas), vc(gas), Tc(gas)), (10^-3)) == -1 && findclosest(Domelist(ϕ(gas), "sr2"), sr(s, Pc(gas), vc(gas), Tc(gas)), (10^-3)) == -1
        
        return "out"
        
    else
        
        if Domelist(ϕ(gas), "hr1")[findclosest(Domelist(ϕ(gas), "sr1"), sr(s, Pc(gas), vc(gas), T(gas)), (10^-3))] <= hr(h, Pc(gas), vc(gas)) <= Domelist(ϕ(gas), "hr2")[findclosest(Domelist(ϕ(gas), "sr2"), sr(s, Pc(gas), vc(gas), Tc(gas)), (10^-3))]
            
            return "in"
            
        else
            
            return "out"
            
        end
        
    end
    
end

function DomeVerification(gas::vdWGas, s::sAmt{Float64,EX,MA}, u::uAmt{Float64,EX,MA})
    
    if findclosest(Domelist(ϕ(gas), "ur1"), ur(u, Pc(gas), vc(gas)), (10^-3)) == -1 && findclosest(Domelist(ϕ(gas), "ur2"), ur(u, Pc(gas), vc(gas)), (10^-3)) == -1
        
        return "out"
        
    else
        
        if Domelist(ϕ(gas), "sr1")[findclosest(Domelist(ϕ(gas), "ur1"), ur(u, Pc(gas), vc(gas)), (10^-3))] <= sr(s, Pc(gas), vc(gas), Tc(gas)) <= Domelist(ϕ(gas), "sr2")[findclosest(Domelist(ϕ(gas), "ur2"), ur(u, Pc(gas), vc(gas)), (10^-3))]
            
            return "in"
            
        else
            
            return "out"
            
        end
        
    end
    
end

# function to find the quality of a saturated mixture when P and T are not given

function FindQ(a::BProperty{Float64,EX,MA}, b::BProperty{Float64,EX,MA}, aArray1::Array, aArray2::Array, bArray1::Array, bArray2::Array)

    i = 1
    
    while i <= points
        
        y = ((a - aArray1[i])/(aArray2[i] - aArray1[i])) - ((b - bArray1[i])/(bArray2[i] - bArray1[i]))
        
        if round(y, digits = 3) == 0
            
            return [((a - aArray1[i])/(aArray2[i] - aArray1[i])), aArray1[i]] 
            
        else
            
            yt(j) = ((a - aArray1[j])/(aArray2[j] - aArray1[j])) - ((b - bArray1[j])/(bArray2[j] - bArray1[j]))
            
            j1 = i + 0.5*points
            
            j2 = i + 0.3*points
            
            j3 = i + 0.1*points
            
            j4 = i + 0.05*points
            
            j5 = i + 0.01*points
            
            j6 = i + 0.001*points
            
            if j1 <= points && yt(j1)*y > 0
                    
                i = j1
                
            elseif j2 <= points && yt(j2)*y > 0
                    
                i = j2   
                    
            elseif j3 <= points && yt(j3)*y > 0
                    
                i = j3
                        
            elseif j4 <= points && yt(j4)*y > 0
                    
                i = j4
                            
            elseif j5 <= points && yt(j5)*y > 0
                    
                i = j5
                
            elseif j6 <= points && yt(j6)*y > 0
                    
                i = j6
                
            else
                
                i = i + 1
                
            end
            
        end
        
    end
    
end

# Now the functions are implemented using the vdWGas as an argument

Tr(Tc::sysT{Float64,EX}, T::sysT{Float64,EX}) = T/Tc

Pr(Pc::sysP{Float64,EX}, P::sysP{Float64,EX}) = P/Pc

vr(vc::vAmt{Float64,EX,MA}, v::vAmt{Float64,EX,MA}) = v/vc

ur(u::uAmt{Float64,EX,MA}, Pc::sysP{Float64,EX}, vc::vAmt{Float64,EX,MA}) = u/(Pc*vc)

hr(h::hAmt{Float64,EX,MA}, Pc::sysP{Float64,EX}, vc::vAmt{Float64,EX,MA}) = h/(Pc*vc)

sr(s::sAmt{Float64,EX,MA}, Pc::sysP{Float64,EX}, vc::vAmt{Float64,EX,MA}, Tc::sysT{Float64,EX}) = (Tc*s)/(Pc*vc)

function P_vdw(gas::vdWGas, T::sysT{Float64,EX}, v::vAmt{Float64,EX,MA}) 
    
    if vr1list[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))] <= vr(vc(gas), v) <= vr2list[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))]
        
        return Pc(gas)*(Pr_sat_list[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))])
    
    else 
        
        return Pc(gas)*(8*Tr(Tc(gas),T)/(3*vr(vc(gas),v) - 1) - 3/(vr(vc(gas),v)^2))
    
    end
    
end

function T_vdw(gas::vdWGas, P::sysP{Float64,EX}, v::vAmt{Float64,EX,MA})
    
    if vr1list[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))] <= vr(vc(gas), v) <= vr2list[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))]
        
        return Tc(gas)*(Tr_sat_list[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))])
        
    else
        
        return Tc(gas)*(((Pr(Pc(gas), P)*(3*vr(vc(gas), v) - 1))/8) + (3/(8*(vr(vc(gas), v)^2))))
        
    end
    
end

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, T::sysT{Float64,EX}, Mol::Bool = False)
    
    vrvdw = roots(Poly([(-3/8),(9/8),(-Tr(Tc(gas), T) - (Pr(Pc(gas), P)/8)),(3*Pr(Pc(gas), P)/8)]))
    
    if vr1list[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))] <= vrvdw <= vr2list[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))]
        
        print("Error, T and P can only define a State outside the dome")
        
    else
        
        Mol ? vf = vc(gas)*vrvdw*M(gas) : vf = vc(gas)*vrvdw
        
        return vf
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, s::sAmt{Float64,EX,MA})
    
    if DomeVerification(v,s) == "in"
        
        Q = FindQ(vr(vc(gas), v), sr(s, Pc(gas), vc(gas), Tc(gas)), vrlist1, vrlist2, Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))[1]
        
        vl = FindQ(vr(vc(gas), v), sr(s, Pc(gas), vc(gas), Tc(gas)), vrlist1, vrlist2, Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))[2]
        
        return Tc(gas)*Tr_sat_list[findclosest(vr1list, vr(vc(gas), vl), (10^-3))]
        
    else
        
        return Tc(gas)*(exp((3/(8*ϕ(gas)))*(sr(s,Pc(gas), vc(gas), Tc(gas)) + C2())))*((3*vr(vc(gas), v) - 1)^(-1/ϕ(gas)))
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, u::uAmt{Float64,EX,MA})
    
    if DomeVerification(v,u) == "in"
        
        Q = FindQ(vr(vc(gas), v), sr(s, Pc(gas), vc(gas), Tc(gas)), vrlist1, vrlist2, Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))[1]
        
        vl = FindQ(vr(vc(gas), v), sr(s, Pc(gas), vc(gas), Tc(gas)), vrlist1, vrlist2, Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))[2]
        
        return Tc(gas)*Tr_sat_list[findclosest(vr1list, vr(vc(gas), vl), (10^-3))]
        
    else
        
        return Tc(gas)*(3/8*ϕ(gas))*(ur(u, Pc(gas), vc(gas)) - C1 + (3/vr(vc(gas), v)))
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, h::hAmt{Float64,EX,MA})
    
    if DomeVerification(gas, v, h) == "in"
        
        Q = FindQ(vr(vc(gas), v), sr(s, Pc(gas), vc(gas), Tc(gas)), vrlist1, vrlist2, Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))[1]
        
        vl = FindQ(vr(vc(gas), v), sr(s, Pc(gas), vc(gas), Tc(gas)), vrlist1, vrlist2, Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))[2]
        
        return Tc(gas)*Tr_sat_list[findclosest(vr1list, vr(vc(gas), vl), (10^-3))]
        
    else
        
        return Tc(gas)*(hr(h, Pc(gas), vc(gas)) - C1() + (6/vr(vc(gas), v)))/((8*ϕ(gas)/3) + (8*vr(vc(gas), v)/(3*vr(vc(gas), v) - 1)))
        
    end
    
end

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, u::uAmt{Float64,EX,MA}, Mol::Bool = False)
    
    if Domelist(ϕ(gas), "ur1")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))] <= ur(u, Pc(gas), vc(gas)) <= Domelist(ϕ(gas), "ur2")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))]
        
        Q = (ur(u ,Pc(gas), vc(gas)) - Domelist(ϕ(gas), "ur1")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))])/(Domelist(ϕ(gas), "ur2")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))] - Domelist(ϕ(gas), "ur1")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))])
        
        vrf1 = vr1list[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))] + Q*(vr2list[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))] - vr1list[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        vrf1 = roots(Poly([-3,(9 - (9/ϕ(gas))),(-Pr(Pc(gas), P) - (ur(u, Pc(gas), vc(gas)) - C1())*(3/ϕ(gas))),3*Pr(Pc(gas), P)]))
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, h::hAmt{Float64,EX,MA}, Mol::Bool = False)
        
    if Domelist(ϕ(gas), "hr1")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))] <= hr(h, Pc(gas), vc(gas)) <= Domelist(ϕ(gas), "hr2")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))]
        
        Q = (hr(h ,Pc(gas), vc(gas)) - Domelist(ϕ(gas), "hr1")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))])/(Domelist(ϕ(gas), "hr2")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))] - Domelist(ϕ(gas), "hr1")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))])
        
        vrf1 = vr1list[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))] + Q*(vr2list[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))] - vr1list[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        vrf1 = roots(Poly([(-ϕ(gas)),(3*ϕ(gas) - 3),(-(ϕ(gas)*Pr(Pc(gas), P)/3) - (hr(h, Pc(gas), vc(gas)) - C1())),(ϕ(gas)*Pr(Pc(gas), P) + Pr(Pc(gas), P))]))
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end 

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, s::sAmt{Float64,EX,MA}, Mol::Bool = False)
        
    if Domelist(ϕ(gas), "sr1")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))] <= sr(s, Pc(gas), vc(gas), Tc(gas)) <= Domelist(ϕ(gas), "sr2")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))]
        
        Q = (sr(s ,Pc(gas), vc(gas), Tc(gas)) - Domelist(ϕ(gas), "sr1")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))])/(Domelist(ϕ(gas), "sr2")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))] - Domelist(ϕ(gas), "sr1")[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))])
        
        vrf1 = vr1list[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))] + Q*(vr2list[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))] - vr1list[findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        f(vr) = ((8*exp((3/(8*ϕ(gas)))*(sr(s, Pc(gas), vc(gas), Tc(gas)) + C2())))/((3*vr - 1)^(1 + (1/ϕ(gas))))) - (3/vr^2) - Pr(Pc(gas), P)
        
        vrf1 = find_zero(f,1,Order1())
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end

function v_vdw(gas::vdWGas, T::sysT{Float64,EX}, u::uAmt{Float64,EX,MA}, Mol::Bool = False)
    
    if Domelist(ϕ(gas), "ur1")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))] <= ur(u, Pc(gas), vc(gas)) <= Domelist(ϕ(gas), "ur2")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))]
        
        Q = (ur(u ,Pc(gas), vc(gas)) - Domelist(ϕ(gas), "ur1")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))])/(Domelist(ϕ(gas), "ur2")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))] - Domelist(ϕ(gas), "ur1")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))])
        
        vrf1 = vr1list[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))] + Q*(vr2list[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))] - vr1list[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        vrf1 = -3/(ur(u, Pc(gas), vc(gas)) - C1() - (8*Tr(Tc(gas), T)*ϕ(gas)/3))
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end

function v_vdw(gas::vdWGas, T::sysT{Float64,EX}, s::sAmt{Float64,EX,MA}, Mol::Bool = False)
    
    if Domelist(ϕ(gas), "sr1")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))] <= sr(s, Pc(gas), vc(gas), Tc(gas)) <= Domelist(ϕ(gas), "sr2")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))]
        
        Q = (sr(s ,Pc(gas), vc(gas), Tc(gas)) - Domelist(ϕ(gas), "sr1")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))])/(Domelist(ϕ(gas), "sr2")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))] - Domelist(ϕ(gas), "sr1")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))])
        
        vrf1 = vr1list[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))] + Q*(vr2list[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))] - vr1list[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        vrf1 = (1/3)*(((Tr(Tc(gas), T)/exp((3/(8*ϕ(gas)))*(sr(s, Pc(gas), vc(gas), Tc(gas)) + C2)))^(-ϕ(gas))) + 1)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end

function v_vdw(gas::vdWGas, u::uAmt{Float64,EX,MA}, s::sAmt{Float64,EX,MA}, Mol::Bool = False)
    
    if DomeVerification(gas, s, u) == "in"
        
        Q = FindQ(sr(s, Pc(gas), vc(gas), Tc(gas)), ur(u, Pc(gas), vc(gas)), DomeList(ϕ(gas), "sr1"), DomeList(ϕ(gas), "sr2"), Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))[1]
        
        sl = FindQ(sr(s, Pc(gas), vc(gas), Tc(gas)), ur(u, Pc(gas), vc(gas)), DomeList(ϕ(gas), "sr1"), DomeList(ϕ(gas), "sr2"), Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))[2]
        
        sv = ((s - sl)/Q) + sl
        
        vlr = vr1list[findclosest(DomeList(ϕ(gas), "sr1"), sr(sl, Pc(gas), vc(gas), Tc(gas)), (10^-3))]
        
        vvr = vr2list[findclosest(DomeList(ϕ(gas), "sr2"), sr(sv, Pc(gas), vc(gas), Tc(gas)), (10^-3))]
        
        vrf1 = vlr + Q*(vvr - vlr)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else
        
        f(vr) = (8*ϕ(gas)/3)*(exp((3/(8*ϕ(gas)))*(sr(s, Pc(gas), vc(gas), Tc(gas)) + C2())))*((3*vr - 1)^(-1/ϕ(gas))) - ur(u, Pc(gas), vc(gas)) - (3/vr) + C1()
    
        vrf1 = find_zero(f,0.5,Order1()) #secant method
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2        
        
    end
    
end

function v_vdw(gas::vdWGas, u::uAmt{Float64,EX,MA}, h::hAmt{Float64,EX,MA}, Mol::Bool = False)
    
    if DomeVerification(gas, h, u) == "in"
        
        Q = FindQ(hr(h, Pc(gas), vc(gas)), ur(u, Pc(gas), vc(gas)), DomeList(ϕ(gas), "hr1"), DomeList(ϕ(gas), "hr2"), Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))[1]
        
        hl = FindQ(hr(h, Pc(gas), vc(gas)), ur(u, Pc(gas), vc(gas)), DomeList(ϕ(gas), "hr1"), DomeList(ϕ(gas), "hr2"), Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))[2]
        
        hv = ((h - hl)/Q) + hl
        
        vlr = vr1list[findclosest(DomeList(ϕ(gas), "hr1"), hr(hl, Pc(gas), vc(gas)), (10^-3))]
        
        vvr = vr2list[findclosest(DomeList(ϕ(gas), "hr2"), hr(hv, Pc(gas), vc(gas)), (10^-3))]
        
        vrf1 = vlr + Q*(vvr - vlr)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else
            
        vrf1 = roots(Poly([ϕ(gas),((-ϕ(gas)/3)*(ur(u, Pc(gas), vc(gas)) - C1()) + 3*ϕ(gas) + 3 + (ϕ(gas)/3)*(hr(h, Pc(gas), vc(gas)) - C1()) - 6*ϕ(gas)),((ϕ(gas) + 1)*(ur(u, Pc(gas), vc(gas)) - C1()) - ϕ(gas)*(hr(h, Pc(gas), vc(gas)) - C1()))]))
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2        
        
    end
    
end

function v_vdw(gas::vdWGas, s::sAmt{Float64,EX,MA}, h::hAmt{Float64,EX,MA}, Mol::Bool = False)
    
    if DomeVerification(gas, h, s) == "in"
        
        Q = FindQ(hr(h, Pc(gas), vc(gas)), sr(s, Pc(gas), vc(gas), Tc(gas)), DomeList(ϕ(gas), "hr1"), DomeList(ϕ(gas), "hr2"), Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))[1]
        
        hl = FindQ(hr(h, Pc(gas), vc(gas)), sr(s, Pc(gas), vc(gas), Tc(gas)), DomeList(ϕ(gas), "hr1"), DomeList(ϕ(gas), "hr2"), Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))[2]
        
        hv = ((h - hl)/Q) + hl
        
        vlr = vr1list[findclosest(DomeList(ϕ(gas), "hr1"), hr(hl, Pc(gas), vc(gas)), (10^-3))]
        
        vvr = vr2list[findclosest(DomeList(ϕ(gas), "hr2"), hr(hv, Pc(gas), vc(gas)), (10^-3))]
        
        vrf1 = vlr + Q*(vvr - vlr)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else
            
        f(vr) = ((8*ϕ(gas)/3) + (8*vr/(3*vr - 1)))*exp((3/(8*ϕ(gas)))*(sr(s, Pc(gas), vc(gas), Tc(gas)) + C2()))*((3*vr - 1)^(-1/ϕ(gas))) - hr(h, Pc(gas), vc(gas)) + C1() - (6/vr)
    
        vrf1 = find_zero(f,0.5,Order1()) #secant method
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2        
        
    end
    
end

function v_vdw(gas::vdWGas, T::sysT{Float64,EX}, h::hAmt{Float64,EX,MA}, Mol::Bool = False)
    
    if Domelist(ϕ(gas), "hr1")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))] <= hr(h, Pc(gas), vc(gas)) <= Domelist(ϕ(gas), "hr2")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))]
        
        Q = (hr(h ,Pc(gas), vc(gas)) - Domelist(ϕ(gas), "hr1")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))])/(Domelist(ϕ(gas), "hr2")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))] - Domelist(ϕ(gas), "hr1")[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))])
        
        vrf1 = vr1list[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))] + Q*(vr2list[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))] - vr1list[findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        vrf1 = roots(Poly([6, ((-8*ϕ(gas)*Tr(Tc(gas), T)/3) + hr(h, Pc(gas), vc(gas)) - C1() - 18), (8*Tr(Tc(gas), T)*(ϕ(gas) + 1) - 3*hr(h, Pc(gas), vc(gas)) + 3*C1())]))
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end