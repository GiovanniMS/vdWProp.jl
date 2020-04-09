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

# function to find the quality of a saturated mixture when P and T are not given, and also find if the pair is inside or outside of the dome.

function FindQ(a::BProperty{Float64,EX,MA}, b::BProperty{Float64,EX,MA}, aArray1::Array, aArray2::Array, bArray1::Array, bArray2::Array)

    i = 1
    
    while i <= points
        
        typeof(a) == vAmt{Float64,EX,MA} ? A = [v(aArray1[i]), v(aArray2[i])] : 
    
        typeof(a) == uAmt{Float64,EX,MA} ? A = [u(aArray1[i]), u(aArray2[i])] : 

        typeof(a) == hAmt{Float64,EX,MA} ? A = [h(aArray1[i]), h(aArray2[i])] : 

        typeof(a) == sAmt{Float64,EX,MA} ? A = [s(aArray1[i]), s(aArray2[i])] : z = 1

        typeof(b) == vAmt{Float64,EX,MA} ? B = [v(aArray1[i]), v(aArray2[i])] : 

        typeof(b) == uAmt{Float64,EX,MA} ? B = [u(aArray1[i]), u(aArray2[i])] :

        typeof(b) == hAmt{Float64,EX,MA} ? B = [h(aArray1[i]), h(aArray2[i])] :

        typeof(b) == sAmt{Float64,EX,MA} ? B = [s(aArray1[i]), s(aArray2[i])] : z = 1
        
        y = ((a - A[1])/(A[2] - A[1])) - ((b - B[1])/(B[2] - B[1]))
        
        if y < AMT(0.001)
            
            return [((a - A[1])/(A[2] - A[1])), A[1]] 
            
        else
            
            function yt(j) 
                
                typeof(a) == vAmt{Float64,EX,MA} ? A = [v(aArray1[j]), v(aArray2[j])] : 
    
                typeof(a) == uAmt{Float64,EX,MA} ? A = [u(aArray1[j]), u(aArray2[j])] : 

                typeof(a) == hAmt{Float64,EX,MA} ? A = [h(aArray1[j]), h(aArray2[j])] : 

                typeof(a) == sAmt{Float64,EX,MA} ? A = [s(aArray1[j]), s(aArray2[j])] : z = 1

                typeof(b) == vAmt{Float64,EX,MA} ? B = [v(aArray1[j]), v(aArray2[j])] : 

                typeof(b) == uAmt{Float64,EX,MA} ? B = [u(aArray1[j]), u(aArray2[j])] :

                typeof(b) == hAmt{Float64,EX,MA} ? B = [h(aArray1[j]), h(aArray2[j])] :

                typeof(b) == sAmt{Float64,EX,MA} ? B = [s(aArray1[j]), s(aArray2[j])] : z = 1
                
                return ((a - A[1])/(A[2] - A[1])) - ((b - B[1])/(B[2] - B[1]))
                
            end
            
            j1 = Integer(i + 0.5*points)
            
            j2 = Integer(i + 0.3*points)
            
            j3 = Integer(i + 0.1*points)
            
            j4 = Integer(i + 0.05*points)
            
            j5 = Integer(i + 0.01*points)
            
            j6 = Integer(i + 0.001*points)
            
            if j1 <= points && yt(j1)*y > AMT(0)
                    
                i = j1
                
            elseif j2 <= points && yt(j2)*y > AMT(0)
                    
                i = j2   
                    
            elseif j3 <= points && yt(j3)*y > AMT(0)
                    
                i = j3
                        
            elseif j4 <= points && yt(j4)*y > AMT(0)
                    
                i = j4
                            
            elseif j5 <= points && yt(j5)*y > AMT(0)
                    
                i = j5
                
            elseif j6 <= points && yt(j6)*y > AMT(0)
                    
                i = j6
                
            else
                
                i = i + 1
                
            end
            
        end
        
    end
    
    return "out"
    
end

#Function to find the real root in the array that will be the result for the 3 degree polinomial

function ImVerification(a::Array)
    
    i = 1
    
    while i > 0 && i < 4
        
        imag(a[i]) == 0 ? break : i = i + 1
        
    end
    
    return Float64(a[i])
    
end

# Now the functions are implemented using the vdWGas as an argument

Tr(Tc::sysT{Float64,EX}, T::sysT{Float64,EX}) = T/Tc

Pr(Pc::sysP{Float64,EX}, P::sysP{Float64,EX}) = P/Pc

vr(vc::vAmt{Float64,EX,MA}, v::vAmt{Float64,EX,MA}) = v/vc

ur(u::uAmt{Float64,EX,MA}, Pc::sysP{Float64,EX}, vc::vAmt{Float64,EX,MA}) = u/(Pc*vc)

hr(h::hAmt{Float64,EX,MA}, Pc::sysP{Float64,EX}, vc::vAmt{Float64,EX,MA}) = h/(Pc*vc)

sr(s::sAmt{Float64,EX,MA}, Pc::sysP{Float64,EX}, vc::vAmt{Float64,EX,MA}, Tc::sysT{Float64,EX}) = (Tc*s)/(Pc*vc)

function P_vdw(gas::vdWGas, T::sysT{Float64,EX}, v::vAmt{Float64,EX,MA}) 
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3)) # Array position for the saturated fluid and gas 
    
    if vr1list[SatP] <= vr(vc(gas), v) <= vr2list[SatP]
        
        return Pc(gas)*(Pr_sat_list[SatP])
    
    else 
        
        return Pc(gas)*(8*Tr(Tc(gas),T)/(3*vr(vc(gas),v) - 1) - 3/(vr(vc(gas),v)^2))
    
    end
    
end

function T_vdw(gas::vdWGas, P::sysP{Float64,EX}, v::vAmt{Float64,EX,MA})
    
    SatP = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))
    
    if vr1list[SatP] <= vr(vc(gas), v) <= vr2list[SatP]
        
        return Tc(gas)*(Tr_sat_list[SatP])
        
    else
        
        return Tc(gas)*(((Pr(Pc(gas), P)*(3*vr(vc(gas), v) - 1))/8) + (3/(8*(vr(vc(gas), v)^2))))
        
    end
    
end

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, T::sysT{Float64,EX}, Mol::Bool = False)
    
    SatP = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))
    
    vrvdw = roots(Poly([(-3/8),(9/8),(-Tr(Tc(gas), T) - (Pr(Pc(gas), P)/8)),(3*Pr(Pc(gas), P)/8)]))
    
    if vr1list[SatP] <= vrvdw <= vr2list[SatP]
        
        print("Error, T and P can only define a State outside the dome")
        
    else
        
        Mol ? vf = vc(gas)*vrvdw*M(gas) : vf = vc(gas)*vrvdw
        
        return vf
        
    end
    
end

function s_vdw(gas::vdWGas, T::sysT{Float64,EX}, v::vAmt{Float64,EX,MA})
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))
    
    if vr1list[SatP] <= vr(vc(gas), v) <= vr2list[SatP]
        
        Q = (vr(vc(gas), v) - vr1list[SatP])/(vr2list[SatP] - vr1list[SatP])
        
        srf1 = Domelist(ϕ(gas), "sr1")[SatP] + Q*(Domelist(ϕ(gas), "sr2")[SatP] - Domelist(ϕ(gas), "sr1")[SatP])
        
        Mol ? srf2 = srf1*M(gas) : srf2 = srf1 
        
        return Pc(gas)*vc(gas)*srf2/Tc(gas)
        
    else 
        
        srf1 = (8/3)*log(3*vr(vc(gas), v) - 1) + (8*ϕ(gas)/3)*log(Tr(Tc(gas), T)) - C2()
        
        Mol ? srf2 = srf1*M(gas) : srf2 = srf1 
        
        return Pc(gas)*vc(gas)*srf2/Tc(gas)
        
    end
    
end

function u_vdw(gas::vdWGas, T::sysT{Float64,EX}, v::vAmt{Float64,EX,MA})
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))
    
    if vr1list[SatP] <= vr(vc(gas), v) <= vr2list[SatP]
        
        Q = (vr(vc(gas), v) - vr1list[SatP])/(vr2list[SatP] - vr1list[SatP])
        
        urf1 = Domelist(ϕ(gas), "ur1")[SatP] + Q*(Domelist(ϕ(gas), "ur2")[SatP] - Domelist(ϕ(gas), "ur1")[SatP])
        
        Mol ? urf2 = urf1*M(gas) : urf2 = urf1 
        
        return Pc(gas)*vc(gas)*urf2
        
    else 
        
        urf1 = C1() + (8*Tr(Tc(gas), T)*ϕ(gas)/3) - (3/vr(vc(gas), v))
        
        Mol ? urf2 = urf1*M(gas) : urf2 = urf1 
        
        return Pc(gas)*vc(gas)*urf2
        
    end
    
end

function h_vdw(gas::vdWGas, T::sysT{Float64,EX}, v::vAmt{Float64,EX,MA})
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))
    
    if vr1list[SatP] <= vr(vc(gas), v) <= vr2list[SatP]
        
        Q = (vr(vc(gas), v) - vr1list[SatP])/(vr2list[SatP] - vr1list[SatP])
        
        hrf1 = Domelist(ϕ(gas), "hr1")[SatP] + Q*(Domelist(ϕ(gas), "hr2")[SatP] - Domelist(ϕ(gas), "hr1")[SatP])
        
        Mol ? hrf2 = hrf1*M(gas) : hrf2 = hrf1 
        
        return Pc(gas)*vc(gas)*hrf2
        
    else 
        
        hrf1 = C1() + (8*Tr(Tc(gas), T)*ϕ(gas)/3) + (8*Tr(Tc(gas), T)*vr(vc(gas), v)/(3*vr(vc(gas), v) - 1)) - (6/vr(vc(gas), v))
        
        Mol ? hrf2 = hrf1*M(gas) : hrf2 = hrf1 
        
        return Pc(gas)*vc(gas)*hrf2
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, s::sAmt{Float64,EX,MA})
    
    FQ = FindQ(vr(vc(gas), v), sr(s, Pc(gas), vc(gas), Tc(gas)), vrlist1, vrlist2, Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))
    
    if FQ == "out"
        
        return Tc(gas)*(exp((3/(8*ϕ(gas)))*(sr(s,Pc(gas), vc(gas), Tc(gas)) + C2())))*((3*vr(vc(gas), v) - 1)^(-1/ϕ(gas)))
        
    else
        
        vl = FQ[2]
        
        return Tc(gas)*Tr_sat_list[findclosest(vr1list, vr(vc(gas), vl), (10^-3))]
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, u::uAmt{Float64,EX,MA})
    
    FQ = FindQ(vr(vc(gas), v), sr(s, Pc(gas), vc(gas), Tc(gas)), vrlist1, vrlist2, Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))
    
    if FQ == "out"
        
        return Tc(gas)*(3/8*ϕ(gas))*(ur(u, Pc(gas), vc(gas)) - C1 + (3/vr(vc(gas), v)))
        
    else
        
        vl = FQ[2]
        
        return Tc(gas)*Tr_sat_list[findclosest(vr1list, vr(vc(gas), vl), (10^-3))]
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, h::hAmt{Float64,EX,MA})
    
    FQ = FindQ(vr(vc(gas), v), sr(s, Pc(gas), vc(gas), Tc(gas)), vrlist1, vrlist2, Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2")) 
    
    if FQ == "out"
        
        return Tc(gas)*(hr(h, Pc(gas), vc(gas)) - C1() + (6/vr(vc(gas), v)))/((8*ϕ(gas)/3) + (8*vr(vc(gas), v)/(3*vr(vc(gas), v) - 1)))
        
    else
        
        vl = FQ[2]
        
        return Tc(gas)*Tr_sat_list[findclosest(vr1list, vr(vc(gas), vl), (10^-3))]
        
    end
    
end

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, u::uAmt{Float64,EX,MA}, Mol::Bool = False)
    
    SatP = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))
    
    if Domelist(ϕ(gas), "ur1")[SatP] <= ur(u, Pc(gas), vc(gas)) <= Domelist(ϕ(gas), "ur2")[SatP]
        
        Q = (ur(u ,Pc(gas), vc(gas)) - Domelist(ϕ(gas), "ur1")[SatP])/(Domelist(ϕ(gas), "ur2")[SatP] - Domelist(ϕ(gas), "ur1")[SatP])
        
        vrf1 = vr1list[SatP] + Q*(vr2list[SatP] - vr1list[SatP])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        vrf1 = roots(Poly([-3,(9 - (9/ϕ(gas))),(-Pr(Pc(gas), P) - (ur(u, Pc(gas), vc(gas)) - C1())*(3/ϕ(gas))),3*Pr(Pc(gas), P)]))
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, h::hAmt{Float64,EX,MA}, Mol::Bool = False)
    
    SatP = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))
        
    if Domelist(ϕ(gas), "hr1")[SatP] <= hr(h, Pc(gas), vc(gas)) <= Domelist(ϕ(gas), "hr2")[SatP]
        
        Q = (hr(h ,Pc(gas), vc(gas)) - Domelist(ϕ(gas), "hr1")[SatP])/(Domelist(ϕ(gas), "hr2")[SatP] - Domelist(ϕ(gas), "hr1")[SatP])
        
        vrf1 = vr1list[SatP] + Q*(vr2list[SatP] - vr1list[SatP])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        vrf1 = roots(Poly([(-ϕ(gas)),(3*ϕ(gas) - 3),(-(ϕ(gas)*Pr(Pc(gas), P)/3) - (hr(h, Pc(gas), vc(gas)) - C1())),(ϕ(gas)*Pr(Pc(gas), P) + Pr(Pc(gas), P))]))
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end 

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, s::sAmt{Float64,EX,MA}, Mol::Bool = False)
    
    SatP = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))
        
    if Domelist(ϕ(gas), "sr1")[SatP] <= sr(s, Pc(gas), vc(gas), Tc(gas)) <= Domelist(ϕ(gas), "sr2")[SatP]
        
        Q = (sr(s ,Pc(gas), vc(gas), Tc(gas)) - Domelist(ϕ(gas), "sr1")[SatP])/(Domelist(ϕ(gas), "sr2")[SatP] - Domelist(ϕ(gas), "sr1")[SatP])
        
        vrf1 = vr1list[SatP] + Q*(vr2list[SatP] - vr1list[SatP])
        
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
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))
    
    if Domelist(ϕ(gas), "ur1")[SatP] <= ur(u, Pc(gas), vc(gas)) <= Domelist(ϕ(gas), "ur2")[SatP]
        
        Q = (ur(u ,Pc(gas), vc(gas)) - Domelist(ϕ(gas), "ur1")[SatP])/(Domelist(ϕ(gas), "ur2")[SatP] - Domelist(ϕ(gas), "ur1")[SatP])
        
        vrf1 = vr1list[SatP] + Q*(vr2list[SatP] - vr1list[SatP])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        vrf1 = -3/(ur(u, Pc(gas), vc(gas)) - C1() - (8*Tr(Tc(gas), T)*ϕ(gas)/3))
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end

function v_vdw(gas::vdWGas, T::sysT{Float64,EX}, s::sAmt{Float64,EX,MA}, Mol::Bool = False)
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))
    
    if Domelist(ϕ(gas), "sr1")[SatP] <= sr(s, Pc(gas), vc(gas), Tc(gas)) <= Domelist(ϕ(gas), "sr2")[SatP]
        
        Q = (sr(s ,Pc(gas), vc(gas), Tc(gas)) - Domelist(ϕ(gas), "sr1")[SatP])/(Domelist(ϕ(gas), "sr2")[SatP] - Domelist(ϕ(gas), "sr1")[SatP])
        
        vrf1 = vr1list[SatP] + Q*(vr2list[SatP] - vr1list[SatP])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        vrf1 = (1/3)*(((Tr(Tc(gas), T)/exp((3/(8*ϕ(gas)))*(sr(s, Pc(gas), vc(gas), Tc(gas)) + C2)))^(-ϕ(gas))) + 1)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end

function v_vdw(gas::vdWGas, u::uAmt{Float64,EX,MA}, s::sAmt{Float64,EX,MA}, Mol::Bool = False)
    
    FQ = FindQ(sr(s, Pc(gas), vc(gas), Tc(gas)), ur(u, Pc(gas), vc(gas)), DomeList(ϕ(gas), "sr1"), DomeList(ϕ(gas), "sr2"), Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))
    
    if FQ == "out"
        
        f(vr) = (8*ϕ(gas)/3)*(exp((3/(8*ϕ(gas)))*(sr(s, Pc(gas), vc(gas), Tc(gas)) + C2())))*((3*vr - 1)^(-1/ϕ(gas))) - ur(u, Pc(gas), vc(gas)) - (3/vr) + C1()
    
        vrf1 = find_zero(f,0.5,Order1()) #metodo da secante
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2         
        
    else
        
        Q = FQ[1]
        
        sl = FQ[2]
        
        sv = ((s - sl)/Q) + sl
        
        vlr = vr1list[findclosest(DomeList(ϕ(gas), "sr1"), sr(sl, Pc(gas), vc(gas), Tc(gas)), (10^-3))]
        
        vvr = vr2list[findclosest(DomeList(ϕ(gas), "sr2"), sr(sv, Pc(gas), vc(gas), Tc(gas)), (10^-3))]
        
        vrf1 = vlr + Q*(vvr - vlr)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2    
        
    end
    
end

function v_vdw(gas::vdWGas, u::uAmt{Float64,EX,MA}, h::hAmt{Float64,EX,MA}, Mol::Bool = False)
    
    FQ = FindQ(hr(h, Pc(gas), vc(gas)), ur(u, Pc(gas), vc(gas)), DomeList(ϕ(gas), "hr1"), DomeList(ϕ(gas), "hr2"), Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))
    
    if FQ == "out"
        
        vrf1 = roots(Poly([ϕ(gas),((-ϕ(gas)/3)*(ur(u, Pc(gas), vc(gas)) - C1()) + 3*ϕ(gas) + 3 + (ϕ(gas)/3)*(hr(h, Pc(gas), vc(gas)) - C1()) - 6*ϕ(gas)),((ϕ(gas) + 1)*(ur(u, Pc(gas), vc(gas)) - C1()) - ϕ(gas)*(hr(h, Pc(gas), vc(gas)) - C1()))]))
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2        
        
    else
        
        Q = FQ[1]
        
        hl = FQ[2]
        
        hv = ((h - hl)/Q) + hl
        
        vlr = vr1list[findclosest(DomeList(ϕ(gas), "hr1"), hr(hl, Pc(gas), vc(gas)), (10^-3))]
        
        vvr = vr2list[findclosest(DomeList(ϕ(gas), "hr2"), hr(hv, Pc(gas), vc(gas)), (10^-3))]
        
        vrf1 = vlr + Q*(vvr - vlr)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2   
        
    end
    
end

function v_vdw(gas::vdWGas, s::sAmt{Float64,EX,MA}, h::hAmt{Float64,EX,MA}, Mol::Bool = False)
    
    FQ = FindQ(hr(h, Pc(gas), vc(gas)), sr(s, Pc(gas), vc(gas), Tc(gas)), DomeList(ϕ(gas), "hr1"), DomeList(ϕ(gas), "hr2"), Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))
    
    if FQ == "out"
        
        f(vr) = ((8*ϕ(gas)/3) + (8*vr/(3*vr - 1)))*exp((3/(8*ϕ(gas)))*(sr(s, Pc(gas), vc(gas), Tc(gas)) + C2()))*((3*vr - 1)^(-1/ϕ(gas))) - hr(h, Pc(gas), vc(gas)) + C1() - (6/vr)
    
        vrf1 = find_zero(f,0.5,Order1()) #metodo da secante
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2 
        
    else
        
        Q = FQ[1]
        
        hl = FQ[2]
        
        hv = ((h - hl)/Q) + hl
        
        vlr = vr1list[findclosest(DomeList(ϕ(gas), "hr1"), hr(hl, Pc(gas), vc(gas)), (10^-3))]
        
        vvr = vr2list[findclosest(DomeList(ϕ(gas), "hr2"), hr(hv, Pc(gas), vc(gas)), (10^-3))]
        
        vrf1 = vlr + Q*(vvr - vlr)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2      
        
    end
    
end

function v_vdw(gas::vdWGas, T::sysT{Float64,EX}, h::hAmt{Float64,EX,MA}, Mol::Bool = False)
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))
    
    if Domelist(ϕ(gas), "hr1")[SatP] <= hr(h, Pc(gas), vc(gas)) <= Domelist(ϕ(gas), "hr2")[SatP]
        
        Q = (hr(h ,Pc(gas), vc(gas)) - Domelist(ϕ(gas), "hr1")[SatP])/(Domelist(ϕ(gas), "hr2")[SatP] - Domelist(ϕ(gas), "hr1")[SatP])
        
        vrf1 = vr1list[SatP] + Q*(vr2list[SatP] - vr1list[SatP])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        vrf1 = roots(Poly([6, ((-8*ϕ(gas)*Tr(Tc(gas), T)/3) + hr(h, Pc(gas), vc(gas)) - C1() - 18), (8*Tr(Tc(gas), T)*(ϕ(gas) + 1) - 3*hr(h, Pc(gas), vc(gas)) + 3*C1())]))
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end

# with all the possible pairs implemented, the next step is implement a function that gets all the six properties when a random pair is given

function State(gas::vdWGas, a::AMOUNTS{Float64,EX}, b::AMOUNTS{Float64,EX}, Mol::Bool = False)
    
    ta = typeof(a)
    
    tb = typeof(b)
    
    if (ta == sysP{Float64,EX} && tb == sysT{Float64,EX}) || (tb == sysP{Float64,EX} && ta == sysT{Float64,EX})
        
        ta == sysP{Float64,EX} ? P = a : P = b
        
        tb == sysT{Float64,EX} ? T = b : T = a
        
        v = v_vdw(gas, P, T)
        
        u = u_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)
        
        s = s_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    elseif (ta == sysP{Float64,EX} && tb == vAmt{Float64,EX,MA}) || (tb == sysP{Float64,EX} && ta == vAmt{Float64,EX,MA})
        
        ta == sysP{Float64,EX} ? P = a : P = b
        
        tb == vAmt{Float64,EX,MA} ? v = b : v = a
        
        T = T_vdw(gas, P, v)
        
        u = u_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)
        
        s = s_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    elseif (ta == sysP{Float64,EX} && tb == uAmt{Float64,EX,MA}) || (tb == sysP{Float64,EX} && ta == uAmt{Float64,EX,MA})
        
        ta == sysP{Float64,EX} ? P = a : P = b
        
        tb == uAmt{Float64,EX,MA} ? u = b : u = a
        
        v = v_vdw(gas, P, u)
        
        T = T_vdw(gas, P, v)
        
        h = h_vdw(gas, T, v)
        
        s = s_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    elseif (ta == sysP{Float64,EX} && tb == hAmt{Float64,EX,MA}) || (tb == sysP{Float64,EX} && ta == hAmt{Float64,EX,MA})
        
        ta == sysP{Float64,EX} ? P = a : P = b
        
        tb == hAmt{Float64,EX,MA} ? h = b : h = a
        
        v = v_vdw(gas, P, h)
        
        T = T_vdw(gas, P, v)
        
        u = u_vdw(gas, T, v)
        
        s = s_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    elseif (ta == sysP{Float64,EX} && tb == sAmt{Float64,EX,MA}) || (tb == sysP{Float64,EX} && ta == sAmt{Float64,EX,MA})
        
        ta == sysP{Float64,EX} ? P = a : P = b
        
        tb == sAmt{Float64,EX,MA} ? s = b : s = a
        
        v = v_vdw(gas, P, s)
        
        T = T_vdw(gas, P, v)
        
        u = u_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    elseif (ta == sysT{Float64,EX} && tb == vAmt{Float64,EX,MA}) || (tb == sysT{Float64,EX} && ta == vAmt{Float64,EX,MA})
        
        ta == sysT{Float64,EX} ? T = a : T = b
        
        tb == vAmt{Float64,EX,MA} ? v = b : v = a
        
        P = P_vdw(gas, T, v)
        
        u = u_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)
        
        s = s_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    elseif (ta == sysT{Float64,EX} && tb == sAmt{Float64,EX,MA}) || (tb == sysT{Float64,EX} && ta == sAmt{Float64,EX,MA})
        
        ta == sysT{Float64,EX} ? T = a : T = b
        
        tb == sAmt{Float64,EX,MA} ? s = b : s = a
        
        v = v_vdw(gas, T, s)
        
        P = P_vdw(gas, T, v)
        
        u = u_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    elseif (ta == sysT{Float64,EX} && tb == hAmt{Float64,EX,MA}) || (tb == sysT{Float64,EX} && ta == hAmt{Float64,EX,MA})
        
        ta == sysT{Float64,EX} ? T = a : T = b
        
        tb == hAmt{Float64,EX,MA} ? h = b : h = a
        
        v = v_vdw(gas, T, h)
        
        P = P_vdw(gas, T, v)
        
        u = u_vdw(gas, T, v)
        
        s = s_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    elseif (ta == sysT{Float64,EX} && tb == uAmt{Float64,EX,MA}) || (tb == sysT{Float64,EX} && ta == uAmt{Float64,EX,MA})
        
        ta == sysT{Float64,EX} ? T = a : T = b
        
        tb == uAmt{Float64,EX,MA} ? u = b : u = a
        
        v = v_vdw(gas, T, u)
        
        P = P_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)
        
        s = s_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    elseif (ta == vAmt{Float64,EX,MA} && tb == sAmt{Float64,EX,MA}) || (tb == vAmt{Float64,EX,MA} && ta == sAmt{Float64,EX,MA})
        
        ta == vAmt{Float64,EX,MA} ? v = a : v = b
        
        tb == sAmt{Float64,EX,MA} ? s = b : s = a
        
        T = T_vdw(gas, v, s)
                
        P = P_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)
        
        u = u_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    elseif (ta == vAmt{Float64,EX,MA} && tb == hAmt{Float64,EX,MA}) || (tb == vAmt{Float64,EX,MA} && ta == hAmt{Float64,EX,MA})
        
        ta == vAmt{Float64,EX,MA} ? v = a : v = b
        
        tb == hAmt{Float64,EX,MA} ? h = b : h = a
        
        T = T_vdw(gas, v, h)
                
        P = P_vdw(gas, T, v)
        
        s = s_vdw(gas, T, v)
        
        u = u_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    elseif (ta == vAmt{Float64,EX,MA} && tb == uAmt{Float64,EX,MA}) || (tb == vAmt{Float64,EX,MA} && ta == uAmt{Float64,EX,MA})
        
        ta == vAmt{Float64,EX,MA} ? v = a : v = b
        
        tb == uAmt{Float64,EX,MA} ? u = b : u = a
        
        T = T_vdw(gas, v, u)
                
        P = P_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)
        
        s = s_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    elseif (ta == sAmt{Float64,EX,MA} && tb == hAmt{Float64,EX,MA}) || (tb == sAmt{Float64,EX,MA} && ta == hAmt{Float64,EX,MA})
        
        ta == sAmt{Float64,EX,MA} ? s = a : s = b
        
        tb == hAmt{Float64,EX,MA} ? h = b : h = a
        
        v = v_vdw(gas, s, h)
        
        T = T_vdw(gas, v, s)
                
        P = P_vdw(gas, T, v)
        
        u = u_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    elseif (ta == sAmt{Float64,EX,MA} && tb == uAmt{Float64,EX,MA}) || (tb == sAmt{Float64,EX,MA} && ta == uAmt{Float64,EX,MA})
        
        ta == sAmt{Float64,EX,MA} ? s = a : s = b
        
        tb == uAmt{Float64,EX,MA} ? u = b : u = a
        
        v = v_vdw(gas, u, s)
        
        T = T_vdw(gas, v, s)
                
        P = P_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    elseif (ta == hAmt{Float64,EX,MA} && tb == uAmt{Float64,EX,MA}) || (tb == hAmt{Float64,EX,MA} && ta == uAmt{Float64,EX,MA})
        
        ta == hAmt{Float64,EX,MA} ? h = a : h = b
        
        tb == uAmt{Float64,EX,MA} ? u = b : u = a
        
        v = v_vdw(gas, u, h)
        
        T = T_vdw(gas, v, h)
                
        P = P_vdw(gas, T, v)
        
        s = s_vdw(gas, T, v)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas)] : St = [P, T, v, u, h, s]
        
        return St
        
    else
        
        println("ERROR, the arguments needs to be properties between P,T,v,h,u,s and the base for the intensive ones needs to be mass.")
        
    end
    
end