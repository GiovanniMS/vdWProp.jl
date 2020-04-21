module vdWProp

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

function FindQ(a::_Amt{Float64,EX}, b::_Amt{Float64,EX}, aArray1::Array, aArray2::Array, bArray1::Array, bArray2::Array)

    i = 1
    
    while i <= points
        
        Q = ((a - AMT(aArray1[i]))/(AMT(aArray2[i] - aArray1[i])))
        
        y = ((a - AMT(aArray1[i]))/(AMT(aArray2[i] - aArray1[i]))) - ((b - AMT(bArray1[i]))/(AMT(bArray2[i] - bArray1[i])))
        
        if abs(AMTConvert(y)) < 0.001 && AMT(0) < Q < AMT(1) 
            
            return [Q, aArray1[i]] 
            
        else
            
            yt(j) = ((a - AMT(aArray1[j]))/(AMT(aArray2[j] - aArray1[j]))) - ((b - AMT(bArray1[j]))/(AMT(bArray2[j] - bArray1[j])))
            
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

#Function to convert a AMT into a Float64

function AMTConvert(amt::_Amt{Float64,EX})
    
    x = amt + AMT(0.0001)
    
    z = "$x"
    
    a = String(SubString(z, 10:14))
    
    b = parse(Float64, a)
    
    ba = AMT(b)
    
    d = amt - ba
    
    z2 = "$d"
    
    a2 = String(SubString(z2, 10:14))
    
    b2 = parse(Float64, a2)
    
    d < AMT(0) ? r = (b - b2) : r = (b + b2)
    
    return r
    
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
    
    if SatP > 0 && AMT(vr1list[SatP]) < vr(vc(gas), v) < AMT(vr2list[SatP])
        
        return Pc(gas)*(Pr_sat_list[SatP])
    
    else 
        
        return Pc(gas)*(8*Tr(Tc(gas),T)/(3*vr(vc(gas),v) - AMT(1)) - AMT(3)/(vr(vc(gas),v)^2))
    
    end
    
end

function T_vdw(gas::vdWGas, P::sysP{Float64,EX}, v::vAmt{Float64,EX,MA})
    
    SatP = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))
    
    if SatP > 0 && AMT(vr1list[SatP]) < vr(vc(gas), v) < AMT(vr2list[SatP])
        
        return Tc(gas)*(Tr_sat_list[SatP])
        
    else
        
        return Tc(gas)*(((Pr(Pc(gas), P)*(3*vr(vc(gas), v) - AMT(1)))/8) + (AMT(3)/(8*(vr(vc(gas), v)^2))))
        
    end
    
end

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, T::sysT{Float64,EX}, Mol::Bool = false)
    
    SatP = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))

    Proots = AMTConvert(Pr(Pc(gas), P))
    
    Troots = AMTConvert(Tr(Tc(gas), T))
    
    vrvdw = roots(Poly([(-3/8),(9/8),(-Troots - (Proots/8)),(3*Proots/8)]))
    
    vrvdw = ImVerification(vrvdw)
        
    if SatP > 0 && vr1list[SatP] < vrvdw < vr2list[SatP]
        
        print("Error, T and P can only define a State outside the dome")
        
    else
        
        Mol ? vf = vc(gas)*vrvdw*M(gas) : vf = vc(gas)*vrvdw
        
        return vf
        
    end
    
end

function s_vdw(gas::vdWGas, T::sysT{Float64,EX}, v::vAmt{Float64,EX,MA}, Mol::Bool = false)
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))
    
    if SatP > 0 && AMT(vr1list[SatP]) < vr(vc(gas), v) < AMT(vr2list[SatP])
        
        vr1 = AMT(vr1list[SatP])
    
        vr2 = AMT(vr2list[SatP])
        
        Q = (vr(vc(gas), v) - vr1)/(vr2 - vr1)
        
        srf1 = AMT(Domelist(ϕ(gas), "sr1")[SatP]) + Q*(AMT(Domelist(ϕ(gas), "sr2")[SatP]) - AMT(Domelist(ϕ(gas), "sr1")[SatP]))
        
        Mol ? srf2 = srf1*M(gas) : srf2 = srf1 
        
        return Pc(gas)*vc(gas)*srf2/Tc(gas)
        
    else 
        
        srf1 = (8/3)*log(3*vr(vc(gas), v) - AMT(1)) + (8*ϕ(gas)/3)*log(Tr(Tc(gas), T)) - AMT(C2())
        
        Mol ? srf2 = srf1*M(gas) : srf2 = srf1 
        
        return Pc(gas)*vc(gas)*srf2/Tc(gas)
        
    end
    
end

function u_vdw(gas::vdWGas, T::sysT{Float64,EX}, v::vAmt{Float64,EX,MA}, Mol::Bool = false)
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))
    
    if SatP > 0 && AMT(vr1list[SatP]) < vr(vc(gas), v) < AMT(vr2list[SatP])
        
        vr1 = AMT(vr1list[SatP])
    
        vr2 = AMT(vr2list[SatP])
        
        Q = (vr(vc(gas), v) - vr1)/(vr2 - vr1)
        
        urf1 = AMT(Domelist(ϕ(gas), "ur1")[SatP]) + Q*(AMT(Domelist(ϕ(gas), "ur2")[SatP]) - AMT(Domelist(ϕ(gas), "ur1")[SatP]))
        
        Mol ? urf2 = urf1*M(gas) : urf2 = urf1 
        
        return Pc(gas)*vc(gas)*urf2
        
    else 
        
        urf1 = AMT(C1()) + (8*Tr(Tc(gas), T)*ϕ(gas)/3) - (AMT(3)/vr(vc(gas), v))
        
        Mol ? urf2 = urf1*M(gas) : urf2 = urf1 
        
        return Pc(gas)*vc(gas)*urf2
        
    end
    
end

function h_vdw(gas::vdWGas, T::sysT{Float64,EX}, v::vAmt{Float64,EX,MA}, Mol::Bool = false)
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))
    
    if SatP > 0 && AMT(vr1list[SatP]) < vr(vc(gas), v) < AMT(vr2list[SatP])
        
        vr1 = AMT(vr1list[SatP])
        
        vr2 = AMT(vr2list[SatP])
        
        Q = (vr(vc(gas), v) - vr1)/(vr2 - vr1)
        
        hrf1 = AMT(Domelist(ϕ(gas), "hr1")[SatP]) + Q*(AMT(Domelist(ϕ(gas), "hr2")[SatP]) - AMT(Domelist(ϕ(gas), "hr1")[SatP]))
        
        Mol ? hrf2 = hrf1*M(gas) : hrf2 = hrf1 
        
        return Pc(gas)*vc(gas)*hrf2
        
    else 
        
        hrf1 = AMT(C1()) + (8*Tr(Tc(gas), T)*ϕ(gas)/3) + (8*Tr(Tc(gas), T)*vr(vc(gas), v)/(3*vr(vc(gas), v) - AMT(1))) - (AMT(6)/vr(vc(gas), v))
        
        Mol ? hrf2 = hrf1*M(gas) : hrf2 = hrf1 
        
        return Pc(gas)*vc(gas)*hrf2
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, s::sAmt{Float64,EX,MA})
    
    FQ = FindQ(vr(vc(gas), v), sr(s, Pc(gas), vc(gas), Tc(gas)), vr1list, vr2list, Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))
    
    if FQ == "out"
        
        return Tc(gas)*(exp((3/(8*ϕ(gas)))*(AMTConvert(sr(s,Pc(gas), vc(gas), Tc(gas))) + C2())))*((3*vr(vc(gas), v) - AMT(1))^(-1/ϕ(gas)))
        
    else
        
        vlr = FQ[2]
        
        return Tc(gas)*Tr_sat_list[findclosest(vr1list, AMT(vlr), (10^-3))]
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, u::uAmt{Float64,EX,MA})
    
    FQ = FindQ(vr(vc(gas), v), ur(u, Pc(gas), vc(gas)), vr1list, vr2list, Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))
    
    if FQ == "out"
        
        return Tc(gas)*(3/(8*ϕ(gas)))*(ur(u, Pc(gas), vc(gas)) - AMT(C1()) + (AMT(3)/vr(vc(gas), v)))
        
    else
        
        vlr = FQ[2]
        
        return Tc(gas)*Tr_sat_list[findclosest(vr1list, AMT(vlr), (10^-3))]
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, h::hAmt{Float64,EX,MA})
    
    FQ = FindQ(vr(vc(gas), v), hr(h, Pc(gas), vc(gas)), vr1list, vr2list, Domelist(ϕ(gas), "hr1"), Domelist(ϕ(gas), "hr2")) 
    
    if FQ == "out"
        
        return Tc(gas)*(hr(h, Pc(gas), vc(gas)) - AMT(C1()) + (AMT(6)/vr(vc(gas), v)))/(AMT(8*ϕ(gas)/3) + (8*vr(vc(gas), v)/(3*vr(vc(gas), v) - AMT(1))))
        
    else
        
        vlr = FQ[2]
        
        return Tc(gas)*Tr_sat_list[findclosest(vr1list, AMT(vlr), (10^-3))]
        
    end
    
end

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, u::uAmt{Float64,EX,MA}, Mol::Bool = false)
    
    SatP = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))
        
    if SatP > 0 && AMT(Domelist(ϕ(gas), "ur1")[SatP]) < ur(u, Pc(gas), vc(gas)) < AMT(Domelist(ϕ(gas), "ur2")[SatP])
        
        ur1 = AMT(Domelist(ϕ(gas), "ur1")[SatP])
        
        ur2 = AMT(Domelist(ϕ(gas), "ur2")[SatP])
        
        Q = (ur(u ,Pc(gas), vc(gas)) - ur1)/(ur2 - ur1)
        
        vrf1 = AMT(vr1list[SatP]) + Q*(vr2list[SatP] - vr1list[SatP])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        Proots = AMTConvert(Pr(Pc(gas), P))
        
        uroots = AMTConvert(ur(u, Pc(gas), vc(gas)))
        
        vrf1 = roots(Poly([-3,(9 - (9/ϕ(gas))),(-Proots - (uroots - C1())*(3/ϕ(gas))),3*Proots]))
        
        vrf1 = ImVerification(vrf1)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, h::hAmt{Float64,EX,MA}, Mol::Bool = false)
    
    SatP = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))
                
    if SatP > 0 && AMT(Domelist(ϕ(gas), "hr1")[SatP]) < hr(h, Pc(gas), vc(gas)) < AMT(Domelist(ϕ(gas), "hr2")[SatP])
        
        hr1 = AMT(Domelist(ϕ(gas), "hr1")[SatP])
        
        hr2 = AMT(Domelist(ϕ(gas), "hr2")[SatP])
        
        Q = (hr(h ,Pc(gas), vc(gas)) - hr1)/(hr2 - hr1)
        
        vrf1 = AMT(vr1list[SatP]) + Q*(vr2list[SatP] - vr1list[SatP])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        Proots = AMTConvert(Pr(Pc(gas), P))
        
        hroots = AMTConvert(hr(h, Pc(gas), vc(gas)))
        
        vrf1 = roots(Poly([(-ϕ(gas)),(3*ϕ(gas) - 3),(-(ϕ(gas)*Proots/3) - (hroots) - C1()),(ϕ(gas)*Proots + Proots)]))
        
        vrf1 = ImVerification(vrf1)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end 

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, s::sAmt{Float64,EX,MA}, Mol::Bool = false)
    
    SatP = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))        
        
    if SatP > 0 && AMT(Domelist(ϕ(gas), "sr1")[SatP]) < sr(s, Pc(gas), vc(gas), Tc(gas)) < AMT(Domelist(ϕ(gas), "sr2")[SatP])
        
        sr1 = AMT(Domelist(ϕ(gas), "sr1")[SatP])
        
        sr2 = AMT(Domelist(ϕ(gas), "sr2")[SatP])
        
        Q = (sr(s ,Pc(gas), vc(gas), Tc(gas)) - sr1)/(sr2 - sr1)
        
        vrf1 = AMT(vr1list[SatP]) + Q*(vr2list[SatP] - vr1list[SatP])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        Proots = AMTConvert(Pr(Pc(gas), P))
        
        sroots = AMTConvert(sr(s, Pc(gas), vc(gas), Tc(gas)))        
        
        f(vr) = (((8*exp((3/(8*ϕ(gas)))*(sroots) + C2())))/((3*vr - 1)^(1 + (1/ϕ(gas))))) - (3/vr^2) - Proots
        
        vrf1 = find_zero(f,0.5,Order1())
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end

function v_vdw(gas::vdWGas, T::sysT{Float64,EX}, u::uAmt{Float64,EX,MA}, Mol::Bool = false)
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))
            
    if SatP > 0 && AMT(Domelist(ϕ(gas), "ur1")[SatP]) < ur(u, Pc(gas), vc(gas)) < AMT(Domelist(ϕ(gas), "ur2")[SatP])
        
        ur1 = AMT(Domelist(ϕ(gas), "ur1")[SatP])

        ur2 = AMT(Domelist(ϕ(gas), "ur2")[SatP])
        
        Q = (ur(u ,Pc(gas), vc(gas)) - ur1)/(ur2 - ur1)
        
        vrf1 = AMT(vr1list[SatP]) + Q*(vr2list[SatP] - vr1list[SatP])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        vrf1 = -AMT(3)/(ur(u, Pc(gas), vc(gas)) - AMT(C1()) - (8*Tr(Tc(gas), T)*ϕ(gas)/3))
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end

function v_vdw(gas::vdWGas, T::sysT{Float64,EX}, s::sAmt{Float64,EX,MA}, Mol::Bool = false)
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))        
    
    if SatP > 0 && AMT(Domelist(ϕ(gas), "sr1")[SatP]) < sr(s, Pc(gas), vc(gas), Tc(gas)) < AMT(Domelist(ϕ(gas), "sr2")[SatP])
        
        sr1 = AMT(Domelist(ϕ(gas), "sr1")[SatP])

        sr2 = AMT(Domelist(ϕ(gas), "sr2")[SatP])
        
        Q = (sr(s ,Pc(gas), vc(gas), Tc(gas)) - sr1)/(sr2 - sr1)
        
        vrf1 = AMT(vr1list[SatP]) + Q*(vr2list[SatP] - vr1list[SatP])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        vrf1 = (1/3)*(((Tr(Tc(gas), T)/exp((3/(8*ϕ(gas)))*(AMTConvert(sr(s, Pc(gas), vc(gas), Tc(gas))) + C2())))^(-ϕ(gas))) + AMT(1))
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end

function v_vdw(gas::vdWGas, u::uAmt{Float64,EX,MA}, s::sAmt{Float64,EX,MA}, Mol::Bool = false)
    
    FQ = FindQ(sr(s, Pc(gas), vc(gas), Tc(gas)), ur(u, Pc(gas), vc(gas)), Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"), Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))
    
    if FQ == "out"
        
        uroots = AMTConvert(ur(u, Pc(gas), vc(gas)))
        
        sroots = AMTConvert(sr(s, Pc(gas), vc(gas), Tc(gas)))
        
        f(vr) = (8*ϕ(gas)/3)*(exp((3/(8*ϕ(gas)))*(sroots + C2())))*((3*vr - 1)^(-1/ϕ(gas))) - uroots - (3/vr) + C1()
    
        vrf1 = find_zero(f,0.5,Order1()) #metodo da secante
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2         
        
    else
        
        Q = FQ[1]
        
        slr = AMT(FQ[2])
        
        svr = ((sr(s, Pc(gas), vc(gas), Tc(gas)) - slr)/Q) + slr
        
        vlr = vr1list[findclosest(Domelist(ϕ(gas), "sr1"), slr, (10^-3))]
        
        vvr = vr2list[findclosest(Domelist(ϕ(gas), "sr2"), svr, (10^-3))]
        
        vrf1 = AMT(vlr) + Q*(vvr - vlr)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2    
        
    end
    
end

function v_vdw(gas::vdWGas, u::uAmt{Float64,EX,MA}, h::hAmt{Float64,EX,MA}, Mol::Bool = false)
    
    FQ = FindQ(hr(h, Pc(gas), vc(gas)), ur(u, Pc(gas), vc(gas)), Domelist(ϕ(gas), "hr1"), Domelist(ϕ(gas), "hr2"), Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))
    
    if FQ == "out"
        
        uroots = AMTConvert(ur(u, Pc(gas), vc(gas)))
        
        hroots = AMTConvert(hr(h, Pc(gas), vc(gas)))
        
        vrf1 = roots(Poly([ϕ(gas),((-ϕ(gas)/3)*(uroots) - C1()) + 3*ϕ(gas) + 3 + (ϕ(gas)/3)*(hroots - C1()) - 6*ϕ(gas),((ϕ(gas) + 1)*(uroots - C1()) - ϕ(gas)*(hroots - C1()))]))[2]
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2        
        
    else
        
        Q = FQ[1]
        
        hlr = AMT(FQ[2])
        
        hvr = ((hr(h, Pc(gas), vc(gas)) - hlr)/Q) + hlr
        
        vlr = vr1list[findclosest(Domelist(ϕ(gas), "hr1"), hlr, (10^-3))]
        
        vvr = vr2list[findclosest(Domelist(ϕ(gas), "hr2"), hvr, (10^-3))]
        
        vrf1 = AMT(vlr) + Q*(vvr - vlr)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2   
        
    end
    
end

function v_vdw(gas::vdWGas, s::sAmt{Float64,EX,MA}, h::hAmt{Float64,EX,MA}, Mol::Bool = false)
    
    FQ = FindQ(hr(h, Pc(gas), vc(gas)), sr(s, Pc(gas), vc(gas), Tc(gas)), Domelist(ϕ(gas), "hr1"), Domelist(ϕ(gas), "hr2"), Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))
    
    if FQ == "out"
        
        hroots = AMTConvert(hr(h, Pc(gas), vc(gas)))
        
        sroots = AMTConvert(sr(s, Pc(gas), vc(gas), Tc(gas)))
        
        f(vr) = ((8*ϕ(gas)/3) + (8*vr/(3*vr - 1)))*exp((3/(8*ϕ(gas)))*(sroots + C2()))*((3*vr - 1)^(-1/ϕ(gas))) - hroots + C1() - (6/vr)
    
        vrf1 = find_zero(f,0.5,Order1()) #metodo da secante
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2 
        
    else
        
        Q = FQ[1]
        
        hlr = AMT(FQ[2])
        
        hvr = ((hr(h, Pc(gas), vc(gas)) - hlr)/Q) + hlr
        
        vlr = vr1list[findclosest(Domelist(ϕ(gas), "hr1"), hlr, (10^-3))]
        
        vvr = vr2list[findclosest(Domelist(ϕ(gas), "hr2"), hvr, (10^-3))]
        
        vrf1 = AMT(vlr) + Q*(vvr - vlr)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2      
        
    end
    
end

function v_vdw(gas::vdWGas, T::sysT{Float64,EX}, h::hAmt{Float64,EX,MA}, Mol::Bool = false)
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))
        
    if  SatP > 0 && AMT(Domelist(ϕ(gas), "hr1")[SatP]) < hr(h, Pc(gas), vc(gas)) < AMT(Domelist(ϕ(gas), "hr2")[SatP])
        
        hr1 = AMT(Domelist(ϕ(gas), "hr1")[SatP])

        hr2 = AMT(Domelist(ϕ(gas), "hr2")[SatP])
        
        Q = (hr(h ,Pc(gas), vc(gas)) - hr1)/(hr2 - hr1)
        
        vrf1 = AMT(vr1list[SatP]) + Q*(vr2list[SatP] - vr1list[SatP])
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    else 
        
        Troots = AMTConvert(Tr(Tc(gas), T))
        
        hroots = AMTConvert(hr(h, Pc(gas), vc(gas)))
        
        vrf1 = roots(Poly([6, ((-8*ϕ(gas)*Troots/3) + hroots) - C1() - 18, (8*Troots*(ϕ(gas) + 1) - 3*hroots + 3*C1())]))[2]
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return vc(gas)*vrf2
        
    end
    
end

# functions of properties not given as arguments

ar(gas::vdWGas, a::aAmt{Float64,EX,MA}) = a/(Pc(gas)*vc(gas))

cpr(gas::vdWGas, cp::cpAmt{Float64,EX,MA}) = (Tc(gas)*cp)/(Pc(gas)*vc(gas))

cvr(gas::vdWGas, cv::cvAmt{Float64,EX,MA}) = (Tc(gas)*cv)/(Pc(gas)*vc(gas))

function a_vdw(gas::vdWGas, T::sysT{Float64,EX}, v::vAmt{Float64,EX,MA}, Mol::Bool = false)
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))        
    
    if SatP > 0 && AMT(vr1list[SatP]) < vr(vc(gas), v) < AMT(vr2list[SatP])
        
        vlr = AMT(vr1list[SatP])

        vvr = AMT(vr2list[SatP])
        
        Q = (vr(vc(gas), v) - vlr)/(vvr - vlr)
       
        arf1 = AMT(Domelist(ϕ(gas), "ar1")[SatP]) + Q*(Domelist(ϕ(gas), "ar2")[SatP] - Domelist(ϕ(gas), "ar1")[SatP])
        
        Mol ? arf2 = arf1*M(gas) : arf2 = arf1 
        
        return Pc(gas)*vc(gas)*arf2
        
    else 
        
        arf1 = AMT(C1()) + Tr(Tc(gas), T)*(C2() - (ϕ(gas)*log(AMTConvert(Tr(Tc(gas), T)))/Zc) + (ϕ(gas)/Zc)) - (8/3)*(Tr(Tc(gas), T)*log(3*AMTConvert(vr(vc(gas), v)) - 1)) - (AMT(3)/vr(vc(gas), v))
        
        Mol ? arf2 = arf1*M(gas) : arf2 = arf1 
        
        return Pc(gas)*vc(gas)*arf2
        
    end
    
end

cv_vdw(gas::vdWGas) = (8/3)*ϕ(gas)

function cp_vdw(gas::vdWGas, T::sysT{Float64,EX}, v::vAmt{Float64,EX,MA}, Mol::Bool = false)
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))        
    
    if SatP > 0 && AMT(vr1list[SatP]) < vr(vc(gas), v) < AMT(vr2list[SatP])
        
        vlr = AMT(vr1list[SatP])

        vvr = AMT(vr2list[SatP])
        
        Q = (vr(vc(gas), v) - vlr)/(vvr - vlr)
       
        cprf1 = AMT(Domelist(ϕ(gas), "cpr1")[SatP]) + Q*(Domelist(ϕ(gas), "cpr2")[SatP] - Domelist(ϕ(gas), "cpr1")[SatP])
        
        Mol ? cprf2 = cprf1*M(gas) : cprf2 = cprf1 
        
        return Pc(gas)*vc(gas)*cprf2/Tc(gas)
        
    else 
        
        cprf1 = (8*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) + ϕ(gas)*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - (3*vr(vc(gas), v) - 1)^2)))/(3*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - (3*vr(vc(gas), v) - 1)^2))
        
        Mol ? cprf2 = cprf1*M(gas) : cprf2 = cprf1 
        
        return Pc(gas)*vc(gas)*cprf2/Tc(gas)
        
    end
    
end

γ(cp::cpAmt{Float64,EX,MA}, cv::cvAmt{Float64,EX,MA}) = cp/cv

β(gas::vdWGas, v::vAmt{Float64,EX,MA}, P::sysP{Float64,EX}) = 8*(vr(vc(gas), v)^2)/(3*Tc(gas)*(Pr(Pc(gas), P)*(vr(vc(gas), v)^3) - 3*vr(vc(gas), v) + AMT(2)))

Ks(gas::vdWGas, v::vAmt{Float64,EX,MA}, T::sysT{Float64,EX}) = (vr(vc(gas), v)^2)*((3*vr(vc(gas), v) - AMT(1))^2)/(6*Pc(gas)*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - 9*(vr(vc(gas),v)^2) + 6*vr(vc(gas), v) - AMT(1)))

Kt(gas::vdWGas, v::vAmt{Float64,EX,MA}, T::sysT{Float64,EX}) = Ks(gas, v, T)*ϕ(gas)*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - (3*vr(vc(gas), v) - AMT(1))^2)/(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) + ϕ(gas)*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - (3*vr(vc(gas), v) - AMT(1))^2))

k(gas::vdWGas, v::vAmt{Float64,EX,MA}, T::sysT{Float64,EX}) = 6*(4*Tr(Tc(gas), T)*ϕ(gas)*(vr(vc(gas), v)^3) + 4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - 9*ϕ(gas)*(vr(vc(gas), v)^2) + 6*ϕ(gas)*vr(vc(gas), v) - ϕ(gas))/(ϕ(gas)*(3*vr(vc(gas), v) - AMT(1))*(8*Tr(Tc(gas), T)*(vr(vc(gas), v)^2) - 9*vr(vc(gas), v) + AMT(3)))

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

end #module