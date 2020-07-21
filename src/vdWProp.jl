module vdWProp

using EngThermBase

using Polynomials

using Roots

include("dome.jl")

include("substances.jl")

# function to find the quality of a saturated mixture when P and T are not given, and also find if the pair is inside or outside of the dome.

function FindQ(a::_Amt{Float64,EX}, b::_Amt{Float64,EX}, aArray1::Array, aArray2::Array, bArray1::Array, bArray2::Array)

    i = 1
    
    while i <= points
        
        Q = ((a - AMT(aArray1[i]))/(AMT(aArray2[i] - aArray1[i])))
        
        y = ((a - AMT(aArray1[i]))/(AMT(aArray2[i] - aArray1[i]))) - ((b - AMT(bArray1[i]))/(AMT(bArray2[i] - bArray1[i])))
        
        Q2 = ((a - AMT(aArray1[i + 1]))/(AMT(aArray2[i + 1] - aArray1[i + 1])))
        
        y2 = ((a - AMT(aArray1[i + 1]))/(AMT(aArray2[i + 1] - aArray1[i + 1]))) - ((b - AMT(bArray1[i + 1]))/(AMT(bArray2[i + 1] - bArray1[i + 1])))
        
        if abs(amt(y).val) < 0.001 && AMT(0) < Q < AMT(1) 
            
            return [Q, aArray1[i]] 
            
        elseif (y*y2 < 0) && (AMT(0) < Q < AMT(1) || AMT(0) < Q2 < AMT(1))
            
            if (AMT(0) < Q < AMT(1)) && (AMT(0) < Q2 < AMT(1))
                
                if abs(y) > abs(y2)
                    
                    return [Q2, aArray1[i + 1]]
                    
                else
                    
                    return [Q, aArray1[i]]
                    
                end
                
            else
                
                if AMT(0) < Q < AMT(1)
                    
                    return [Q, aArray1[i]] 
                
                else
                
                    return [Q2, aArray1[i + 1]] 
                    
                end
                
            end
            
        else
            
            yt(j) = ((a - AMT(aArray1[j]))/(AMT(aArray2[j] - aArray1[j]))) - ((b - AMT(bArray1[j]))/(AMT(bArray2[j] - bArray1[j])))
            
            j1 = Integer(round(i + 0.5*points, digits = 0))
            
            j2 = Integer(round(i + 0.3*points, digits = 0))
            
            j3 = Integer(round(i + 0.1*points, digits = 0))
            
            j4 = Integer(round(i + 0.05*points, digits = 0))
            
            j5 = Integer(round(i + 0.01*points, digits = 0))
            
            j6 = Integer(round(i + 0.001*points, digits = 0))
            
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

# Function to find the point in the array that a property gets the given quality

function FindWithQ(pr::Number, Q::Number, Array1::Array, Array2::Array)
    
    if 0 <= Q <= 1
    
        i = 1

        Eq(i) = round(Q, digits = 3) - round((pr - Array1[i])/(Array2[i] - Array1[i]), digits = 3)

        i1 = Integer(round((points/2), digits = 0))

        i2 = Integer(round((points/3), digits = 0))

        i3 = Integer(round((points/4), digits = 0))

        i4 = Integer(round((points/5), digits = 0))

        i5 = Integer(round((points/10), digits = 0))

        i6 = Integer(round((points/100), digits = 0))

        i7 = Integer(round((points/1000), digits = 0))

        while i <= points

            if Eq(i) == 0

                break
                
            elseif abs(Eq(i)) < abs(Eq(i + 1))
                    
                break

            else

                if (i + points/2) < points && Eq(i + i1)*Eq(i) > 0

                    i = (i + i1)

                elseif (i + points/3) < points && Eq(i + i2)*Eq(i) > 0

                    i = (i + i2)  

                elseif (i + points/4) < points && Eq(i + i3)*Eq(i) > 0

                    i = (i + i3)  

                elseif (i + points/5) < points && Eq(i + i4)*Eq(i) > 0

                    i = (i + i4)  

                elseif (i + points/10) < points && Eq(i + i5)*Eq(i) > 0

                    i = (i + i5)  

                elseif (i + points/100) < points && Eq(i + i6)*Eq(i) > 0

                    i = (i + i6)  

                elseif (i + points/1000) < points && Eq(i + i7)*Eq(i) > 0

                    i = (i + i7)  

                else

                    i = i + 1

                end
                
            end

        end

        return i
        
    else
        
        println("Quality must be between 0 and 1")
        
    end
    
end

# Function to find the real root in the array that will be the result for the 3 degree polinomial

function ImVerification(a::Array)
    
    i = 1
    
    while i > 0 && i < 4
        
        imag(a[i]) == 0 ? break : i = i + 1
        
    end
    
    return Float64(a[i])
    
end

# Units

T1 = T(1)

P1 = P(1)

v1 = v(1)

u1 = u(1)

s1 = s(1)

h1 = h(1)

a1 = a(1)

# Constants

C1 = 0

C2 = 0

# Function to find the closest number in an array where the numbers in sequence only increase

function findclosest(array::Array,x::AMOUNTS{Float64,EX},p::Number)
    
    x = amt(x).val

    for i in 1:points
    
        y = x - array[i]
        
        if maximum(array) < x || minimum(array) > x
            
            return -1
    
        elseif abs(y) < p
            
            return i
            
        elseif i > 1 && abs(y) > abs(x - array[i - 1])
            
            return i - 1
            
        end
        
    end    
    
    return points

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
        
        return Tc(gas)*(((Pr(Pc(gas), P)*(3*vr(vc(gas), v) - AMT(1)))/8) + ((AMT(3)*(3*vr(vc(gas), v) - AMT(1)))/(8*(vr(vc(gas), v)^2))))
        
    end
    
end

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, T::sysT{Float64,EX}, Mol::Bool = false)
    
    SatP = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))

    Proots = amt(Pr(Pc(gas), P)).val
    
    Troots = amt(Tr(Tc(gas), T)).val
    
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
        
        return [s(amt((Pc(gas)*vc(gas)*srf2/Tc(gas))/s1).val), Q]
        
    else 
        
        srf1 = (8/3)*log(3*vr(vc(gas), v) - AMT(1)) + (8*ϕ(gas)/3)*log(Tr(Tc(gas), T)) - AMT(C2)
        
        Mol ? srf2 = srf1*M(gas) : srf2 = srf1 
        
        return [s(amt((Pc(gas)*vc(gas)*srf2/Tc(gas))/s1).val), "out"]
        
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
        
        return [u(amt((Pc(gas)*vc(gas)*urf2)/u1).val), Q]
        
    else
            
        urf1 = AMT(C1) + (8*Tr(Tc(gas), T)*ϕ(gas)/3) - (AMT(3)/vr(vc(gas), v))
        
        Mol ? urf2 = urf1*M(gas) : urf2 = urf1 
        
        return [u(amt((Pc(gas)*vc(gas)*urf2)/u1).val), "out"]
        
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
        
        return [h(amt((Pc(gas)*vc(gas)*hrf2)/h1).val), Q]
        
    else 
        
        hrf1 = AMT(C1) + (8*Tr(Tc(gas), T)*ϕ(gas)/3) + (8*Tr(Tc(gas), T)*vr(vc(gas), v)/(3*vr(vc(gas), v) - AMT(1))) - (AMT(6)/vr(vc(gas), v))
        
        Mol ? hrf2 = hrf1*M(gas) : hrf2 = hrf1 
        
        return [h(amt((Pc(gas)*vc(gas)*hrf2)/h1).val), "out"]
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, s::sAmt{Float64,EX,MA})
    
    FQ = FindQ(vr(vc(gas), v), sr(s, Pc(gas), vc(gas), Tc(gas)), vr1list, vr2list, Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))
    
    if FQ == "out"
        
        return Tc(gas)*(exp((3/(8*ϕ(gas)))*(amt(sr(s,Pc(gas), vc(gas), Tc(gas))).val + C2)))*((3*vr(vc(gas), v) - AMT(1))^(-1/ϕ(gas)))
        
    else
        
        vlr = FQ[2]
        
        return Tc(gas)*Tr_sat_list[findclosest(vr1list, AMT(vlr), (10^-3))]
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, u::uAmt{Float64,EX,MA})
    
    FQ = FindQ(vr(vc(gas), v), ur(u, Pc(gas), vc(gas)), vr1list, vr2list, Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))
    
    if FQ == "out"
        
        return Tc(gas)*(3/(8*ϕ(gas)))*(ur(u, Pc(gas), vc(gas)) - AMT(C1) + (AMT(3)/vr(vc(gas), v)))
        
    else
        
        vlr = FQ[2]
        
        return Tc(gas)*Tr_sat_list[findclosest(vr1list, AMT(vlr), (10^-3))]
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, h::hAmt{Float64,EX,MA})
    
    FQ = FindQ(vr(vc(gas), v), hr(h, Pc(gas), vc(gas)), vr1list, vr2list, Domelist(ϕ(gas), "hr1"), Domelist(ϕ(gas), "hr2")) 
    
    if FQ == "out"
        
        return Tc(gas)*(hr(h, Pc(gas), vc(gas)) - AMT(C1) + (AMT(6)/vr(vc(gas), v)))/(AMT(8*ϕ(gas)/3) + (8*vr(vc(gas), v)/(3*vr(vc(gas), v) - AMT(1))))
        
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
        
        return [vc(gas)*vrf2, Q]
        
    else 
        
        Proots = amt(Pr(Pc(gas), P)).val
        
        uroots = amt(ur(u, Pc(gas), vc(gas))).val
        
        vrf1 = roots(Poly([-3,(9 - (9/ϕ(gas))),(-Proots - (uroots - C1)*(3/ϕ(gas))),3*Proots]))
        
        vrf1 = ImVerification(vrf1)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return [vc(gas)*vrf2, "out"]
        
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
        
        return [vc(gas)*vrf2, Q]
        
    else 
        
        Proots = amt(Pr(Pc(gas), P)).val
        
        hroots = amt(hr(h, Pc(gas), vc(gas))).val
        
        vrf1 = roots(Poly([(-ϕ(gas)),(3*ϕ(gas) - 3),(-(ϕ(gas)*Proots/3) - (hroots) - C1),(ϕ(gas)*Proots + Proots)]))
        
        vrf1 = ImVerification(vrf1)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return [vc(gas)*vrf2, "out"]
        
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
        
        return [vc(gas)*vrf2, Q]
        
    else 
        
        Proots = amt(Pr(Pc(gas), P)).val
        
        sroots = amt(sr(s, Pc(gas), vc(gas), Tc(gas))).val        
        
        f(vr) = (((8*exp((3/(8*ϕ(gas)))*(sroots) + C2)))/((3*vr - 1)^(1 + (1/ϕ(gas))))) - (3/vr^2) - Proots
        
        vrf1 = find_zero(f,0.5,Order1())
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return [vc(gas)*vrf2, "out"]
        
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
        
        return [vc(gas)*vrf2, Q]
        
    else 
        
        p1 = amt(ur(u, Pc(gas), vc(gas))).val
        
        p2 = (8*amt(Tr(Tc(gas), T)).val*ϕ(gas)/3)
        
        vrf1 = (-3)/(p1 - C1 - p2)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return [vc(gas)*vrf2, "out"]
        
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
        
        return [vc(gas)*vrf2, Q]
        
    else 
        
        vrf1 = (1/3)*(((Tr(Tc(gas), T)/exp((3/(8*ϕ(gas)))*(amt(sr(s, Pc(gas), vc(gas), Tc(gas))).val + C2)))^(-ϕ(gas))) + AMT(1))
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return [vc(gas)*vrf2, "out"]
        
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
        
        return [vc(gas)*vrf2, Q]
        
    else 
        
        Troots = amt(Tr(Tc(gas), T)).val
        
        hroots = amt(hr(h, Pc(gas), vc(gas))).val
        
        vrf1 = roots(Poly([6, ((-8*ϕ(gas)*Troots/3) + hroots) - C1 - 18, (8*Troots*(ϕ(gas) + 1) - 3*hroots + 3*C1)]))[2]
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return [vc(gas)*vrf2, "out"]
        
    end
    
end

# Pairs with only u/s/h

function v_vdw(gas::vdWGas, u::uAmt{Float64,EX,MA}, s::sAmt{Float64,EX,MA}, Mol::Bool = false)
    
    FQ = FindQ(sr(s, Pc(gas), vc(gas), Tc(gas)), ur(u, Pc(gas), vc(gas)), Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"), Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))
    
    if FQ == "out"
        
        uroots = amt(ur(u, Pc(gas), vc(gas))).val
        
        sroots = amt(sr(s, Pc(gas), vc(gas), Tc(gas))).val
        
        f(vr) = (8*ϕ(gas)/3)*(exp((3/(8*ϕ(gas)))*(sroots + C2)))*((3*vr - 1)^(-1/ϕ(gas))) - uroots - (3/vr) + C1
    
        vrf1 = find_zero(f,0.5,Order1()) #metodo da secante
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return [vc(gas)*vrf2, "out"]         
        
    else
        
        Q = FQ[1]
        
        slr = AMT(FQ[2])
        
        svr = ((sr(s, Pc(gas), vc(gas), Tc(gas)) - slr)/Q) + slr
        
        vlr = vr1list[findclosest(Domelist(ϕ(gas), "sr1"), slr, (10^-3))]
        
        vvr = vr2list[findclosest(Domelist(ϕ(gas), "sr2"), svr, (10^-3))]
        
        vrf1 = AMT(vlr) + Q*(vvr - vlr)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return [vc(gas)*vrf2, Q]    
        
    end
    
end

function v_vdw(gas::vdWGas, u::uAmt{Float64,EX,MA}, h::hAmt{Float64,EX,MA}, Mol::Bool = false)
    
    FQ = FindQ(hr(h, Pc(gas), vc(gas)), ur(u, Pc(gas), vc(gas)), Domelist(ϕ(gas), "hr1"), Domelist(ϕ(gas), "hr2"), Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))
    
    if FQ == "out"
        
        uroots = amt(ur(u, Pc(gas), vc(gas))).val
        
        hroots = amt(hr(h, Pc(gas), vc(gas))).val
        
        vrf0 = roots(Poly([-3, (-uroots + (9/ϕ(gas)) + hroots - 9), (3*uroots + (3*uroots/ϕ(gas)) - (3*C1/ϕ(gas)) - 3*hroots)]))
        
        vrf1 = vrf0[2]
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return [vc(gas)*vrf2, "out"]        
        
    else
        
        Q = FQ[1]
        
        hlr = AMT(FQ[2])
        
        hvr = ((hr(h, Pc(gas), vc(gas)) - hlr)/Q) + hlr
        
        vlr = vr1list[findclosest(Domelist(ϕ(gas), "hr1"), hlr, (10^-3))]
        
        vvr = vr2list[findclosest(Domelist(ϕ(gas), "hr2"), hvr, (10^-3))]
        
        vrf1 = AMT(vlr) + Q*(vvr - vlr)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return [vc(gas)*vrf2, Q]   
        
    end
    
end

function v_vdw(gas::vdWGas, s::sAmt{Float64,EX,MA}, h::hAmt{Float64,EX,MA}, Mol::Bool = false)
    
    FQ = FindQ(sr(s, Pc(gas), vc(gas), Tc(gas)), hr(h, Pc(gas), vc(gas)), Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"), Domelist(ϕ(gas), "hr1"), Domelist(ϕ(gas), "hr2"))
    
    if FQ == "out"
        
        hroots = amt(hr(h, Pc(gas), vc(gas))).val
        
        sroots = amt(sr(s, Pc(gas), vc(gas), Tc(gas))).val
        
        f(vr) = ((8*ϕ(gas)/3) + (8*vr/(3*vr - 1)))*exp((3/(8*ϕ(gas)))*(sroots + C2))*((3*vr - 1)^(-1/ϕ(gas))) - hroots + C1 - (6/vr)
    
        vrf1 = find_zero(f,0.5,Order1()) #metodo da secante
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return [vc(gas)*vrf2, "out"] 
        
    else
        
        Q = FQ[1]
        
        slr = AMT(FQ[2])
        
        svr = ((sr(s, Pc(gas), vc(gas), Tc(gas)) - slr)/Q) + slr
        
        vlr = vr1list[findclosest(Domelist(ϕ(gas), "sr1"), slr, (10^-3))]
        
        vvr = vr2list[findclosest(Domelist(ϕ(gas), "sr2"), svr, (10^-3))]
        
        vrf1 = AMT(vlr) + Q*(vvr - vlr)
        
        Mol ? vrf2 = vrf1*M(gas) : vrf2 = vrf1 
        
        return [vc(gas)*vrf2, Q]      
        
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
        
        return a(amt((Pc(gas)*vc(gas)*arf2)/a1).val)
        
    else 
        
        arf1 = AMT(C1) + Tr(Tc(gas), T)*(C2 - (ϕ(gas)*log(amt(Tr(Tc(gas), T)).val)/Zc) + (ϕ(gas)/Zc)) - (8/3)*(Tr(Tc(gas), T)*log(3*amt(vr(vc(gas), v)).val - 1)) - (AMT(3)/vr(vc(gas), v))
        
        Mol ? arf2 = arf1*M(gas) : arf2 = arf1 
        
        return a(amt((Pc(gas)*vc(gas)*arf2)/a1).val)
        
    end
    
end

cv_vdw(gas::vdWGas) = cv((8/3)*ϕ(gas))

function cp_vdw(gas::vdWGas, T::sysT{Float64,EX}, v::vAmt{Float64,EX,MA}, Mol::Bool = false)
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))        
    
    if SatP > 0 && AMT(vr1list[SatP]) < vr(vc(gas), v) < AMT(vr2list[SatP])
        
        vlr = AMT(vr1list[SatP])

        vvr = AMT(vr2list[SatP])
        
        Q = (vr(vc(gas), v) - vlr)/(vvr - vlr)
       
        cprf1 = AMT(Domelist(ϕ(gas), "cpr1")[SatP]) + Q*(Domelist(ϕ(gas), "cpr2")[SatP] - Domelist(ϕ(gas), "cpr1")[SatP])
        
        Mol ? cprf2 = cprf1*M(gas) : cprf2 = cprf1 
        
        return cp(amt((Pc(gas)*vc(gas)*cprf2/Tc(gas))/s1).val)
        
    else 
        
        cprf1 = (8*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) + ϕ(gas)*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - (3*vr(vc(gas), v) - AMT(1))^2)))/(3*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - (3*vr(vc(gas), v) - AMT(1))^2))
        
        Mol ? cprf2 = cprf1*M(gas) : cprf2 = cprf1 
        
        return cp(amt((Pc(gas)*vc(gas)*cprf2/Tc(gas))/s1).val)
        
    end
    
end

gamma(cp::cpAmt{Float64,EX,MA}, cv::cvAmt{Float64,EX,MA}) = cp/cv

beta(gas::vdWGas, v::vAmt{Float64,EX,MA}, P::sysP{Float64,EX}) = 8*(vr(vc(gas), v)^2)/(3*Tc(gas)*(Pr(Pc(gas), P)*(vr(vc(gas), v)^3) - 3*vr(vc(gas), v) + AMT(2)))

ks(gas::vdWGas, v::vAmt{Float64,EX,MA}, T::sysT{Float64,EX}) = (vr(vc(gas), v)^2)*((3*vr(vc(gas), v) - AMT(1))^2)/(6*Pc(gas)*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - 9*(vr(vc(gas),v)^2) + 6*vr(vc(gas), v) - AMT(1)))

kt(gas::vdWGas, v::vAmt{Float64,EX,MA}, T::sysT{Float64,EX}) = ks(gas, v, T)*ϕ(gas)*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - (3*vr(vc(gas), v) - AMT(1))^2)/(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) + ϕ(gas)*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - (3*vr(vc(gas), v) - AMT(1))^2))

k_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, T::sysT{Float64,EX}) = 6*(4*Tr(Tc(gas), T)*ϕ(gas)*(vr(vc(gas), v)^3) + 4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - 9*ϕ(gas)*(vr(vc(gas), v)^2) + 6*ϕ(gas)*vr(vc(gas), v) - AMT(ϕ(gas)))/(ϕ(gas)*(3*vr(vc(gas), v) - AMT(1))*(8*Tr(Tc(gas), T)*(vr(vc(gas), v)^2) - 9*vr(vc(gas), v) + AMT(3)))

# Function with Quality as Argument

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, Q::_Amt{Float64,EX}, Mol::Bool = false)
    
    if 0 <= amt(Q).val <= 1
   
        SatP = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))

        vrl = vr1list[SatP]

        vrv = vr2list[SatP]

        vr = AMT(vrl) + Q*(vrv - vrl)

        Mol ? vrf = vr*M(gas) : vrf = vr

        return vrf*vc(gas)
        
    else
        
        println("Q must be between 0 and 1")
        
    end
    
end

function v_vdw(gas::vdWGas, T::sysT{Float64,EX}, Q::_Amt{Float64,EX}, Mol::Bool = false)
    
    if 0 <= amt(Q).val <= 1
        
        SatT = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))

        vrl = vr1list[SatT]

        vrv = vr2list[SatT]

        vr = AMT(vrl) + Q*(vrv - vrl)

        Mol ? vrf = vr*M(gas) : vrf = vr

        return vrf*vc(gas)
        
    else
        
        println("Q must be between 0 and 1")
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, Q::_Amt{Float64,EX})
   
    vre = amt(vr(vc(gas), v)).val
    
    Q = amt(Q).val
    
    Point = FindWithQ(vre, Q, vr1list, vr2list)
        
    return Tr_sat_list[Point]*Tc(gas)
    
end
    
function T_vdw(gas::vdWGas, u::uAmt{Float64,EX,MA}, Q::_Amt{Float64,EX})
   
    ure = amt(ur(u, Pc(gas), vc(gas))).val
    
    Q = amt(Q).val
    
    Point = FindWithQ(ure, Q, Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))
        
    return Tr_sat_list[Point]*Tc(gas)
    
end
    
function T_vdw(gas::vdWGas, h::hAmt{Float64,EX,MA}, Q::_Amt{Float64,EX})
   
    hre = amt(hr(h, Pc(gas), vc(gas))).val
    
    Q = amt(Q).val
    
    Point = FindWithQ(hre, Q, Domelist(ϕ(gas), "hr1"), Domelist(ϕ(gas), "hr2"))
        
    return Tr_sat_list[Point]*Tc(gas)
    
end

function T_vdw(gas::vdWGas, s::sAmt{Float64,EX,MA}, Q::_Amt{Float64,EX})
   
    sre = amt(sr(s, Pc(gas), vc(gas), Tc(gas))).val
    
    Q = amt(Q).val
    
    Point = FindWithQ(sre, Q, Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))
        
    return Tr_sat_list[Point]*Tc(gas)
    
end
    
# with all the possible pairs implemented, the next step is implement a function that gets all the six properties when a random pair is given

function State(gas::vdWGas, a::AMOUNTS{Float64,EX}, b::AMOUNTS{Float64,EX}, Mol::Bool = false)
    
    ta = typeof(a)
    
    tb = typeof(b)
    
    if (ta == sysP{Float64,EX} && tb == sysT{Float64,EX}) || (tb == sysP{Float64,EX} && ta == sysT{Float64,EX})
        
        ta == sysP{Float64,EX} ? P = a : P = b
        
        tb == sysT{Float64,EX} ? T = b : T = a
        
        v = v_vdw(gas, P, T)
        
        u = u_vdw(gas, T, v)[1]
        
        h = h_vdw(gas, T, v)[1]
        
        s = s_vdw(gas, T, v)[1]
        
        Q = u_vdw(gas, T, v)[2]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == sysP{Float64,EX} && tb == vAmt{Float64,EX,MA}) || (tb == sysP{Float64,EX} && ta == vAmt{Float64,EX,MA})
        
        ta == sysP{Float64,EX} ? P = a : P = b
        
        tb == vAmt{Float64,EX,MA} ? v = b : v = a
        
        T = T_vdw(gas, P, v)
        
        u = u_vdw(gas, T, v)[1]
        
        h = h_vdw(gas, T, v)[1]
        
        s = s_vdw(gas, T, v)[1]
        
        Q = u_vdw(gas, T, v)[2]

        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == sysP{Float64,EX} && tb == uAmt{Float64,EX,MA}) || (tb == sysP{Float64,EX} && ta == uAmt{Float64,EX,MA})
        
        ta == sysP{Float64,EX} ? P = a : P = b
        
        tb == uAmt{Float64,EX,MA} ? u = b : u = a
        
        v = v_vdw(gas, P, u)[1]
        
        T = T_vdw(gas, P, v)
        
        h = h_vdw(gas, T, v)[1]
        
        s = s_vdw(gas, T, v)[1]
        
        Q = v_vdw(gas, P, u)[2]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == sysP{Float64,EX} && tb == hAmt{Float64,EX,MA}) || (tb == sysP{Float64,EX} && ta == hAmt{Float64,EX,MA})
        
        ta == sysP{Float64,EX} ? P = a : P = b
        
        tb == hAmt{Float64,EX,MA} ? h = b : h = a
        
        v = v_vdw(gas, P, h)[1]
        
        T = T_vdw(gas, P, v)
        
        u = u_vdw(gas, T, v)[1]
        
        s = s_vdw(gas, T, v)[1]
        
        Q = v_vdw(gas, P, h)[2]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == sysP{Float64,EX} && tb == sAmt{Float64,EX,MA}) || (tb == sysP{Float64,EX} && ta == sAmt{Float64,EX,MA})
        
        ta == sysP{Float64,EX} ? P = a : P = b
        
        tb == sAmt{Float64,EX,MA} ? s = b : s = a
        
        v = v_vdw(gas, P, s)[1]
        
        T = T_vdw(gas, P, v)
        
        u = u_vdw(gas, T, v)[1]
        
        h = h_vdw(gas, T, v)[1]
        
        Q = v_vdw(gas, P, s)[2]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == sysT{Float64,EX} && tb == vAmt{Float64,EX,MA}) || (tb == sysT{Float64,EX} && ta == vAmt{Float64,EX,MA})
        
        ta == sysT{Float64,EX} ? T = a : T = b
        
        tb == vAmt{Float64,EX,MA} ? v = b : v = a
        
        P = P_vdw(gas, T, v)
        
        u = u_vdw(gas, T, v)[1]
        
        h = h_vdw(gas, T, v)[1]
        
        s = s_vdw(gas, T, v)[1]
        
        Q = u_vdw(gas, T, v)[2]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == sysT{Float64,EX} && tb == sAmt{Float64,EX,MA}) || (tb == sysT{Float64,EX} && ta == sAmt{Float64,EX,MA})
        
        ta == sysT{Float64,EX} ? T = a : T = b
        
        tb == sAmt{Float64,EX,MA} ? s = b : s = a
        
        v = v_vdw(gas, T, s)[1]
        
        P = P_vdw(gas, T, v)
        
        u = u_vdw(gas, T, v)[1]
        
        h = h_vdw(gas, T, v)[1]
        
        Q = v_vdw(gas, T, s)[2]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == sysT{Float64,EX} && tb == hAmt{Float64,EX,MA}) || (tb == sysT{Float64,EX} && ta == hAmt{Float64,EX,MA})
        
        ta == sysT{Float64,EX} ? T = a : T = b
        
        tb == hAmt{Float64,EX,MA} ? h = b : h = a
        
        v = v_vdw(gas, T, h)[1]
        
        P = P_vdw(gas, T, v)
        
        u = u_vdw(gas, T, v)[1]
        
        s = s_vdw(gas, T, v)[1]
        
        Q = v_vdw(gas, T, h)[2]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == sysT{Float64,EX} && tb == uAmt{Float64,EX,MA}) || (tb == sysT{Float64,EX} && ta == uAmt{Float64,EX,MA})
        
        ta == sysT{Float64,EX} ? T = a : T = b
        
        tb == uAmt{Float64,EX,MA} ? u = b : u = a
        
        v = v_vdw(gas, T, u)[1]
        
        P = P_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)[1]
        
        s = s_vdw(gas, T, v)[1]
        
        Q = v_vdw(gas, T, u)[2]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == vAmt{Float64,EX,MA} && tb == sAmt{Float64,EX,MA}) || (tb == vAmt{Float64,EX,MA} && ta == sAmt{Float64,EX,MA})
        
        ta == vAmt{Float64,EX,MA} ? v = a : v = b
        
        tb == sAmt{Float64,EX,MA} ? s = b : s = a
        
        T = T_vdw(gas, v, s)
                
        P = P_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)[1]
        
        u = u_vdw(gas, T, v)[1]
        
        Q = h_vdw(gas, T, v)[2]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == vAmt{Float64,EX,MA} && tb == hAmt{Float64,EX,MA}) || (tb == vAmt{Float64,EX,MA} && ta == hAmt{Float64,EX,MA})
        
        ta == vAmt{Float64,EX,MA} ? v = a : v = b
        
        tb == hAmt{Float64,EX,MA} ? h = b : h = a
        
        T = T_vdw(gas, v, h)
                
        P = P_vdw(gas, T, v)
        
        s = s_vdw(gas, T, v)[1]
        
        u = u_vdw(gas, T, v)[1]
        
        Q = s_vdw(gas, T, v)[2]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == vAmt{Float64,EX,MA} && tb == uAmt{Float64,EX,MA}) || (tb == vAmt{Float64,EX,MA} && ta == uAmt{Float64,EX,MA})
        
        ta == vAmt{Float64,EX,MA} ? v = a : v = b
        
        tb == uAmt{Float64,EX,MA} ? u = b : u = a
        
        T = T_vdw(gas, v, u)
                
        P = P_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)[1]
        
        s = s_vdw(gas, T, v)[1]
        
        Q = h_vdw(gas, T, v)[2]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == sAmt{Float64,EX,MA} && tb == hAmt{Float64,EX,MA}) || (tb == sAmt{Float64,EX,MA} && ta == hAmt{Float64,EX,MA})
        
        ta == sAmt{Float64,EX,MA} ? s = a : s = b
        
        tb == hAmt{Float64,EX,MA} ? h = b : h = a
        
        v = v_vdw(gas, s, h)[1]
        
        T = T_vdw(gas, v, s)
                
        P = P_vdw(gas, T, v)
        
        u = u_vdw(gas, T, v)[1]
        
        Q = v_vdw(gas, s, h)[2]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == sAmt{Float64,EX,MA} && tb == uAmt{Float64,EX,MA}) || (tb == sAmt{Float64,EX,MA} && ta == uAmt{Float64,EX,MA})
        
        ta == sAmt{Float64,EX,MA} ? s = a : s = b
        
        tb == uAmt{Float64,EX,MA} ? u = b : u = a
        
        v = v_vdw(gas, u, s)[1]
        
        T = T_vdw(gas, v, s)
                
        P = P_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)[1]
        
        Q = v_vdw(gas, u, s)[2]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == hAmt{Float64,EX,MA} && tb == uAmt{Float64,EX,MA}) || (tb == hAmt{Float64,EX,MA} && ta == uAmt{Float64,EX,MA})
        
        ta == hAmt{Float64,EX,MA} ? h = a : h = b
        
        tb == uAmt{Float64,EX,MA} ? u = b : u = a
        
        v = v_vdw(gas, u, h)[1]
        
        T = T_vdw(gas, v, h)
                
        P = P_vdw(gas, T, v)
        
        s = s_vdw(gas, T, v)[1]
        
        Q = v_vdw(gas, u, h)[2]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == sysT{Float64,EX} && tb == _Amt{Float64,EX}) || (tb == sysT{Float64,EX} && ta == _Amt{Float64,EX})
        
        ta == sysT{Float64,EX} ? T = a : T = b
        
        tb == _Amt{Float64,EX} ? Q = b : Q = a
        
        v = v_vdw(gas, T, Q)
                
        P = P_vdw(gas, T, v)
        
        u = u_vdw(gas, T, v)[1]
        
        h = h_vdw(gas, T, v)[1]
        
        s = s_vdw(gas, T, v)[1]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == sysP{Float64,EX} && tb == _Amt{Float64,EX}) || (tb == sysP{Float64,EX} && ta == _Amt{Float64,EX})
        
        ta == sysP{Float64,EX} ? P = a : P = b
        
        tb == _Amt{Float64,EX} ? Q = b : Q = a
        
        v = v_vdw(gas, P, Q)
                
        T = T_vdw(gas, P, v)
        
        u = u_vdw(gas, T, v)[1]
        
        h = h_vdw(gas, T, v)[1]
        
        s = s_vdw(gas, T, v)[1]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == vAmt{Float64,EX,MA} && tb == _Amt{Float64,EX}) || (tb == vAmt{Float64,EX,MA} && ta == _Amt{Float64,EX})
        
        ta == vAmt{Float64,EX,MA} ? v = a : v = b
        
        tb == _Amt{Float64,EX} ? Q = b : Q = a
        
        T = T_vdw(gas, v, Q)
                
        P = P_vdw(gas, T, v)
        
        u = u_vdw(gas, T, v)[1]
        
        h = h_vdw(gas, T, v)[1]
        
        s = s_vdw(gas, T, v)[1]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == uAmt{Float64,EX,MA} && tb == _Amt{Float64,EX}) || (tb == uAmt{Float64,EX,MA} && ta == _Amt{Float64,EX})
        
        ta == uAmt{Float64,EX,MA} ? u = a : u = b
        
        tb == _Amt{Float64,EX} ? Q = b : Q = a
        
        T = T_vdw(gas, u, Q)
        
        v = v_vdw(gas, T, u)[1]
                
        P = P_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)[1]
        
        s = s_vdw(gas, T, v)[1]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == hAmt{Float64,EX,MA} && tb == _Amt{Float64,EX}) || (tb == hAmt{Float64,EX,MA} && ta == _Amt{Float64,EX})
        
        ta == hAmt{Float64,EX,MA} ? h = a : h = b
        
        tb == _Amt{Float64,EX} ? Q = b : Q = a
        
        T = T_vdw(gas, h, Q)
        
        v = v_vdw(gas, T, h)[1]
                
        P = P_vdw(gas, T, v)
        
        u = u_vdw(gas, T, v)[1]
        
        s = s_vdw(gas, T, v)[1]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    elseif (ta == sAmt{Float64,EX,MA} && tb == _Amt{Float64,EX}) || (tb == sAmt{Float64,EX,MA} && ta == _Amt{Float64,EX})
        
        ta == sAmt{Float64,EX,MA} ? s = a : s = b
        
        tb == _Amt{Float64,EX} ? Q = b : Q = a
        
        T = T_vdw(gas, s, Q)
        
        v = v_vdw(gas, T, s)[1]
                
        P = P_vdw(gas, T, v)
        
        h = h_vdw(gas, T, v)[1]
        
        u = u_vdw(gas, T, v)[1]
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P)
        
        Ks = ks(gas, v, T)
        
        Kt = kt(gas, v, T)
        
        k = k_vdw(gas, v, T)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, Q]
        
        return St
        
    else
        
        println("ERROR, the arguments needs to be properties between P,T,v,h,u,s,Q and the base for the intensive ones needs to be mass.")
        
    end
    
end

# Function to calculate the initial and final State of a process where there is a property that doesn't change

function IsoProp(gas::vdWGas, a1::AMOUNTS{Float64,EX}, b::AMOUNTS{Float64,EX}, a2::AMOUNTS{Float64,EX},iso::String, Mol::Bool = false)
    
    if Mol == false

    St1 = State(gas, a1, b, Mol)
    
    iso == "P" ? St2 = State(gas, a2, St1[1], Mol) : 
    
    iso == "T" ? St2 = State(gas, a2, St1[2], Mol) : 
    
    iso == "v" ? St2 = State(gas, a2, St1[3], Mol) :
    
    iso == "u" ? St2 = State(gas, a2, St1[4], Mol) :
    
    iso == "h" ? St2 = State(gas, a2, St1[5], Mol) :
    
    iso == "s" ? St2 = State(gas, a2, St1[6], Mol) : println("The supported iso properties are P,T,v,u,h,s.")
    
    return hcat(St1, St2)
        
    else
        
    St1 = State(gas, a1, b, Mol)
    
    iso == "P" ? St2 = State(gas, a2, St1[1], Mol) : 
    
    iso == "T" ? St2 = State(gas, a2, St1[2], Mol) : 
    
    iso == "v" ? St2 = State(gas, a2, St1[3]/M(gas), Mol) :
    
    iso == "u" ? St2 = State(gas, a2, St1[4]/M(gas), Mol) :
    
    iso == "h" ? St2 = State(gas, a2, St1[5]/M(gas), Mol) :
    
    iso == "s" ? St2 = State(gas, a2, St1[6]/M(gas), Mol) : println("The supported iso properties are P,T,v,u,h,s.")
    
    return hcat(St1, St2)  
        
    end
    
end

# Functions for tests

function isoh(gas::vdWGas, T1::sysT{Float64,EX}, v1::vAmt{Float64,EX,MA}, T2::sysT{Float64,EX}, v2::vAmt{Float64,EX,MA})
    
    Tr1 = amt(Tr(Tc(gas), T1)).val
    
    Tr2 = amt(Tr(Tc(gas), T2)).val
    
    vr1 = amt(vr(vc(gas), v1)).val
    
    vr2 = amt(vr(vc(gas), v2)).val    
    
    test1 = ((8*ϕ(gas)/3)*(Tr2 - Tr1))
    
    test2 = 8*((Tr2*vr2/(3*vr2 - 1)) - (Tr1*vr1/(3*vr1 - 1)))
    
    test3 = 6*((1/vr2) - (1/vr1))
    
    test4 = test1 + test2 - test3
    
    return test4

end

function isou(gas::vdWGas, T1::sysT{Float64,EX}, v1::vAmt{Float64,EX,MA}, T2::sysT{Float64,EX}, v2::vAmt{Float64,EX,MA})
    
    Tr1 = amt(Tr(Tc(gas), T1)).val
    
    Tr2 = amt(Tr(Tc(gas), T2)).val
    
    vr1 = amt(vr(vc(gas), v1)).val
    
    vr2 = amt(vr(vc(gas), v2)).val 
        
    test1 = (8*ϕ(gas)/9)*(Tr2 - Tr1)
        
    test2 = ((1/vr2) - (1/vr1))
        
    test3 = test1 - test2
    
    return test3

end
    
function isos(gas::vdWGas, T1::sysT{Float64,EX}, v1::vAmt{Float64,EX,MA}, T2::sysT{Float64,EX}, v2::vAmt{Float64,EX,MA})
    
    Tr1 = amt(Tr(Tc(gas), T1)).val
    
    Tr2 = amt(Tr(Tc(gas), T2)).val
    
    vr1 = amt(vr(vc(gas), v1)).val
    
    vr2 = amt(vr(vc(gas), v2)).val 
        
    test1 = (Tr2/Tr1)^(ϕ(gas))
        
    test2 = (3*vr1 - 1)/(3*vr2 - 1)
        
    test3 = test1 - test2
    
    return test3

end
    
function isoP(gas::vdWGas, T1::sysT{Float64,EX}, v1::vAmt{Float64,EX,MA}, T2::sysT{Float64,EX}, v2::vAmt{Float64,EX,MA})
    
    Tr1 = amt(Tr(Tc(gas), T1)).val
    
    Tr2 = amt(Tr(Tc(gas), T2)).val
    
    vr1 = amt(vr(vc(gas), v1)).val
    
    vr2 = amt(vr(vc(gas), v2)).val 
        
    test1 = (Tr2/(3*vr2 - 1)) - (Tr1/(3*vr1 - 1))
        
    test2 = (3/8)*((1/(vr2^2)) - (1/(vr1^2)))
        
    test3 = test1 - test2
    
    return test3

end
    
function isoT(gas::vdWGas, P1::sysP{Float64,EX}, v1::vAmt{Float64,EX,MA}, P2::sysP{Float64,EX}, v2::vAmt{Float64,EX,MA})
    
    Pr1 = amt(Pr(Pc(gas), P1)).val
    
    Pr2 = amt(Pr(Pc(gas), P2)).val
    
    vr1 = amt(vr(vc(gas), v1)).val
    
    vr2 = amt(vr(vc(gas), v2)).val 
        
    test1 = (Pr2 - Pr1)
        
    test2 = 3*((Pr1*vr1 - Pr2*vr2) + (((3*vr1 - 1)/(vr1^2)) - ((3*vr2 - 1)/(vr2^2))))
        
    test3 = test1 + test2
    
    return test3

end

export State

export IsoProp

end #module
