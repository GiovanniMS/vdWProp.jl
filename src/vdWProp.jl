module vdWProp

using EngThermBase

using Polynomials

using Polynomials.PolyCompat

using Roots

include("dome.jl")

include("substances.jl")

# function to find the quality of a saturated mixture when P and T are not given, and also find if the pair is inside or outside of the dome.

#function FindQ(a::_Amt{Float64,EX}, b::_Amt{Float64,EX}, aArray1::Array, aArray2::Array, bArray1::Array, bArray2::Array)

#    i = 1
    
#    while i <= points
        
        #Q = ((a - AMT(aArray1[i]))/(AMT(aArray2[i] - aArray1[i])))
        
#        y = ((a - AMT(aArray1[i]))/(AMT(aArray2[i] - aArray1[i]))) - ((b - AMT(bArray1[i]))/(AMT(bArray2[i] - bArray1[i])))
        
        #if abs(amt(y).val) < 0.001 && AMT(0) < Q < AMT(1) 
            
        #    return [Q, aArray1[i]] 
            
        #else
            
#            yt(j) = ((a - AMT(aArray1[j]))/(AMT(aArray2[j] - aArray1[j]))) - ((b - AMT(bArray1[j]))/(AMT(bArray2[j] - bArray1[j])))
            
#            j1 = Integer(round(i + 0.5*points, digits = 0))
            
#            j2 = Integer(round(i + 0.3*points, digits = 0))
            
#            j3 = Integer(round(i + 0.1*points, digits = 0))
            
#            j4 = Integer(round(i + 0.05*points, digits = 0))
            
#            j5 = Integer(round(i + 0.01*points, digits = 0))
            
#            j6 = Integer(round(i + 0.001*points, digits = 0))
            
#            if j1 <= points && yt(j1)*y > AMT(0)
                    
#                i = j1
                
#            elseif j2 <= points && yt(j2)*y > AMT(0)
                    
#                i = j2   
                    
#            elseif j3 <= points && yt(j3)*y > AMT(0)
                    
#                i = j3
                        
#            elseif j4 <= points && yt(j4)*y > AMT(0)
                    
#                i = j4
                            
#            elseif j5 <= points && yt(j5)*y > AMT(0)
                    
#                i = j5
                
#            elseif j6 <= points && yt(j6)*y > AMT(0)
                    
#                i = j6
                
#            elseif (i + 1) <= points && yt(i + 1)*y < AMT(0)
                    
#                if abs(amt(y).val) < abs(amt(yt(i + 1)).val)
                
#                    Q = ((a - AMT(aArray1[i]))/(AMT(aArray2[i] - aArray1[i])))
                
#                    println(i)
                
#                    return [Q, aArray1[i]] 
                    
#                else
                    
#                    Q = ((a - AMT(aArray1[i + 1]))/(AMT(aArray2[i + 1] - aArray1[i + 1])))
                
#                    println(i)
                    
#                    return [Q, aArray1[i + 1]] 
                    
#                end
            
#            else
                
#                i = i + 1
                
#            end
            
        #end
        
#    end
    
#    return "out"
    
#end

function FindQ(a::_Amt{Float64,EX}, b::_Amt{Float64,EX}, aArray1::Array, aArray2::Array, bArray1::Array, bArray2::Array)

    i = 1
    
    Qarray = []
    
    yarray = []
    
    while i <= points
        
        Q = ((a - AMT(aArray1[i]))/(AMT(aArray2[i] - aArray1[i])))
        
        y = ((a - AMT(aArray1[i]))/(AMT(aArray2[i] - aArray1[i]))) - ((b - AMT(bArray1[i]))/(AMT(bArray2[i] - bArray1[i])))
        
        append!(Qarray, amt(Q).val)
        
        append!(yarray, abs(amt(y).val))
        
        i = i + 1
        
    end
    
    min = minimum(yarray)
    
    ic = findall(yarray .== min)[1]
    
    if min < (10^-2) && 0 <= round(Qarray[ic], digits = 3) <= 1

        return [AMT(Qarray[ic]), bArray1[ic]]
        
    else
    
        return "out"
        
    end
    
end

# Function to find the point in the array that a property gets the given quality

function FindWithQ(pr::Number, Q::Number, Array1::Array, Array2::Array)
    
    if 0 <= Q <= 1
        
        Eqarray1 = []
        
        #Eqarray2 = []
    
        i = 1

        Eq(i) = Q - (pr - Array1[i])/(Array2[i] - Array1[i])

        # rev 4 while i <= round(0.2647*points, digits = 0)
        
        #while i <= round(0.1872*points, digits = 0)

        #    te  = abs(Eq(i))
            
        #    append!(Eqarray1, te)
            
        #    i = i + 1

        #end
        
        while i <= points

            te  = abs(Eq(i))
            
            append!(Eqarray1, te)
            
            i = i + 1

        end
        
        min1 = minimum(Eqarray1)
        
        #min2 = minimum(Eqarray2)
    
        if min1 < (10^-2)

            ic = findall(Eqarray1 .== min1)[1]

            return ic
            
        #elseif min2 < (10^-2)

        #    ic = findall(Eqarray2 .== min2)[1] + length(Eqarray1)

        #    return ic

        else
            
            println("State not supported")
            
        end
        
    else
        
        println("Quality must be between 0 and 1")
        
    end
    
end

#function FindWithQ(pr::Number, Q::Number, Array1::Array, Array2::Array)
    
#    if 0 <= Q <= 1
        
#        Eqarray1 = []
        
#        Eqarray2 = []
        
#        Eqarray3 = []
        
#        Eqarray4 = []
    
#        i = 1

#        Eq(i) = Q - (pr - Array1[i])/(Array2[i] - Array1[i])

#        while i <= round(points*0.7, digits = 0)

#            te  = abs(Eq(i))
            
#            append!(Eqarray1, te)
            
#            i = i + 1

#        end
        
#        while i <= round(points*0.85, digits = 0)

#            te  = abs(Eq(i))
            
#            append!(Eqarray2, te)
            
#            i = i + 1

#        end
        
#        while i <= round(points*0.97, digits = 0)

#            te  = abs(Eq(i))
            
#            append!(Eqarray3, te)
            
#            i = i + 1

#        end
        
#        while i <= round(points, digits = 0)

#            te  = abs(Eq(i))
            
#            append!(Eqarray4, te)
            
#            i = i + 1

#        end
        
#        min4 = minimum(Eqarray4)
        
#        min3 = minimum(Eqarray3)
        
#        min2 = minimum(Eqarray2)
        
#        min1 = minimum(Eqarray1)
        
#        println(findall(Eqarray1 .== min1)[1])
        
#        println(findall(Eqarray2 .== min2)[1])
        
#        println(findall(Eqarray3 .== min3)[1])
        
#        println(findall(Eqarray4 .== min4)[1])
    
#        if min1 < (10^-2)

#            ic = findall(Eqarray1 .== min1)[1]

#            return ic
            
#        elseif min2 < (10^-7)

#            ic = findall(Eqarray2 .== min2)[1] + length(Eqarray1)

#            return ic
                
#        elseif min3 < (10^-7)

#            ic = findall(Eqarray3 .== min3)[1] + length(Eqarray1) + length(Eqarray2)

#            return ic
            
#        elseif min4 < (10^-7)

#            ic = findall(Eqarray4 .== min4)[1] + length(Eqarray1) + length(Eqarray2) + length(Eqarray3)

#            return ic

#        else
            
#            println("State not supported")
            
#        end
        
#    else
        
#        println("Quality must be between 0 and 1")
        
#    end
    
#end

#function FindWithQ(pr::Number, Q::Number, Array1::Array, Array2::Array)
    
#    if 0 <= Q <= 1
    
#        i = 1

#        Eq(i) = Q - (pr - Array1[i])/(Array2[i] - Array1[i])

#        i1 = Integer(round((points/2), digits = 0))

#        i2 = Integer(round((points/3), digits = 0))

#        i3 = Integer(round((points/4), digits = 0))

#        i4 = Integer(round((points/5), digits = 0))

#        i5 = Integer(round((points/10), digits = 0))

#        i6 = Integer(round((points/100), digits = 0))

#        i7 = Integer(round((points/1000), digits = 0))

#        while i <= points

#            if Eq(i) == 0

#                break
                
#            elseif i < points && abs(Eq(i)) < abs(Eq(i + 1))
                    
#                break

#            else

#                if (i + round((points/2), digits = 0)) < points && Eq(i + i1)*Eq(i) > 0

#                    i = (i + i1)

#                elseif (i + round((points/3), digits = 0)) < points && Eq(i + i2)*Eq(i) > 0

#                    i = (i + i2)  

#                elseif (i + round((points/4), digits = 0)) < points && Eq(i + i3)*Eq(i) > 0

#                    i = (i + i3)  

#                elseif (i + round((points/5), digits = 0)) < points && Eq(i + i4)*Eq(i) > 0

#                    i = (i + i4)  

#                elseif (i + round((points/10), digits = 0)) < points && Eq(i + i5)*Eq(i) > 0

#                    i = (i + i5)  

#                elseif (i + round((points/100), digits = 0)) < points && Eq(i + i6)*Eq(i) > 0

#                    i = (i + i6)  

#                elseif (i + round((points/1000), digits = 0)) < points && Eq(i + i7)*Eq(i) > 0

#                    i = (i + i7)  

#                else

#                    i = i + 1

#                end
                
#            end

#        end
        
#        if i <= points

#            return i
            
#        elseif i == (points + 1) && abs(Eq(points)) < 0.001
            
#            return (i - 1)
            
#        else
            
#            println("State not supported")
            
#        end
        
#    else
        
#        println("Quality must be between 0 and 1")
        
#    end
    
#end

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
    
    xi = Integer(round(x, digits = 0))
    
    yarray = []

    for i in 1:points
    
        y = x - array[i]
        
        append!(yarray, abs(y))
        
    end
    
    min = minimum(yarray)
        
    if (min < 10^(length("$xi"))*10^(-3)) || (x > 10 && min < 10^(length("$xi"))*10^(-2))

        ic = findall(yarray .== min)[1]

        return ic

    else
            
        return -1
            
    end

end

#function findclosest(array::Array,x::AMOUNTS{Float64,EX},p::Number)
    
#    x = amt(x).val

#    for i in 1:points
    
#        y = x - array[i]
        
#        if i > 1 && abs(y) > abs(x - array[i - 1])
            
#            return i - 1
        
#        #if abs(y) < p
            
#            #return i
        
#        elseif maximum(array) < x || minimum(array) > x
            
#            return -1
            
#        end
        
#    end    
    
#    return points

#end

#function findclosest(array::Array,x::AMOUNTS{Float64,EX},p::Number)
    
#    yarray = []
    
#    x = amt(x).val
    
#    for i in 1:points
        
#        y = x - array[i]
        
#        append!(yarray, abs(y))
        
#    end
    
#    min = minimum(yarray)
    
#    if min < (10^-3)
    
#        ic = findall(yarray .== min)[1]

#        return ic
        
#    else
    
#        return "out"
        
#    end

#end

# Now the functions are implemented using the vdWGas as an argument

Tr(Tc::sysT{Float64,EX}, T::sysT{Float64,EX}) = T/Tc

Pr(Pc::sysP{Float64,EX}, P::sysP{Float64,EX}) = P/Pc

vr(vc::vAmt{Float64,EX,MA}, v::vAmt{Float64,EX,MA}) = v/vc

ur(u::uAmt{Float64,EX,MA}, Pc::sysP{Float64,EX}, vc::vAmt{Float64,EX,MA}) = u/(Pc*vc)

hr(h::hAmt{Float64,EX,MA}, Pc::sysP{Float64,EX}, vc::vAmt{Float64,EX,MA}) = h/(Pc*vc)

sr(s::sAmt{Float64,EX,MA}, Pc::sysP{Float64,EX}, vc::vAmt{Float64,EX,MA}, Tc::sysT{Float64,EX}) = (Tc*s)/(Pc*vc)

function P_vdw(gas::vdWGas, T::sysT{Float64,EX}, v::vAmt{Float64,EX,MA}) 
    
    SatP = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3)) # Array position for the saturated fluid and gas 
    
    if SatP > 0 && ((SatP < points && AMT(vr1list[SatP + 1]) < vr(vc(gas), v) < AMT(vr2list[SatP])) || AMT(vr1list[SatP]) < vr(vc(gas), v) < AMT(vr2list[SatP]))# || amt(AMT(vr1list[SatP])).val - amt(vr(vc(gas), v)).val < 10^-5)
        
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
    
    vap = amt(vr(vc(gas), v)).val
    
    if SatP > 0 && (AMT(vr1list[SatP]) < vr(vc(gas), v) < AMT(vr2list[SatP]))# || vap <= round(vr2list[SatP], sigdigits = 3))
        
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
        
        slr = FQ[2]
        
        return Tc(gas)*Tr_sat_list[findclosest(Domelist(ϕ(gas), "sr1"), AMT(slr), (10^-3))]
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, u::uAmt{Float64,EX,MA})
    
    FQ = FindQ(vr(vc(gas), v), ur(u, Pc(gas), vc(gas)), vr1list, vr2list, Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))
    
    if FQ == "out"
        
        return Tc(gas)*(3/(8*ϕ(gas)))*(ur(u, Pc(gas), vc(gas)) - AMT(C1) + (AMT(3)/vr(vc(gas), v)))
        
    else
        
        ulr = FQ[2]
        
        return Tc(gas)*Tr_sat_list[findclosest(Domelist(ϕ(gas), "ur1"), AMT(ulr), (10^-3))]
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, h::hAmt{Float64,EX,MA})
    
    FQ = FindQ(vr(vc(gas), v), hr(h, Pc(gas), vc(gas)), vr1list, vr2list, Domelist(ϕ(gas), "hr1"), Domelist(ϕ(gas), "hr2")) 
    
    if FQ == "out"
        
        return Tc(gas)*(hr(h, Pc(gas), vc(gas)) - AMT(C1) + (AMT(6)/vr(vc(gas), v)))/(AMT(8*ϕ(gas)/3) + (8*vr(vc(gas), v)/(3*vr(vc(gas), v) - AMT(1))))
        
    else
        
        hlr = FQ[2]
        
        return Tc(gas)*Tr_sat_list[findclosest(Domelist(ϕ(gas), "hr1"), AMT(hlr), (10^-3))]
        
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
    
    uap = amt(ur(u, Pc(gas), vc(gas))).val
            
    if SatP > 0 && AMT(Domelist(ϕ(gas), "ur1")[SatP]) < (ur(u, Pc(gas), vc(gas))) < AMT(Domelist(ϕ(gas), "ur2")[SatP])# || uap <= 1.001*(Domelist(ϕ(gas), "ur2")[SatP]))
        
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
    
    if SatP > 0 && AMT(Domelist(ϕ(gas), "sr1")[SatP]) < (sr(s, Pc(gas), vc(gas), Tc(gas))) < AMT(Domelist(ϕ(gas), "sr2")[SatP])

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
    
    #FQ = FindQ(sr(s, Pc(gas), vc(gas), Tc(gas)), ur(u, Pc(gas), vc(gas)), Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"), Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))
    
    FQ = FindQ(ur(u, Pc(gas), vc(gas)), sr(s, Pc(gas), vc(gas), Tc(gas)),Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"), Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))
    
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
    
    #FQ = FindQ(hr(h, Pc(gas), vc(gas)), ur(u, Pc(gas), vc(gas)), Domelist(ϕ(gas), "hr1"), Domelist(ϕ(gas), "hr2"), Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))
    
    FQ = FindQ(ur(u, Pc(gas), vc(gas)), hr(h, Pc(gas), vc(gas)), Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"), Domelist(ϕ(gas), "hr1"), Domelist(ϕ(gas), "hr2"))
    
    if FQ == "out"
        
        uroots = amt(ur(u, Pc(gas), vc(gas))).val
        
        hroots = amt(hr(h, Pc(gas), vc(gas))).val
        
        vrf0 = roots(Poly([ϕ(gas), ((-ϕ(gas)/3)*(uroots - C1) + 3*ϕ(gas) + 3 + (ϕ(gas)/3)*(hroots - C1) - 6*ϕ(gas)), ((ϕ(gas) + 1)*(uroots - C1) - ϕ(gas)*(hroots - C1))]))
        
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
    
    #FQ = FindQ(sr(s, Pc(gas), vc(gas), Tc(gas)), hr(h, Pc(gas), vc(gas)), Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"), Domelist(ϕ(gas), "hr1"), Domelist(ϕ(gas), "hr2"))
    
    FQ = FindQ(hr(h, Pc(gas), vc(gas)), sr(s, Pc(gas), vc(gas), Tc(gas)), Domelist(ϕ(gas), "hr1"), Domelist(ϕ(gas), "hr2"), Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))
    
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
    
    vap = amt(vr(vc(gas), v)).val
    
    if SatP > 0 && (AMT(vr1list[SatP]) < vr(vc(gas), v) < AMT(vr2list[SatP]))# || vap <= round(vr2list[SatP], sigdigits = 3))
        
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

cv_vdw(gas::vdWGas) = cv(amt((Pc(gas)*vc(gas)*((8/3)*ϕ(gas))/Tc(gas))/s1).val)

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

beta1(gas::vdWGas, v::vAmt{Float64,EX,MA}, P::sysP{Float64,EX}) = 8*(vr(vc(gas), v)^2)/(3*Tc(gas)*(Pr(Pc(gas), P)*(vr(vc(gas), v)^3) - 3*vr(vc(gas), v) + AMT(2)))

function beta(gas::vdWGas, v::vAmt{Float64,EX,MA}, P::sysP{Float64,EX}, Q)  
    
    if Q == "out"

        return beta1(gas, v, P)
        
    else
        
        ic = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))
        
        vl = vr1list[ic]*vc(gas)
        
        vv = vr2list[ic]*vc(gas)
        
        return beta1(gas, vl, P) + Q*(beta1(gas, vv, P) - beta1(gas, vl, P))
        
    end
    
end

ks1(gas::vdWGas, v::vAmt{Float64,EX,MA}, T::sysT{Float64,EX}) = (vr(vc(gas), v)^2)*((3*vr(vc(gas), v) - AMT(1))^2)/(6*Pc(gas)*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - 9*(vr(vc(gas),v)^2) + 6*vr(vc(gas), v) - AMT(1)))
    
function ks(gas::vdWGas, v::vAmt{Float64,EX,MA}, T::sysT{Float64,EX}, Q)
    
    if Q == "out"

        return ks1(gas, v, T)
        
    else
        
        ic = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))
        
        vl = vr1list[ic]*vc(gas)
        
        vv = vr2list[ic]*vc(gas)
        
        return ks1(gas, vl, T) + Q*(ks1(gas, vv, T) - ks1(gas, vl, T))
        
    end
        
end

kt1(gas::vdWGas, v::vAmt{Float64,EX,MA}, T::sysT{Float64,EX}, Q) = ks(gas, v, T, Q)*ϕ(gas)*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - (3*vr(vc(gas), v) - AMT(1))^2)/(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) + ϕ(gas)*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - (3*vr(vc(gas), v) - AMT(1))^2))

function kt(gas::vdWGas, v::vAmt{Float64,EX,MA}, T::sysT{Float64,EX}, Q)
    
    if Q == "out"

        return kt1(gas, v, T, Q)
        
    else
        
        ic = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))
        
        vl = vr1list[ic]*vc(gas)
        
        vv = vr2list[ic]*vc(gas)
        
        return kt1(gas, vl, T, Q) + Q*(kt1(gas, vv, T, Q) - kt1(gas, vl, T, Q))
        
    end
    
end

k1_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, T::sysT{Float64,EX}) = 6*(4*Tr(Tc(gas), T)*ϕ(gas)*(vr(vc(gas), v)^3) + 4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - 9*ϕ(gas)*(vr(vc(gas), v)^2) + 6*ϕ(gas)*vr(vc(gas), v) - AMT(ϕ(gas)))/(ϕ(gas)*(3*vr(vc(gas), v) - AMT(1))*(8*Tr(Tc(gas), T)*(vr(vc(gas), v)^2) - 9*vr(vc(gas), v) + AMT(3)))

function k_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, T::sysT{Float64,EX}, Q)

    if Q == "out"

        return k1_vdw(gas, v, T)
        
    else
        
        ic = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))
        
        vl = vr1list[ic]*vc(gas)
        
        vv = vr2list[ic]*vc(gas)
        
        return k1_vdw(gas, vl, T) + Q*(k1_vdw(gas, vv, T) - k1_vdw(gas, vl, T))
        
    end

end

c1_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, T::sysT{Float64,EX}) = (((AMT(1000*sqrt(6)))/(sqrt(ϕ(gas))*sqrt(vr(vc(gas), v))*abs(3*vr(vc(gas),v) - AMT(1))))*(sqrt(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) + ϕ(gas)*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - ((3*vr(vc(gas), v) - AMT(1))^2)))))*(AMT(1u"m/s"))

function c_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, T::sysT{Float64,EX}, Q) 
    
    if amt((4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) + ϕ(gas)*(4*Tr(Tc(gas), T)*(vr(vc(gas), v)^3) - ((3*vr(vc(gas), v) - AMT(1))^2)))).val < 0
        
        return "Speed of Sound out of Domain"
        
    else
        
        
        if Q == "out"

            return c1_vdw(gas, v, T)
        
        else
        
            ic = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))

            vl = vr1list[ic]*vc(gas)

            vv = vr2list[ic]*vc(gas)

            return c1_vdw(gas, vl, T) + Q*(c1_vdw(gas, vv, T) - c1_vdw(gas, vl, T))
        
        end
        
        
    end
    
end
    
   
# Function with Quality as Argument

function v_vdw(gas::vdWGas, P::sysP{Float64,EX}, Q::_Amt{Float64,EX}, Mol::Bool = false)
    
    if 0 <= amt(Q).val <= 1
        
        if round(amt(Q).val, digits = 2) == 0 #rev4
        
            ic = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))

            return [Pr_sat_list[ic]*Pc(gas), Tr_sat_list[ic]*Tc(gas), vr1list[ic]*vc(gas), u1*amt(Domelist(ϕ(gas), "ur1")[ic]*Pc(gas)*vc(gas)).val, h1*amt(Domelist(ϕ(gas), "hr1")[ic]*Pc(gas)*vc(gas)).val, Domelist(ϕ(gas), "sr1")[ic]*Pc(gas)*vc(gas)/Tc(gas)]
        
        elseif amt(Q).val == 1

            ic = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))

            return [Pr_sat_list[ic]*Pc(gas), Tr_sat_list[ic]*Tc(gas), vr2list[ic]*vc(gas), u1*amt(Domelist(ϕ(gas), "ur2")[ic]*Pc(gas)*vc(gas)).val, h1*amt(Domelist(ϕ(gas), "hr2")[ic]*Pc(gas)*vc(gas)).val, Domelist(ϕ(gas), "sr2")[ic]*Pc(gas)*vc(gas)/Tc(gas)]
            
        else
   
            SatP = findclosest(Pr_sat_list, Pr(Pc(gas), P), (10^-3))

            vrl = vr1list[SatP]

            vrv = vr2list[SatP]

            vr = AMT(vrl) + Q*(vrv - vrl)

            Mol ? vrf = vr*M(gas) : vrf = vr

            return vrf*vc(gas)
            
        end
        
    else
        
        println("Q must be between 0 and 1")
        
    end
    
end

function v_vdw(gas::vdWGas, T::sysT{Float64,EX}, Q::_Amt{Float64,EX}, Mol::Bool = false)
    
    if 0 <= amt(Q).val <= 1
        
        if round(amt(Q).val, digits = 2) == 0 #rev4
        
            ic = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))

            return [Pr_sat_list[ic]*Pc(gas), Tr_sat_list[ic]*Tc(gas), vr1list[ic]*vc(gas), u1*amt(Domelist(ϕ(gas), "ur1")[ic]*Pc(gas)*vc(gas)).val, h1*amt(Domelist(ϕ(gas), "hr1")[ic]*Pc(gas)*vc(gas)).val, Domelist(ϕ(gas), "sr1")[ic]*Pc(gas)*vc(gas)/Tc(gas)]
        
        elseif amt(Q).val == 1

            ic = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))

            return [Pr_sat_list[ic]*Pc(gas), Tr_sat_list[ic]*Tc(gas), vr2list[ic]*vc(gas), u1*amt(Domelist(ϕ(gas), "ur2")[ic]*Pc(gas)*vc(gas)).val, h1*amt(Domelist(ϕ(gas), "hr2")[ic]*Pc(gas)*vc(gas)).val, Domelist(ϕ(gas), "sr2")[ic]*Pc(gas)*vc(gas)/Tc(gas)]
            
        else
        
        SatT = findclosest(Tr_sat_list, Tr(Tc(gas), T), (10^-3))

        vrl = vr1list[SatT]

        vrv = vr2list[SatT]

        vr = AMT(vrl) + Q*(vrv - vrl)

        Mol ? vrf = vr*M(gas) : vrf = vr

        return vrf*vc(gas)
            
        end
        
    else
        
        println("Q must be between 0 and 1")
        
    end
    
end

function T_vdw(gas::vdWGas, v::vAmt{Float64,EX,MA}, Q::_Amt{Float64,EX})
    
    if round(amt(Q).val, digits = 2) == 0 #rev4
        
        ic = findclosest(vr1list, vr(vc(gas), v), (10^-3))
        
        return [Pr_sat_list[ic]*Pc(gas), Tr_sat_list[ic]*Tc(gas), vr1list[ic]*vc(gas), u1*amt(Domelist(ϕ(gas), "ur1")[ic]*Pc(gas)*vc(gas)).val, h1*amt(Domelist(ϕ(gas), "hr1")[ic]*Pc(gas)*vc(gas)).val, Domelist(ϕ(gas), "sr1")[ic]*Pc(gas)*vc(gas)/Tc(gas)]
        
    elseif amt(Q).val == 1
        
        ic = findclosest(vr2list, vr(vc(gas), v), (10^-3))
        
        return [Pr_sat_list[ic]*Pc(gas), Tr_sat_list[ic]*Tc(gas), vr2list[ic]*vc(gas), u1*amt(Domelist(ϕ(gas), "ur2")[ic]*Pc(gas)*vc(gas)).val, h1*amt(Domelist(ϕ(gas), "hr2")[ic]*Pc(gas)*vc(gas)).val, Domelist(ϕ(gas), "sr2")[ic]*Pc(gas)*vc(gas)/Tc(gas)]
        
    else
   
        vre = amt(vr(vc(gas), v)).val

        Q = amt(Q).val

        Point = FindWithQ(vre, Q, vr1list, vr2list)

        return Tr_sat_list[Point]*Tc(gas)
        
    end
    
end
    
function T_vdw(gas::vdWGas, u::uAmt{Float64,EX,MA}, Q::_Amt{Float64,EX})
    
    if round(amt(Q).val, digits = 2) == 0 #rev4
        
        ic = findclosest(Domelist(ϕ(gas), "ur1"), ur(u, Pc(gas), vc(gas)), (10^-3))
        
        return [Pr_sat_list[ic]*Pc(gas), Tr_sat_list[ic]*Tc(gas), vr1list[ic]*vc(gas), u1*amt(Domelist(ϕ(gas), "ur1")[ic]*Pc(gas)*vc(gas)).val, h1*amt(Domelist(ϕ(gas), "hr1")[ic]*Pc(gas)*vc(gas)).val, Domelist(ϕ(gas), "sr1")[ic]*Pc(gas)*vc(gas)/Tc(gas)]
        
    elseif amt(Q).val == 1
        
        ic = findclosest(Domelist(ϕ(gas), "ur2"), ur(u, Pc(gas), vc(gas)), (10^-3))
        
        return [Pr_sat_list[ic]*Pc(gas), Tr_sat_list[ic]*Tc(gas), vr2list[ic]*vc(gas), u1*amt(Domelist(ϕ(gas), "ur2")[ic]*Pc(gas)*vc(gas)).val, h1*amt(Domelist(ϕ(gas), "hr2")[ic]*Pc(gas)*vc(gas)).val, Domelist(ϕ(gas), "sr2")[ic]*Pc(gas)*vc(gas)/Tc(gas)]
        
    else
   
        ure = amt(ur(u, Pc(gas), vc(gas))).val

        Q = amt(Q).val

        Point = FindWithQ(ure, Q, Domelist(ϕ(gas), "ur1"), Domelist(ϕ(gas), "ur2"))

        return Tr_sat_list[Point]*Tc(gas)
        
    end
    
end
    
function T_vdw(gas::vdWGas, h::hAmt{Float64,EX,MA}, Q::_Amt{Float64,EX})
    
    if round(amt(Q).val, digits = 2) == 0 #rev4
        
        ic = findclosest(Domelist(ϕ(gas), "hr1"), hr(h, Pc(gas), vc(gas)), (10^-3))
        
        return [Pr_sat_list[ic]*Pc(gas), Tr_sat_list[ic]*Tc(gas), vr1list[ic]*vc(gas), u1*amt(Domelist(ϕ(gas), "ur1")[ic]*Pc(gas)*vc(gas)).val, h1*amt(Domelist(ϕ(gas), "hr1")[ic]*Pc(gas)*vc(gas)).val, Domelist(ϕ(gas), "sr1")[ic]*Pc(gas)*vc(gas)/Tc(gas)]
        
    elseif amt(Q).val == 1
        
        ic = findclosest(Domelist(ϕ(gas), "hr2"), hr(h, Pc(gas), vc(gas)), (10^-3))
        
        return [Pr_sat_list[ic]*Pc(gas), Tr_sat_list[ic]*Tc(gas), vr2list[ic]*vc(gas), u1*amt(Domelist(ϕ(gas), "ur2")[ic]*Pc(gas)*vc(gas)).val, h1*amt(Domelist(ϕ(gas), "hr2")[ic]*Pc(gas)*vc(gas)).val, Domelist(ϕ(gas), "sr2")[ic]*Pc(gas)*vc(gas)/Tc(gas)]
        
    else
   
        hre = amt(hr(h, Pc(gas), vc(gas))).val

        Q = amt(Q).val

        Point = FindWithQ(hre, Q, Domelist(ϕ(gas), "hr1"), Domelist(ϕ(gas), "hr2"))

        return Tr_sat_list[Point]*Tc(gas)
        
    end
    
end

function T_vdw(gas::vdWGas, s::sAmt{Float64,EX,MA}, Q::_Amt{Float64,EX})
    
    if round(amt(Q).val, digits = 2) == 0 #rev4
        
        ic = findclosest(Domelist(ϕ(gas), "sr1"), sr(s, Pc(gas), vc(gas), Tc(gas)), (10^-3))
        
        return [Pr_sat_list[ic]*Pc(gas), Tr_sat_list[ic]*Tc(gas), vr1list[ic]*vc(gas), u1*amt(Domelist(ϕ(gas), "ur1")[ic]*Pc(gas)*vc(gas)).val, h1*amt(Domelist(ϕ(gas), "hr1")[ic]*Pc(gas)*vc(gas)).val, Domelist(ϕ(gas), "sr1")[ic]*Pc(gas)*vc(gas)/Tc(gas)]
        
    elseif amt(Q).val == 1
        
        ic = findclosest(Domelist(ϕ(gas), "sr2"), sr(s, Pc(gas), vc(gas), Tc(gas)), (10^-3))
        
        return [Pr_sat_list[ic]*Pc(gas), Tr_sat_list[ic]*Tc(gas), vr2list[ic]*vc(gas), u1*amt(Domelist(ϕ(gas), "ur2")[ic]*Pc(gas)*vc(gas)).val, h1*amt(Domelist(ϕ(gas), "hr2")[ic]*Pc(gas)*vc(gas)).val, Domelist(ϕ(gas), "sr2")[ic]*Pc(gas)*vc(gas)/Tc(gas)]
        
    else
   
        sre = amt(sr(s, Pc(gas), vc(gas), Tc(gas))).val

        Q = amt(Q).val

        Point = FindWithQ(sre, Q, Domelist(ϕ(gas), "sr1"), Domelist(ϕ(gas), "sr2"))

        return Tr_sat_list[Point]*Tc(gas)
        
    end
    
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
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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
        
       β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
        return St
        
    elseif (ta == sysT{Float64,EX} && tb == _Amt{Float64,EX}) || (tb == sysT{Float64,EX} && ta == _Amt{Float64,EX})
        
        ta == sysT{Float64,EX} ? T = a : T = b
        
        tb == _Amt{Float64,EX} ? Q = b : Q = a
        
        if round(amt(Q).val, digits = 2) == 0 || amt(Q).val == 1
            
            P = v_vdw(gas, T, Q)[1]
            
            v = v_vdw(gas, T, Q)[3]
            
            u = v_vdw(gas, T, Q)[4]
            
            h = v_vdw(gas, T, Q)[5]
            
            s = v_vdw(gas, T, Q)[6]
            
        else 
        
            v = v_vdw(gas, T, Q)

            P = P_vdw(gas, T, v)

            u = u_vdw(gas, T, v)[1]

            h = h_vdw(gas, T, v)[1]

            s = s_vdw(gas, T, v)[1]
            
        end
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
        return St
        
    elseif (ta == sysP{Float64,EX} && tb == _Amt{Float64,EX}) || (tb == sysP{Float64,EX} && ta == _Amt{Float64,EX})
        
        ta == sysP{Float64,EX} ? P = a : P = b
        
        tb == _Amt{Float64,EX} ? Q = b : Q = a
        
        if round(amt(Q).val, digits = 2) == 0 || amt(Q).val == 1
            
            T = v_vdw(gas, P, Q)[2]
            
            v = v_vdw(gas, P, Q)[3]
            
            u = v_vdw(gas, P, Q)[4]
            
            h = v_vdw(gas, P, Q)[5]
            
            s = v_vdw(gas, P, Q)[6]
            
        else 
        
            v = v_vdw(gas, P, Q)

            T = T_vdw(gas, P, v)

            u = u_vdw(gas, T, v)[1]

            h = h_vdw(gas, T, v)[1]

            s = s_vdw(gas, T, v)[1]
            
        end
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
        return St
        
    elseif (ta == vAmt{Float64,EX,MA} && tb == _Amt{Float64,EX}) || (tb == vAmt{Float64,EX,MA} && ta == _Amt{Float64,EX})
        
        ta == vAmt{Float64,EX,MA} ? v = a : v = b
        
        tb == _Amt{Float64,EX} ? Q = b : Q = a
        
        if round(amt(Q).val, digits = 2) == 0 || amt(Q).val == 1
            
            P = T_vdw(gas, v, Q)[1]
            
            T = T_vdw(gas, v, Q)[2]
            
            u = T_vdw(gas, v, Q)[4]
            
            h = T_vdw(gas, v, Q)[5]
            
            s = T_vdw(gas, v, Q)[6]
            
        else 
        
            T = T_vdw(gas, v, Q)

            P = P_vdw(gas, T, v)

            u = u_vdw(gas, T, v)[1]

            h = h_vdw(gas, T, v)[1]

            s = s_vdw(gas, T, v)[1]
            
        end
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
        return St
        
    elseif (ta == uAmt{Float64,EX,MA} && tb == _Amt{Float64,EX}) || (tb == uAmt{Float64,EX,MA} && ta == _Amt{Float64,EX})
        
        ta == uAmt{Float64,EX,MA} ? u = a : u = b
        
        tb == _Amt{Float64,EX} ? Q = b : Q = a
        
        if round(amt(Q).val, digits = 2) == 0 || amt(Q).val == 1
            
            P = T_vdw(gas, u, Q)[1]
            
            T = T_vdw(gas, u, Q)[2]
            
            v = T_vdw(gas, u, Q)[3]
            
            h = T_vdw(gas, u, Q)[5]
            
            s = T_vdw(gas, u, Q)[6]
            
        else 
        
            T = T_vdw(gas, u, Q)

            v = v_vdw(gas, T, u)[1]

            P = P_vdw(gas, T, v)

            h = h_vdw(gas, T, v)[1]

            s = s_vdw(gas, T, v)[1]
            
        end
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
        return St
        
    elseif (ta == hAmt{Float64,EX,MA} && tb == _Amt{Float64,EX}) || (tb == hAmt{Float64,EX,MA} && ta == _Amt{Float64,EX})
        
        ta == hAmt{Float64,EX,MA} ? h = a : h = b
        
        tb == _Amt{Float64,EX} ? Q = b : Q = a
        
        if round(amt(Q).val, digits = 2) == 0 || amt(Q).val == 1
            
            P = T_vdw(gas, h, Q)[1]
            
            T = T_vdw(gas, h, Q)[2]
            
            v = T_vdw(gas, h, Q)[3]
            
            u = T_vdw(gas, h, Q)[4]
            
            s = T_vdw(gas, h, Q)[6]
            
        else 
        
            T = T_vdw(gas, h, Q)

            v = v_vdw(gas, T, h)[1]

            P = P_vdw(gas, T, v)

            u = u_vdw(gas, T, v)[1]

            s = s_vdw(gas, T, v)[1]
            
        end
        
        a = a_vdw(gas, T, v)
        
        cv = cv_vdw(gas)
        
        cp = cp_vdw(gas, T, v)
        
        γ = gamma(cp, cv)
        
        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
        return St
        
    elseif (ta == sAmt{Float64,EX,MA} && tb == _Amt{Float64,EX}) || (tb == sAmt{Float64,EX,MA} && ta == _Amt{Float64,EX})
        
        ta == sAmt{Float64,EX,MA} ? s = a : s = b
        
        tb == _Amt{Float64,EX} ? Q = b : Q = a
        
        if round(amt(Q).val, digits = 2) == 0 || amt(Q).val == 1
            
            P = T_vdw(gas, s, Q)[1]
            
            T = T_vdw(gas, s, Q)[2]
            
            v = T_vdw(gas, s, Q)[3]
            
            u = T_vdw(gas, s, Q)[4]
            
            h = T_vdw(gas, s, Q)[5]
            
        else    
        
            T = T_vdw(gas, s, Q)

            v = v_vdw(gas, T, s)[1]

            P = P_vdw(gas, T, v)

            h = h_vdw(gas, T, v)[1]

            u = u_vdw(gas, T, v)[1]
            
        end

        a = a_vdw(gas, T, v)

        cv = cv_vdw(gas)

        cp = cp_vdw(gas, T, v)

        γ = gamma(cp, cv)

        β = beta(gas, v, P, Q)
        
        Ks = ks(gas, v, T, Q)
        
        Kt = kt(gas, v, T, Q)
        
        k = k_vdw(gas, v, T, Q)
        
        c = c_vdw(gas, v, T, Q)
        
        Mol ? St = [P, T, v*M(gas), u*M(gas), h*M(gas), s*M(gas), a*M(gas), cv*M(gas), cp*M(gas), γ, β, Ks, Kt, k, c, Q] : St = [P, T, v, u, h, s, a, cv, cp, γ, β, Ks, Kt, k, c, Q]
        
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

#Plots

using Plots

function PrD(vr,n)
    
    if vr < vr1list[n] || vr > vr2list[n]
        
        return 8*Tr_sat_list[n]/(3*vr - 1) - 3/(vr^2)
        
    else
        
        return Pr_sat_list[n]
        
    end
    
end

PrND(vr,Tr) = 8*Tr/(3*vr - 1) - 3/(vr^2)

function TrD(vr,n)
    
    if vr < vr1list[n] || vr > vr2list[n]
        
        return ((Pr_sat_list[n]*(3*vr - 1))/8) + ((3*(3*vr - 1))/(8*(vr^2)))
        
    else
        
        return Tr_sat_list[n]
        
    end
    
end

TrND(vr,Pr) = ((Pr*(3*vr - 1))/8) + ((3*(3*vr - 1))/(8*(vr^2)))

function PlotDome(Tconst::Array = [], T::Bool = false, vpoints::Number = 9.3672917718184e7 ,str::String = "log")
    
    labels = []
    
    x = [vr1list, vr2list]
    
    xlog = [log.(vr1list), log.(vr2list)]
    
    if T == false
        
        Ts = "Tᵣ - "

        data = [Pr_sat_list, Pr_sat_list]
        
        datalog = [log.(Pr_sat_list), log.(Pr_sat_list)]
        
        is = []
        
        if vpoints <= 3000
            
            vx = Array(range(0.3496149073373103, vpoints; length = 50000))
            
        else 
        
            vx = vcat(Array(range(0.3496149073373103, 3000; length = 50000)), Array(range(3000, vpoints; length = 50000)))
            
        end
            
        for n in 1:length(Tconst)
            
            Tv = string(Tconst[n])
            
            if n == 1
            
                labels = [string(Ts,Tv)]
                
            else
                
                labels = hcat(labels, [string(Ts,Tv)])
                
            end

            append!(is, findclosest(Tr_sat_list, AMT(Tconst[n]), (10^-3)))
            
            x = vcat(x, [vx])
            
            xlog = vcat(xlog, [log.(vx)])
            
            if is[n] == -1
            
                data = vcat(data, [PrND.(vx, Tconst[n])])
                
                datalog = vcat(datalog, [log.(PrND.(vx, Tconst[n]))])
                
            else
                
                data = vcat(data, [PrD.(vx, is[n])])
                
                datalog = vcat(datalog, [log.(PrD.(vx, is[n]))])
                
            end

        end
        
        labels == [] ? labels = nothing : nothing
        
        str == "original" ? plot(x[3:end], data[3:end], title = "Pᵣ x vᵣ", label = labels,xlabel = "vᵣ", ylabel = "Pᵣ", markersize = 3, tickfontsize = 6, guidefontsize = 8, legendfontsize = 5, titlefontsize = 9, width = 1, thickness_scaling = 2) :
        
        str == "log" ? plot(xlog[3:end], datalog[3:end], title = "log(Pᵣ) x log(vᵣ)", label = labels, xlabel = "log(vᵣ)", ylabel = "log(Pᵣ)", markersize = 3, tickfontsize = 6, guidefontsize = 8, legendfontsize = 5, titlefontsize = 9, width = 1, thickness_scaling = 2) :
        
        nothing
        
        str == "original" ? plot!(x[1:2], data[1:2], width = 2, color = :black, label = ["Saturation Dome" nothing]) :
        
        str == "log" ? plot!(xlog[1:2], datalog[1:2], width = 2, color = :black, label = ["Saturation Dome" nothing]) :
        
        println("For log(Pᵣ) x log(vᵣ) there is no need for string argument / For Pᵣ x vᵣ use the string original")
        
    else
        
        Ps = "Pᵣ - "
        
        data = [Tr_sat_list, Tr_sat_list]
        
        datalog = [log.(Tr_sat_list), log.(Tr_sat_list)]
        
        is = []
        
        if vpoints <= 3000
            
            vx = Array(range(0.3496149073373103, vpoints; length = 50000))
            
        else 
        
            vx = vcat(Array(range(0.3496149073373103, 3000; length = 50000)), Array(range(3000, vpoints; length = 50000)))
            
        end
            
        for n in 1:length(Tconst)
            
            Pv = string(Tconst[n])
            
            if n == 1
                
                labels = [string(Ps,Pv)]
                
            else
            
                labels = hcat(labels, [string(Ps,Pv)])
                
            end

            append!(is, findclosest(Pr_sat_list, AMT(Tconst[n]), (10^-3)))
            
            x = vcat(x, [vx])
            
            xlog = vcat(xlog, [log.(vx)])
            
            if is[n] == -1
            
                data = vcat(data, [TrND.(vx, Tconst[n])])
                
                datalog = vcat(datalog, [log.(TrND.(vx, Tconst[n]))])
                
            else
                
                data = vcat(data, [TrD.(vx, is[n])])
                
                datalog = vcat(datalog, [log.(TrD.(vx, is[n]))])
                
            end

        end
        
        labels == [] ? labels = nothing : nothing
        
        str == "original" ? plot(x[3:end], data[3:end], title = "Tᵣ x vᵣ",label = labels, xlabel = "vᵣ", ylabel = "Tᵣ", markersize = 3, tickfontsize = 6, guidefontsize = 8, legendfontsize = 5, titlefontsize = 9, width = 1, thickness_scaling = 2) :
        
        str == "log" ? plot(xlog[3:end], datalog[3:end], title = "log(Tᵣ) x log(vᵣ)",label = labels, xlabel = "log(vᵣ)", ylabel = "log(Tᵣ)", markersize = 3, tickfontsize = 6, guidefontsize = 8, legendfontsize = 5, titlefontsize = 9, width = 1, thickness_scaling = 2) :
        
        nothing
        
        str == "original" ? plot!(x[1:2], data[1:2], width = 2, color = :black, label = ["Saturation Dome" nothing]) :
        
        str == "log" ? plot!(xlog[1:2], datalog[1:2], width = 2, color = :black, label = ["Saturation Dome" nothing]) :
        
        println("For log(Pᵣ) x log(vᵣ) there is no need for string argument / For Pᵣ x vᵣ use the string original")
        
    end

end    

function PlotDome(gas::vdWGas, Tconst::Array = [], T::Bool = false, vpoints::Number = 9.3672917718184e7 ,str::String = "log")
    
    vcr = amt(vc(gas)).val
    
    Tcr = amt(Tc(gas)).val
    
    Pcr = amt(Pc(gas)).val
    
    v1 = vr1list*vcr
    
    v2 = vr2list*vcr
    
    Psat = Pr_sat_list*Pcr
    
    Tsat = Tr_sat_list*Tcr

    labels = []
    
    x = [v1, v2]
    
    xlog = [log.(v1), log.(v2)]
    
    if T == false
        
        Ts = "T - "

        data = [Psat, Psat]
        
        datalog = [log.(Psat), log.(Psat)]
        
        is = []
        
        if vpoints <= 3000
            
            vx = Array(range(0.3496149073373103, vpoints; length = 50000))
            
        else 
        
            vx = vcat(Array(range(0.3496149073373103, 3000; length = 50000)), Array(range(3000, vpoints; length = 50000)))
            
        end
            
        for n in 1:length(Tconst)
            
            Tv = string(Tconst[n])
            
            if n == 1
            
                labels = [string(Ts,Tv,"K")]
                
            else
                
                labels = hcat(labels, [string(Ts,Tv,"K")])
                
            end

            append!(is, vdWProp.findclosest(Tr_sat_list, AMT(Tconst[n]/Tcr), (10^-3)))
            
            x = vcat(x, [vx*vcr])
            
            xlog = vcat(xlog, [log.(vx*vcr)])
            
            if is[n] == -1
            
                data = vcat(data, [Pcr*PrND.(vx, Tconst[n]/Tcr)])
                
                datalog = vcat(datalog, [log.(Pcr*PrND.(vx, Tconst[n]/Tcr))])
                
            else
                
                data = vcat(data, [Pcr*PrD.(vx, is[n])])
                
                datalog = vcat(datalog, [log.(Pcr*PrD.(vx, is[n]))])
                
            end

        end
        
        labels == [] ? labels = nothing : nothing
        
        str == "original" ? plot(x[3:end], data[3:end], title = "P x v", label = labels,xlabel = "v[m³/kg]", ylabel = "P[kPa]", markersize = 3, tickfontsize = 6, guidefontsize = 8, legendfontsize = 5, titlefontsize = 9, width = 1, thickness_scaling = 2) :
        
        str == "log" ? plot(xlog[3:end], datalog[3:end], title = "log(P) x log(v)", label = labels, xlabel = "log(v[m³/kg])", ylabel = "log(P[kPa])", markersize = 3, tickfontsize = 6, guidefontsize = 8, legendfontsize = 5, titlefontsize = 9, width = 1, thickness_scaling = 2) :
        
        nothing
        
        str == "original" ? plot!(x[1:2], data[1:2], width = 2, color = :black, label = ["Saturation Dome" nothing]) :
        
        str == "log" ? plot!(xlog[1:2], datalog[1:2], width = 2, color = :black, label = ["Saturation Dome" nothing]) :
        
        println("For log(P) x log(v) there is no need for string argument / For P x v use the string original")
        
    else
        
        Ps = "P - "
        
        data = [Tsat, Tsat]
        
        datalog = [log.(Tsat), log.(Tsat)]
        
        is = []
        
        #0.3697942021696152
        
        if vpoints <= 3000
            
            vx = Array(range(0.3496149073373103, vpoints; length = 50000))
            
        else 
        
            vx = vcat(Array(range(0.3496149073373103, 3000; length = 50000)), Array(range(3000, vpoints; length = 50000)))
            
        end
            
        for n in 1:length(Tconst)
            
            Pv = string(Tconst[n])
            
            if n == 1
                
                labels = [string(Ps,Pv,"kPa")]
                
            else
            
                labels = hcat(labels, [string(Ps,Pv,"kPa")])
                
            end

            append!(is, vdWProp.findclosest(Pr_sat_list, AMT(Tconst[n]/Pcr), (10^-3)))
            
            x = vcat(x, [vx*vcr])
            
            xlog = vcat(xlog, [log.(vcr*vx)])
            
            if is[n] == -1
            
                data = vcat(data, [Tcr*TrND.(vx, Tconst[n]/Pcr)])
                
                datalog = vcat(datalog, [log.(Tcr*TrND.(vx, Tconst[n]/Pcr))])
                
            else
                
                data = vcat(data, [Tcr*TrD.(vx, is[n])])
                
                datalog = vcat(datalog, [log.(Tcr*TrD.(vx, is[n]))])
                
            end

        end
        
        labels == [] ? labels = nothing : nothing
        
        str == "original" ? plot(x[3:end], data[3:end], title = "T x v",label = labels, xlabel = "v[m³/kg]", ylabel = "T[K]", markersize = 3, tickfontsize = 6, guidefontsize = 8, legendfontsize = 5, titlefontsize = 9, width = 1, thickness_scaling = 2) :
        
        str == "log" ? plot(xlog[3:end], datalog[3:end], title = "log(T) x log(v)",label = labels, xlabel = "log(v[m³/kg])", ylabel = "log(T[K])", markersize = 3, tickfontsize = 6, guidefontsize = 8, legendfontsize = 5, titlefontsize = 9, width = 1, thickness_scaling = 2) :
        
        nothing
        
        str == "original" ? plot!(x[1:2], data[1:2], width = 2, color = :black, label = ["Saturation Dome" nothing]) :
        
        str == "log" ? plot!(xlog[1:2], datalog[1:2], width = 2, color = :black, label = ["Saturation Dome" nothing]) :
        
        println("For log(T) x log(v) there is no need for string argument / For T x v use the string original")
        
    end
    
end    

export PlotDome

end #module