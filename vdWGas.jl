using EngThermBase

# criacao da estrutura da substancia que sera utilizada como argumento

struct vdWGas
    nome::String
    Pc::AbstractFloat
    Tc::AbstractFloat
    α::AbstractFloat
    b::AbstractFloat
    M::AbstractFloat
end

nome(gas::vdWGas) = gas.nome
Pc(gas::vdWGas) = gas.Pc
Tc(gas::vdWGas) = gas.Tc
α(gas::vdWGas) = gas.α
b(gas::vdWGas) = gas.b
M(gas::vdWGas) = gas.M

# Funcoes de conversao de unidades

convertP(P::Float64) = AbstractFloat(P*101.325)
convertM(M::Float64) = AbstractFloat(M*0.001)
convertv(v::Float64,M::AbstractFloat) = AbstractFloat(v*(10^(-6))/M)
convertα(α::Float64,M::AbstractFloat) = AbstractFloat(α*101.325*0.000001/(M^2))
convertb(b::Float64,M::AbstractFloat) = AbstractFloat(b*(10^(-2))*0.001/M)

# Funcao para dar a estrutura nas unidades convertidas

function vdWGasC(nome::String,Pc::AbstractFloat,Tc::AbstractFloat,α::AbstractFloat,b::AbstractFloat,M::AbstractFloat)
    
    Pcon = convertP(Pc)
    
    Mcon = convertM(M)
    
    αcon = convertα(α,Mcon)
    
    bcon = convertb(b,Mcon)
    
    return vdWGas(nome,Pcon,Tc,αcon,bcon,Mcon)
    
end