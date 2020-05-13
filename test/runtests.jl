using Test, vdWProp, EngThermBase, Polynomials, Roots

@testset "vdWProp" begin
    
    pP = vdWProp.IsoProp(vdWProp.Hg, P(1000), T(1500), T(1650), "P")
    
    pT = vdWProp.IsoProp(vdWProp.Hg, P(1000), T(1500), P(1200), "T")
    
    pv = vdWProp.IsoProp(vdWProp.Hg, P(1000), T(1500), T(1650), "v")
    
    pu = vdWProp.IsoProp(vdWProp.Hg, P(1000), T(1500), T(1650), "u")
    
    ph = vdWProp.IsoProp(vdWProp.Hg, P(1000), T(1500), T(1650), "h")
    
    ps = vdWProp.IsoProp(vdWProp.Hg, P(1000), T(1500), T(1650), "s")
    
    @test vdWProp.isoP(vdWProp.Hg, pP[2,1], pP[3,1], pP[2,2], pP[3,2]) < 0.000001
    
    @test vdWProp.isoT(vdWProp.Hg, pT[1,1], pT[3,1], pT[1,2], pT[3,2]) < 0.000001
    
    @test vdWProp.isou(vdWProp.Hg, pu[2,1], pu[3,1], pu[2,2], pu[3,2]) < 0.000001
    
    @test vdWProp.isoh(vdWProp.Hg, ph[2,1], ph[3,1], ph[2,2], ph[3,2]) < 0.000001
    
    @test vdWProp.isos(vdWProp.Hg, ps[2,1], ps[3,1], ps[2,2], ps[3,2]) < 0.000001

end