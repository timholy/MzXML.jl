using mzXML, Unitful
using Test

@testset "mzXML" begin
    scans, info = mzXML.load("test32.mzXML")
    @test info[:msModel] == "API 3000"
    @test eltype(scans) == mzXML.MSscan{Float32,Float32}
    @test length(scans) == 2
    scan = scans[1]
    @test scan.polarity == '-'
    @test scan.msLevel == 2
    @test scan.basePeakMz == 321.4
    @test scan.totIonCurrent == 1.637827475594e06
    @test scan.retentionTime == 0.004*u"s"
    @test scans[2].retentionTime == 0.809*u"s"

    idx = mzXML.index("test32.mzXML")
    @test idx == [1231, 4242]

    scans, info = mzXML.load("test64.mzXML")
    @test eltype(scans) == mzXML.MSscan{Float64,Float64}
end
