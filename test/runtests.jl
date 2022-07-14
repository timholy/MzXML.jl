using MzXML, Unitful
using MzXML.IntervalSets
using MzCore
using Test

# The test96.mzXML file comes from https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000089747

@testset "MzXML" begin
    scans, info = MzXML.load("test32.mzXML")
    @test info[:msModel] == "API 3000"
    @test eltype(scans) == MzXML.MSscan{Float32,Float32}
    @test length(scans) == 2
    scan = scans[1]
    @test scan.polarity == '-'
    @test scan.msLevel == 2
    @test scan.basePeakMz == 321.4
    @test scan.totIonCurrent == 1.637827475594e06
    @test scan.retentionTime == 0.004*u"s"
    @test scans[2].retentionTime == 0.809*u"s"

    # MzCore traits
    @test intensitytype(scans[1]) === Float32
    @test mztype(scans[1]) === Float32
    axmz, axt = limits(scans)
    @test minimum(axmz.val) === 110.9f0
    @test maximum(axmz.val) === 587.6f0
    @test minimum(axt.val) === 0.004*u"s"
    @test maximum(axt.val) === 0.809*u"s"

    scans, info = MzXML.load("test32.mzXML"; timeinterval=0.2u"s" .. Inf*u"s")
    @test length(scans) == 1 && scans[1].retentionTime == 0.809*u"s"

    io = IOBuffer()
    show(io, scan)
    @test occursin("└─basePeak 321.4", String(take!(io)))

    scans, info = MzXML.load("test32.mzXML"; timeinterval=10.2u"s" .. Inf*u"s", timeshift=10u"s")
    @test length(scans) == 1 && scans[1].retentionTime == 10.809*u"s"

    idx = MzXML.index("test32.mzXML")
    @test idx == [1231, 4242]

    scans, info = MzXML.load("test64.mzXML")
    @test eltype(scans) == MzXML.MSscan{Float64,Float64}
    @test intensitytype(scans[1]) === Float64
    @test mztype(scans[1]) === Float64

    scans, info = MzXML.load("test96.mzXML")
    @test eltype(scans) == MzXML.MSscan{Float32,Float32}
    @test intensitytype(scans[1]) === Float32
    @test mztype(scans[1]) === Float32
end
