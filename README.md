# MzXML

[![codecov.io](http://codecov.io/github/timholy/MzXML.jl/coverage.svg?branch=master)](http://codecov.io/github/timholy/MzXML.jl?branch=master)

A Julia package for reading mass spectrometry [mzXML files](https://en.wikipedia.org/wiki/Mass_spectrometry_data_format).

Example:

```julia
julia> using MzXML

julia> cd(joinpath(dirname(dirname(pathof(MzXML))), "test"))  # location of a test file

julia> scans, info = MzXML.load("test32.mzXML");

julia> info
Dict{Symbol,Any} with 3 entries:
  :endTime   => 24.156 s
  :msModel   => "API 3000"
  :startTime => 0.004 s

julia> scans
2-element Array{MzXML.MSscan{Float32,Float32},1}:
 └─basePeak 321.4: 0.004 s: Pair{Float32,Float32}[111.3 => 6251.25, 111.4 => 0.0, 166.7 => 0.0, 166.8 => 18753.75, 166.9 => 0.0, 189.1 => 0.0, 189.2 => 12502.5, 189.3 => 0.0, 191.7 => 0.0, 191.8 => 37507.5  …  421.1 => 12502.5, 421.2 => 0.0, 437.9 => 0.0, 438.0 => 6251.25, 438.1 => 0.0, 438.5 => 0.0, 438.6 => 6251.25, 438.7 => 0.0, 481.0 => 0.0, 481.1 => 18753.75]

 └─basePeak 321.0: 0.809 s: Pair{Float32,Float32}[110.9 => 12502.5, 111.0 => 0.0, 186.9 => 0.0, 187.0 => 6251.25, 187.1 => 0.0, 199.5 => 0.0, 199.6 => 12502.5, 199.7 => 18753.75, 199.8 => 50010.0, 199.9 => 6251.25  …  507.4 => 18753.75, 507.5 => 0.0, 526.4 => 0.0, 526.5 => 6251.25, 526.6 => 0.0, 545.6 => 0.0, 545.7 => 6251.25, 545.8 => 0.0, 587.5 => 0.0, 587.6 => 18753.75]
```

The "└─" indicates MS² or higher on the specified base peak. Also shown is the retention time (in seconds), and `mz=>ion count` pairs.
