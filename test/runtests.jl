using PaperUnfittedFEMDarcy
using Test

@testset "FIG8-1_h-convergence-on-a-cut-square.jl" begin
    include("FIG8-1_h-convergence-on-a-cut-square.jl")
end

@testset "FIG8-2_sensitivity-with-respect-to-gamma.jl" begin
    include("FIG8-2_sensitivity-with-respect-to-gamma.jl")
end

@testset "FIG8-3_pressure-robustness.jl" begin
    include("FIG8-3_pressure-robustness.jl")
end

@testset "FIG8-4_conditioning-with-respect-to-cut-cell-length.jl" begin
    include("FIG8-4_conditioning-with-respect-to-cut-cell-length.jl")
end


