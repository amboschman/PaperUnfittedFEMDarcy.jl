module PaperUnfittedFEMDarcy

# Gridap 
using Gridap
using GridapEmbedded
using GridapPETSc  #v0.5.5

# Manufactured solutions
include("./DarcySolutions.jl")
export SmoothSolenoidalDarcy2D_u, SmoothSolenoidalDarcy2D_p, SmoothSolenoidalDarcy2D_f, SmoothSolenoidalDarcy2D_g, SmoothSolenoidalDarcy2D_kinv, SmoothSolenoidalDarcy2D_name
export PresRobustDarcy2D_u, PresRobustDarcy2D_p, PresRobustDarcy2D_f, PresRobustDarcy2D_g, PresRobustDarcy2D_kinv, PresRobustDarcy2D_name

# Methods
include("./Methods.jl")
export MethodType
export UnfittedMethod
export Unfitted, UnfittedRT

# Stabilization
include("./StabilizationMethods.jl")
export StabilizationType
export BulkGhostPenaltyStabilization
export BulkGhostPenalty, NoBGP, BGP, ALBGP

# Problems
include("./DarcyProblems.jl")
export ProblemType
export DarcyProblem
export DarcyOnCutSquare
export SmoothSolenoidalDarcyOnCutSquare, PresRobustDarcyOnCutSquare
export setup_geometry

# Postprocessing/Analysis 
using LinearAlgebra: cond
include("./Analysis.jl")
export AnalysisMethod
export UnfittedAnalysis
export setup_results, compute_linf_mass_err, compute_linf_res

# Drivers
include("./Drivers.jl")
export main

end
