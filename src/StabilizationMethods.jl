#======================================================================# 
#                          Stabilization types                         #
#======================================================================#

abstract type StabilizationType end

#----------------------------------------------------------------------# 
#                         Bulk Ghost Penalty Stabilization             #
#----------------------------------------------------------------------#

abstract type BulkGhostPenaltyStabilization <: StabilizationType end

"""
BulkGhostPenaltyStabilization()
...
# Fields
- `name::String`: name of the stabilization.
- `ξ::Float64`:  global scaling parameter for all (non-zero) stabilization terms
- `ξᵤ::Float64`: parameter weighting the u-stabilization term. 
- `ξₚ::Float64`: parameter weighting the p-stabilization term.
- `ξₛ::Float64`: parameter weighting the (∇⋅u)-stabilization term.
- `ξₘ::Float64`: parameter weighting the mixed-stabilization term.
- `ξ_AL::Float64`: parameter weighting the Augmented-Lagrangian stabilization term.
- `ζ::Bool`: form parameter: test functions in the stabilization terms are used in difference form (operation_v(v) - proj_op_v) (ζ=true) or full form (operation_v(v)) (ζ=false).
...
Note: Use `0.0` for the parameters ξᵤ, ξₚ, ξₛ and ξₘ if that particular stabilization term needs to be switched off. 
"""

struct BulkGhostPenalty <: BulkGhostPenaltyStabilization
    name::String
    ξ::Float64
    ξᵤ::Float64
    ξₚ::Float64
    ξₛ::Float64
    ξₘ::Float64
    ξ_AL::Float64
    ζ::Bool
end

function NoBGP(;
    name="NoBGP",
    ξ=0.0,
    ξᵤ=0.0,
    ξₚ=0.0,
    ξₛ=0.0,
    ξₘ=0.0,
    ξ_AL=0.0,
    ζ=true) 
    BulkGhostPenalty(name,ξ,ξᵤ,ξₚ,ξₛ,ξₘ,ξ_AL,ζ)
end

# Mix-based stabilization methods:
function BGP(;
    name="BGP",
    ξ=1.0,
    ξᵤ=1.0,
    ξₚ=0.0,
    ξₛ=0.0,
    ξₘ=1.0,
    ξ_AL=0.0,
    ζ=true) 
    BulkGhostPenalty(name,ξ,ξᵤ,ξₚ,ξₛ,ξₘ,ξ_AL,ζ)
end

function ALBGP(;
    name="ALBGP",
    ξ=1.0,
    ξᵤ=1.0,
    ξₚ=0.0,
    ξₛ=0.0,
    ξₘ=1.0,
    ξ_AL=1.0,
    ζ=true) 
    BulkGhostPenalty(name,ξ,ξᵤ,ξₚ,ξₛ,ξₘ,ξ_AL,ζ)
end