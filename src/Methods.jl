#======================================================================# 
#                          Method types                                #
#======================================================================#

abstract type MethodType end

#----------------------------------------------------------------------# 
#                         Unfitted Methods                             #
#----------------------------------------------------------------------#

abstract type UnfittedMethod <: MethodType end

"""
Unfitted()
...
# Fields
- `name::String`: name of the method.
- `shape_u::DataType`: shape function for velocity field u.
- `shape_p::DataType`: shape function for pressure field p.
- `k::Int`: polynomial order for velocity field basis.
- `degree::Int`: degree for integration.
- `topo::Polytope`: Polytope used for discretization.
- `contra_variant_piola_map_type::Gridap.ReferenceFEs.ContraVariantPiolaMapType`: type of contravariant Piola map
- `PETSc::Bool`: switch to use GridapPETSc for the linear solver
- `m::Int`: symmetry parameter.
- `γ₀::Float`: unscaled value of Nitsche parameter.
- `s::Float`: scaling of Nitsche parameter with edge length to the power `s`.
...
Note that `m` can be used to distinguish between two different cases, that is:
(N.1) Symmetric Nitsche formulation for `m=1`. This is the default option; 
(N.2) Non-symmetric Nitsche formulation for `m=0`.
"""
struct Unfitted <: UnfittedMethod
    name::String
    shape_u::DataType
    shape_p::DataType
    k::Int
    degree::Int
    topo::Polytope
    contra_variant_piola_map_type::Gridap.ReferenceFEs.ContraVariantPiolaMapType
    PETSc::Bool
    m::Int
    γ₀::Float64
    s::Int
end

function UnfittedRT(;
    name    = "UnfittedRT",
    shape_u = RaviartThomas,
    shape_p = Lagrangian,
    k       = 0,
    degree  = 2*2*(k+1),
    topo    = QUAD,
    contra_variant_piola_map_type = Gridap.ReferenceFEs.ContraVariantPiolaMap(),
    PETSc   = false,
    m       = 0,
    γ₀      = 1e3,
    s       = -1)
    Unfitted(name,shape_u,shape_p,k,degree,topo,contra_variant_piola_map_type,PETSc,m,γ₀,s)
end