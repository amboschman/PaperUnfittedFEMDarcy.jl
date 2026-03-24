module FIG8_2_sensitivity_with_respect_to_gamma

#======================================================================# 
#               FIG 8.2 :   Sensitivity with respect to γ₀             #
#======================================================================#
#= 
This script reproduces the sensitivity analysis with respect to the penalty parameter γ on a cut square using Γ = Γᵤ as shown in Figure 8.2 of the paper. The used parameter settings are further detailed in the comments. To run the sensitivity analysis, the user has to select (i) the level of mesh refinement, (ii) the ratio h_cut/h and (iii) the penalty parameter γ₀ (referred to as γ in the paper). 
=#

using PaperUnfittedFEMDarcy
using Printf

#----------------------------------------------------------------------# 
#                       Fixed Parameter Selection                      #
#----------------------------------------------------------------------#

# Problem parameters
bnd_u, bnd_p = [5,6,7,8], []   # boundary conditions: here Γ = Γᵤ

# Method parameters
topo = Main.PaperUnfittedFEMDarcy.Gridap.QUAD   # mesh topology: paper results use QUAD
k = 0          # polynomial degree: paper results uses 0
contra_variant_piola_map_type=Main.PaperUnfittedFEMDarcy.Gridap.ReferenceFEs.ScaledContraVariantPiolaMap()  # rescaling of the Piola map to improve conditioning: paper results use this option.
PETSc = false  # paper results uses true. Machine-specific installation of PETSc may be required to use this option. If PETSc is not available, set to false to use the default linear solver in Gridap.
m = 0          # symmetry parameter for the problem formulation: paper results use 0 (non-symmetric formulation). 
s = -1         # exponent of penalty parameter γ=γ₀⋅h^s: paper results use -1, i.e., scaling with 1/h. Because bnd_u =[], the chosen value does not affect the solution. 

# Stabilisation parameters
ζ = true    # stabilisation uses the difference of the test function and its projection: paper results use true. 
stab_std = BGP(ξ=0.0,ξᵤ=0.0,ξₚ=0.0,ξₛ=0.0,ξₘ=0.0,ξ_AL=0.0,ζ=ζ)   # std method: uses ξ=0.0 to switch off all stabilisation terms.
stab_BGP = BGP(ξ=1.0e0,ξᵤ=1.0,ξₚ=0.0,ξₛ=1.0,ξₘ=1.0,ξ_AL=0.0,ζ=ζ) # BGP method. Note that τ₀= ξᵤ = ξₛ and τ_d = ξₘ in the paper.
stab_ALBGP = BGP(ξ=1.0e0,ξᵤ=1.0,ξₚ=0.0,ξₛ=1.0,ξₘ=1.0,ξ_AL=1.0,ζ=ζ) # AL-BGP method. Note that τ₀= ξᵤ = ξₛ, τ_d = ξₘ and τ_AL = ξ_AL in the paper.

# Analysis parameters
analysis = UnfittedAnalysis(linf_mass=false,linf_res=false)

#----------------------------------------------------------------------# 
#                       Run h-convergence Test                         #
#----------------------------------------------------------------------#

# h-convergence test parameters
ns = [8,16] # number of elements: paper results use [8,16,32,64,128,256,512]
ϕ  = 0.5    # fraction h_cut/h: paper results use [0.5,0.5e-6]
γ₀ = 1.0e0  # unscaled value of penalty parameter: paper results use [1.0e0,1.0e2,1.0e4] (and refers to this as γ).

# Run 'std' unfitted method (no stabilisation)
println(" \n --- 'std' results using h_cut/h = $ϕ and γ₀ = $γ₀ --- \n \n| \t h \t \t | \t ||u-u_h||_L2 \t | \t ||p-p_h||_L2 \t | \t ||div(u-u_h)||_L2 \t |")
for n in ns 
    ε = ϕ*(1.0/(n-2))  # cut cell size
    problem = SmoothSolenoidalDarcyOnCutSquare(n=n,bnd_u=bnd_u,bnd_p=bnd_p,ε=ε)
    method  = UnfittedRT(topo=topo,k=k,contra_variant_piola_map_type=contra_variant_piola_map_type,PETSc=PETSc,m=m,γ₀=γ₀,s=s)
    results_std = main(problem,method,stab_std,analysis) 
    println(@sprintf("| \t %.2e \t | \t %.2e \t | \t %.2e \t | \t %.2e \t \t |", results_std[:h], results_std[:err_u_L2], results_std[:err_p_L2], results_std[:err_divu_L2]))
end

# Run "BGP" method
println(" \n --- 'BGP' results using h_cut/h = $ϕ and γ₀ = $γ₀ --- \n \n| \t h \t \t | \t ||u-u_h||_L2 \t | \t ||p-p_h||_L2 \t | \t ||div(u-u_h)||_L2 \t |")
for n in ns 
    ε = ϕ*(1.0/(n-2))  # cut cell size
    problem = SmoothSolenoidalDarcyOnCutSquare(n=n,bnd_u=bnd_u,bnd_p=bnd_p,ε=ε)
    method  = UnfittedRT(topo=topo,k=k,contra_variant_piola_map_type=contra_variant_piola_map_type,PETSc=PETSc,m=m,γ₀=γ₀,s=s)
    results_BGP = main(problem,method,stab_BGP,analysis) 
    println(@sprintf("| \t %.2e \t | \t %.2e \t | \t %.2e \t | \t %.2e \t \t |", results_BGP[:h], results_BGP[:err_u_L2], results_BGP[:err_p_L2], results_BGP[:err_divu_L2]))
end

# Run "AL-BGP" method
println(" \n --- 'AL-BGP' results using h_cut/h = $ϕ and γ₀ = $γ₀ --- \n \n| \t h \t \t | \t ||u-u_h||_L2 \t | \t ||p-p_h||_L2 \t | \t ||div(u-u_h)||_L2 \t |")
for n in ns 
    ε = ϕ*(1.0/(n-2))  # cut cell size
    problem = SmoothSolenoidalDarcyOnCutSquare(n=n,bnd_u=bnd_u,bnd_p=bnd_p,ε=ε)
    method  = UnfittedRT(topo=topo,k=k,contra_variant_piola_map_type=contra_variant_piola_map_type,PETSc=PETSc,m=m,γ₀=γ₀,s=s)
    results_ALBGP = main(problem,method,stab_ALBGP,analysis) 
    println(@sprintf("| \t %.2e \t | \t %.2e \t | \t %.2e \t | \t %.2e \t \t |", results_ALBGP[:h], results_ALBGP[:err_u_L2], results_ALBGP[:err_p_L2], results_ALBGP[:err_divu_L2]))
end

end # module