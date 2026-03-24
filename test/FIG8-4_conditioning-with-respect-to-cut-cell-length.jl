module FIG8_4_conditioning_with_respect_to_cut_cell_length

#======================================================================# 
#     FIG 8.4 :   Conditioning with respect to cut cell length         # 
#======================================================================#
#= 
This script reproduces the conditioning analysis with respect to the cut-cell ratio h_cut/h using Γ = Γᵤ∪Γₚ as shown in Figure 8.4 of the paper. The used parameter settings are further detailed in the comments. To run the conditioning analysis, the user has to specify the cut-cell ratios h_cut/h and (ii) the mesh refinement level. 
=#

using PaperUnfittedFEMDarcy
using Printf

#----------------------------------------------------------------------# 
#                       Fixed Parameter Selection                      #
#----------------------------------------------------------------------#

# Problem parameters
bnd_u, bnd_p = [7,8], [5,6]   # boundary conditions: here Γ = Γᵤ∪Γₚ
n            = 8              # number of uncut elements fixing the mesh refinement level: paper results use 32. Here a lower value is selected to reduce the runtime of the test, as computing the 1-cond number can be expensive for (i) fine meshes and (ii) in case no stabilisation is applied as a result of small cut cells.

# Method parameters
topo = Main.PaperUnfittedFEMDarcy.Gridap.QUAD   # mesh topology: paper results use QUAD
k = 0          # polynomial degree: paper results uses 0
contra_variant_piola_map_type=Main.PaperUnfittedFEMDarcy.Gridap.ReferenceFEs.ScaledContraVariantPiolaMap()  # rescaling of the Piola map to improve conditioning: paper results use this option.
PETSc = false  # paper results uses true. Machine-specific installation of PETSc may be required to use this option. If PETSc is not available, set to false to use the default linear solver in Gridap.
γ₀ = 1.0e0      # unscaled value of penalty parameter: paper results use 0.0e0 (and refers to this as γ). Because bnd_u =[], the chosen value does not affect the solution. 
m = 0          # symmetry parameter for the problem formulation: paper results use 0 (non-symmetric formulation). 
s = -1         # exponent of penalty parameter γ=γ₀⋅h^s: paper results use -1, i.e., scaling with 1/h. Because bnd_u =[], the chosen value does not affect the solution. 
method  = UnfittedRT(topo=topo,k=k,contra_variant_piola_map_type=contra_variant_piola_map_type,PETSc=PETSc,m=m,γ₀=γ₀,s=s)

# Stabilisation parameters
ζ = true    # stabilisation uses the difference of the test function and its projection: paper results use true. 
stab_std = BGP(ξ=0.0,ξᵤ=0.0,ξₚ=0.0,ξₛ=0.0,ξₘ=0.0,ξ_AL=0.0,ζ=ζ)   # std method: uses ξ=0.0 to switch off all stabilisation terms.
stab_BGP = BGP(ξ=1.0e0,ξᵤ=1.0,ξₚ=0.0,ξₛ=1.0,ξₘ=1.0,ξ_AL=0.0,ζ=ζ) # BGP method. Note that τ₀= ξᵤ = ξₛ and τ_d = ξₘ in the paper.
stab_ALBGP = BGP(ξ=1.0e0,ξᵤ=1.0,ξₚ=0.0,ξₛ=1.0,ξₘ=1.0,ξ_AL=1.0,ζ=ζ) # AL-BGP method. Note that τ₀= ξᵤ = ξₛ, τ_d = ξₘ and τ_AL = ξ_AL in the paper.

# Analysis parameters
analysis = UnfittedAnalysis(cond1_A=true,linf_mass=false,linf_res=false)

#----------------------------------------------------------------------# 
#                       Run h-convergence Test                         #
#----------------------------------------------------------------------#

# h-convergence test parameters
ϕs = [0.5e0,0.5e-1,0.5e-9] # fraction h_cut/h: paper results use [0.5e0,0.5e-1,0.5e-2,0.5e-3,0.5e-4,0.5e-5,0.5e-6,0.5e-7,0.5e-8,0.5e-9]

# Run 'std' unfitted method (no stabilisation)
println(" \n --- 'std' results --- \n \n| \t h_cut/h \t | \t Cond. nr \t | \t ||u-u_h||_L2 \t | \t ||p-p_h||_L2 \t | \t ||div(u-u_h)||_L2 \t |")
for ϕ in ϕs 
    ε = ϕ*(1.0/(n-2))  # cut cell size
    problem = SmoothSolenoidalDarcyOnCutSquare(n=n,bnd_u=bnd_u,bnd_p=bnd_p,ε=ε)
    results_std = main(problem,method,stab_std,analysis) 
    println(@sprintf("| \t %.1e \t | \t %.2e \t | \t %.2e \t | \t %.2e \t | \t %.2e \t \t |", ϕ, results_std[:cond1_A], results_std[:err_u_L2], results_std[:err_p_L2], results_std[:err_divu_L2]))
end

# Run "BGP" method
println(" \n --- 'BGP' results --- \n \n| \t h_cut/h \t | \t Cond. nr \t | \t ||u-u_h||_L2 \t | \t ||p-p_h||_L2 \t | \t ||div(u-u_h)||_L2 \t |")
for ϕ in ϕs 
    ε = ϕ*(1.0/(n-2))  # cut cell size
    problem = SmoothSolenoidalDarcyOnCutSquare(n=n,bnd_u=bnd_u,bnd_p=bnd_p,ε=ε)
    results_BGP = main(problem,method,stab_BGP,analysis) 
    println(@sprintf("| \t %.1e \t | \t %.2e \t | \t %.2e \t | \t %.2e \t | \t %.2e \t \t |", ϕ, results_BGP[:cond1_A], results_BGP[:err_u_L2], results_BGP[:err_p_L2], results_BGP[:err_divu_L2]))
end

# Run "AL-BGP" method
println(" \n --- 'AL-BGP' results --- \n \n| \t h_cut/h \t | \t Cond. nr \t | \t ||u-u_h||_L2 \t | \t ||p-p_h||_L2 \t | \t ||div(u-u_h)||_L2 \t |")
for ϕ in ϕs 
    ε = ϕ*(1.0/(n-2))  # cut cell size
    problem = SmoothSolenoidalDarcyOnCutSquare(n=n,bnd_u=bnd_u,bnd_p=bnd_p,ε=ε)
    results_ALBGP = main(problem,method,stab_ALBGP,analysis) 
    println(@sprintf("| \t %.1e \t | \t %.2e \t | \t %.2e \t | \t %.2e \t | \t %.2e \t \t |", ϕ, results_ALBGP[:cond1_A], results_ALBGP[:err_u_L2], results_ALBGP[:err_p_L2], results_ALBGP[:err_divu_L2]))
end

end # module