#======================================================================# 
#                            Analysis                                  #
#======================================================================#

abstract type AnalysisMethod end

# Unfitted analysis
struct UnfittedAnalysis <: AnalysisMethod
    name::String
    geometry::Bool
    vtu::Bool
    results::Bool
    linf_mass::Bool
    linf_res::Bool
    cond2_A::Bool
    cond1_A::Bool
end

function UnfittedAnalysis(;
    name     = "UnfittedAnalysis",
    geometry = false,
    vtu      = false, 
    results  = true,
    linf_mass = true,
    linf_res  = true,
    cond2_A   = false,
    cond1_A   = false
  )
  UnfittedAnalysis(name,geometry,vtu,results,linf_mass,linf_res,cond2_A,cond1_A)
end

function setup_results(problem::DarcyProblem,method::Unfitted,Ω,Γᵤ,uh,ph,γ,h,t)
    dΩ     = Measure(Ω,method.degree)
    dΓᵤ    = Measure(Γᵤ,method.degree)
    nᵤ     = get_normal_vector(Γᵤ)

    err_u    = problem.u - uh
    err_uN   = problem.u⋅nᵤ - uh⋅nᵤ
    err_p    = problem.p - ph
    err_divu = problem.g + (∇⋅(uh))

    results               = Dict{Symbol,Any}()
    results[:uh_L2]       = sqrt(sum(∫( uh⋅uh )dΩ))
    results[:ph_L2]       = sqrt(sum(∫( ph⋅ph )dΩ))
    results[:uh_Hdiv]     = sqrt(sum(∫( uh⋅uh + (∇⋅(uh))*(∇⋅(uh)) )dΩ))
    results[:divuh_L2]    = sqrt(sum(∫( (∇⋅(uh))*(∇⋅(uh)) )dΩ))

    results[:err_u_L2]    = sqrt(sum(∫( err_u⋅err_u )dΩ))
    results[:err_p_L2]    = sqrt(sum(∫( err_p⋅err_p )dΩ))
    results[:err_divu_L2] = sqrt(sum(∫( err_divu⋅err_divu )dΩ))
    results[:err_u_Hdiv]  = sqrt(sum(∫( err_u⋅err_u + (∇⋅(err_u))*(∇⋅(err_u)) )dΩ))

    if !isempty(problem.bnd_u)
        results[:uhN_L2]      = sqrt(sum(γ*∫( (uh⋅nᵤ)*(uh⋅nᵤ) )dΓᵤ))
        results[:err_uN_L2]   = sqrt(sum(γ*∫( err_uN⋅err_uN )dΓᵤ))
        results[:err_u_Vh]    = sqrt(sum(∫( err_u⋅err_u + (∇⋅(err_u))*(∇⋅(err_u)) )dΩ + γ*∫( err_uN⋅err_uN )dΓᵤ))
    else
        results[:uhN_L2]      = 0.0
        results[:err_uN_L2]   = 0.0
        results[:err_u_Vh]    = sqrt(sum(∫( err_u⋅err_u + (∇⋅(err_u))*(∇⋅(err_u)) )dΩ ))
    end
    
    results[:k] = method.k
    results[:h] = h
    results[:runtime] = time() - t
    results
end

function compute_linf_mass_err(model,dΩ,err,refFEₚ)
    
    Q       = FESpace(model,refFEₚ;conformity=:L2)    
    P       = TrialFESpace(Q)

    a(p,q) = ∫(q*p)dΩ
    l(q) = ∫(q*err)dΩ

    operator = AffineFEOperator(a,l,P,Q)
    projQh_err = solve(operator)

    norm(get_cell_dof_values(projQh_err),Inf)
end

function compute_linf_res(operator,u,p,X,xh)
    
    residual(x) = Gridap.FESpaces.residual(operator,x)
    linf_res    = norm(residual(xh),Inf)
    xh_ex       = interpolate_everywhere([u p],X)
    linf_res_ex = norm(residual(xh_ex),Inf)
    
    linf_res, linf_res_ex
end