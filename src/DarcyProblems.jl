#======================================================================# 
#                         Problem types                                #
#======================================================================#

abstract type ProblemType end

abstract type DarcyProblem <: ProblemType end

#----------------------------------------------------------------------#
#                         Darcy Problems                               #
#----------------------------------------------------------------------#
# Problem : DarcyOnCutSquare
"""
DarcyOnCutSquare()
...
# Arguments
- `name::String`: name of the problem.
- `u::Function`: the exact solution of the velocity field.
- `p::Function`: the exact solution of the pressure field.
- `f::Function`: the forcing for the momentum equation.
- `g::Function`: the forcing for the compressibility equation.
- `kinv::TensorValue` the inverse porosity matrix.
- `n::Integer`: the number of elements in each spatial direction.
- `pmin::Point{2,Float64}`: the bottom left point of the (uncut) interior of the domain of interest.
- `pmax::Point{2,Float64}`: the upper right point of the (uncut) interior of the domain of interest.
- `bnd_u::Vector{Int64}`: vector storing tags for velocity boundary condition [5=bottom,6=top,7=left,8=right]
- `bnd_p::Vector{Int64}`: vector storing tags for pressure boundary condition [5=bottom,6=top,7=left,8=right]
- `ε::Float64: cut length
...
"""
struct DarcyOnCutSquare <: DarcyProblem
    name::String
    u::Function
    p::Function
    f::Function
    g::Function
    kinv::TensorValue{2,2,Float64,4}
    n::Int
    pmin::Point{2,Float64}
    pmax::Point{2,Float64}
    bnd_u::Vector{Int64}
    bnd_p::Vector{Int64}
    ε::Float64
end

function SmoothSolenoidalDarcyOnCutSquare(;
    name     = SmoothSolenoidalDarcy2D_name*"OnCutSquare",
    u        = SmoothSolenoidalDarcy2D_u,
    p        = SmoothSolenoidalDarcy2D_p,
    f        = SmoothSolenoidalDarcy2D_f,
    g        = SmoothSolenoidalDarcy2D_g,
    kinv     = SmoothSolenoidalDarcy2D_kinv,
    n        = 10,
    pmin     = Point(0.0,0.0),
    pmax     = Point(1.0,1.0),
    bnd_u    = [7,8],
    bnd_p    = [5,6],
    ε        = 1.0e-1)
    DarcyOnCutSquare(name,u,p,f,g,kinv,n,pmin,pmax,bnd_u,bnd_p,ε)
end

function PresRobustDarcyOnCutSquare(;
    name     = PresRobustDarcy2D_name*"OnCutSquare",
    u        = PresRobustDarcy2D_u,
    p        = PresRobustDarcy2D_p,
    f        = PresRobustDarcy2D_f,
    g        = PresRobustDarcy2D_g,
    kinv     = PresRobustDarcy2D_kinv,
    n        = 10,
    pmin     = Point(0.0,0.0),
    pmax     = Point(1.0,1.0),
    bnd_u    = [7,8],
    bnd_p    = [5,6],
    ε        = 1.0e-1)
    DarcyOnCutSquare(name,u,p,f,g,kinv,n,pmin,pmax,bnd_u,bnd_p,ε)
end

function setup_geometry(problem::DarcyOnCutSquare,method::Unfitted)
    n     = problem.n
    pmin  = problem.pmin
    pmax  = problem.pmax
    ε     = problem.ε
    dp    = pmax - pmin     # dp[1] here is the interior (uncut) length of square
    @assert !is_simplex(method.topo)
    
    hint    = dp[1]/(n-2)     # element-size based on interior length
    bgpmin  = pmin - Point(hint,hint)
    bgpmax  = pmax + Point(hint,hint)
    bgdp    = bgpmax - bgpmin
    partition = (n,n)
    bgmodel = CartesianDiscreteModel(bgpmin,bgpmax,partition)
    h       = bgdp[1]/n
    @assert hint ≈ h

    e1 = VectorValue(1,0)
    e2 = VectorValue(0,1)
    x0 = pmin + Point(0.5*dp[1],0.5*dp[2])
    L1 = dp[1] + 2*ε
    L2 = dp[2] + 2*ε
    @assert (ε/h)<1.0

    plane1 = plane(x0=x0-0.5*L2*e2,v=-e2,name="bottom")
    plane2 = plane(x0=x0+0.5*L1*e1,v= e1,name="right")
    plane3 = plane(x0=x0+0.5*L2*e2,v= e2,name="top")
    plane4 = plane(x0=x0-0.5*L1*e1,v=-e1,name="left")
 
    geo12 = intersect(plane1,plane2)
    geo34 = intersect(plane3,plane4)
 
    square = intersect(geo12,geo34)
    cutgeo = cut(bgmodel, square)

    bnd_to_plane = Dict(5=>plane1,8=>plane2,6=>plane3,7=>plane4)

    function setup_embedded_boundary(bnd, bnd_to_plane, square, cutgeo)
        if !isempty(bnd)
            intersection = bnd_to_plane[bnd[1]]
            for b in bnd[2:end]
                intersection = intersect(bnd_to_plane[b],intersection)
            end
            Γ = EmbeddedBoundary(cutgeo,setdiff(square,intersection))
        else
            Γ = EmbeddedBoundary(cutgeo,bnd_to_plane[5], bnd_to_plane[6]) # this is empty!
        end
    end    

    Γᵤ = setup_embedded_boundary(problem.bnd_u, bnd_to_plane, square, cutgeo)
    Γₚ = setup_embedded_boundary(problem.bnd_p, bnd_to_plane, square, cutgeo)

    bgmodel, cutgeo, h, Γᵤ, Γₚ
end
