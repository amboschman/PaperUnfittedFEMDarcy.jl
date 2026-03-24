#======================================================================#
#                               Drivers                                #
#======================================================================#

# Driver: Bulk Ghost Penalty stabilized unfitted Darcy Problem
"""
main(problem::DarcyProblem,method::Unfitted,stabilization::BulkGhostPenalty,analysis::UnfittedAnalysis)
...
# Arguments
- `problem::DarcyProblem`: Darcy problem with p, u or mixed boundary conditions.
- `method::Unfitted`: Unfitted method.
- `stabilization::BulkGhostPenaltyStabilization`: BulkGhostPenalty stabilization.
- `analysis::UnfittedAnalysis`: unfitted analysis method
...
"""
function main(problem::DarcyProblem,method::Unfitted,stabilization::BulkGhostPenaltyStabilization,analysis::UnfittedAnalysis)

    t = time()

    # Problem parameters
    u    = problem.u
    p    = problem.p
    f    = problem.f
    g    = problem.g
    kinv = problem.kinv

    # Boundary conditions
    u_Γᵤ = u
    p_Γₚ = p

    # Geometry
    @assert !is_simplex(method.topo)
    bgmodel, cutgeo, h, Γᵤ, Γₚ = setup_geometry(problem, method)

    # Triangulations
    Ωbg  = Triangulation(bgmodel)
    Ωact = Triangulation(cutgeo,ACTIVE)

    # Physical domain
    Ω = Triangulation(cutgeo,PHYSICAL)

    # Numerical integration
    degree=method.degree
    dΩ    = Measure(Ω,degree)
    dΓᵤ   = Measure(Γᵤ,degree)
    dΓₚ   = Measure(Γₚ,degree)

    # Retrieve stabilization + AL parameters
    ξ    = stabilization.ξ              # If non-zero, create aggregates and add (non-zero) ξᵤ-, ξₚ-, ξₛ- and ξₘ-weighted stabilization terms
    ξ_AL = stabilization.ξ_AL           # If non-zero, add ξ_AL-weighted Augmented Lagrangian stabilization terms to the weak form

    # Setup aggregates
    if !isapprox(ξ,0.0)
        # Setup aggregates
        strategy  = AggregateAllCutCells()
        aggregates= aggregate(strategy,cutgeo)
        aggregate_to_cells=setup_aggregate_to_cells(aggregates)
        aggregates_bounding_box_model=setup_aggregates_bounding_box_model(bgmodel,aggregate_to_cells)
        # Setup agg cells and triangulation
        agg_cells    =flatten(aggregate_to_cells) 
        @assert length(agg_cells)>0 "No aggregated cells"
        Ωbg_agg_cells=view(Ωbg,agg_cells)
        # Numerical integration on agg cells
        dΩbg_agg_cells = Measure(Ωbg_agg_cells,degree)
    end
   
    # Export geometry
    if analysis.geometry
        export_geometry("trian-act-$(problem.name)-$(method.name)-n$(problem.n)", Ωact)
        export_geometry("trian-phys-$(problem.name)-$(method.name)-n$(problem.n)", Ω)
        export_geometry("trian-bg-$(problem.name)-$(method.name)-n$(problem.n)", Ωbg)
        if !isapprox(ξ,0.0)
        writevtk(Ωbg,"trian-bg-incl-agg-$(problem.name)-$(method.name)-n$(problem.n)",celldata=["aggregate"=>aggregates,"color_agg"=>color_aggregates(aggregates,bgmodel)])
        end
        writevtk(Γₚ,"trian-Gp-$(problem.name)-$(method.name)-n$(problem.n)",cellfields=["n"=>get_normal_vector(Γₚ)])
        writevtk(Γᵤ,"trian-Gu-$(problem.name)-$(method.name)-n$(problem.n)",cellfields=["n"=>get_normal_vector(Γᵤ)])
    end

    # FE Spaces
    k       = method.k
    contra_variant_piola_map_type = method.contra_variant_piola_map_type
    refFEᵤ = ReferenceFE(method.shape_u(),Float64,k)
    if isa(method.shape_u(),BDM)
        @assert k>0
        refFEₚ = ReferenceFE(method.shape_p(),Float64,k-1)
    elseif isa(method.shape_u(),RaviartThomas)
        refFEₚ = ReferenceFE(method.shape_p(),Float64,k)
    end
    V       = FESpace(Ωact,refFEᵤ;conformity=:Hdiv,contra_variant_piola_map_type=contra_variant_piola_map_type)
    if isempty(problem.bnd_p)
        # Set up zero-mean pressure space, with fixed interior dof
        Qnzm = FESpace(Ωact, ReferenceFE(lagrangian,Float64,k), conformity=:L2)
        int_cells = restrict_cells(cutgeo,IN) # global (background mesh') cell identifiers of interior cells
        if !isapprox(ξ,0.0)
            nonagg_int_cell = int_cells[findfirst(!in(agg_cells),int_cells)]   # cell identifier (background mesh) of first interior cell not in aggregate
            @assert nonagg_int_cell!==Nothing "No interior cells that are not part of an aggregate"
            local_id = findfirst(isequal(nonagg_int_cell),Ωact.tface_to_mface) # local (active mesh') cell id of interior cell not in aggregate
        else 
            local_id = findfirst(isequal(int_cells[1]),Ωact.tface_to_mface) # local (active mesh') cell id of interior cell
        end
        dof_to_fix = get_cell_dof_ids(Qnzm)[local_id][1]
        spaceWithConstantFixed = Gridap.FESpaces.FESpaceWithConstantFixed(Qnzm,true,Int64(dof_to_fix))
        Qzm_vol_i = assemble_vector(v->∫(v)*dΩ,Qnzm)
        Qzm_vol = sum(Qzm_vol_i)
        Q     = Gridap.FESpaces.ZeroMeanFESpace(spaceWithConstantFixed,Qzm_vol_i,Qzm_vol)
    else
        Q     = FESpace(Ωact,refFEₚ;conformity=:L2)
    end
    U = TrialFESpace(V)
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U, P])
    Y = MultiFieldFESpace([V, Q])
    dx    = get_trial_fe_basis(X)
    dy    = get_fe_basis(Y)
    du,dp = dx
    dv,dq = dy

    if !isapprox(ξ,0.0)
        # ref_agg_cell_to_agg_cell_map: \hat{K} -> K
        ref_agg_cell_to_agg_cell_map=get_cell_map(Ωbg_agg_cells)
        agg_cells_to_aggregate      =setup_cells_to_aggregate(aggregate_to_cells)
        ref_agg_cell_to_ref_bb_map  =setup_ref_agg_cell_to_ref_bb_map(aggregates_bounding_box_model,
                                                                    agg_cells_to_aggregate,ref_agg_cell_to_agg_cell_map)
        # Spaces on bounding boxes
        Qbb=FESpace(aggregates_bounding_box_model,refFEₚ,conformity=:L2) # We need a DG space to represent the L2 projection
        Pbb=TrialFESpace(Qbb)
        pbb=get_trial_fe_basis(Pbb)
        qbb=get_fe_basis(Qbb)
        Vbb=FESpace(aggregates_bounding_box_model,refFEᵤ,conformity=:L2,contra_variant_piola_map_type=contra_variant_piola_map_type)
        Ubb=TrialFESpace(Vbb)
        ubb=get_trial_fe_basis(Ubb)
        vbb=get_fe_basis(Vbb)
    end

    # Weak form
    m     = method.m
    γ₀    = method.γ₀
    s     = method.s
    γ     = γ₀*h.^(s)
    nᵤ = get_normal_vector(Γᵤ)
    nₚ = get_normal_vector(Γₚ)
    # We rely here on the fact that SparseArrays will ignore zero entries when assembling the global matrix, so no need to split this into different cases.
    a((u,p), (v,q)) = ∫(v⋅kinv⋅u - (∇⋅v)*p - (∇⋅u)*q)dΩ + ∫(γ*(u⋅nᵤ)*(v⋅nᵤ))dΓᵤ + ∫((v⋅nᵤ)*p)dΓᵤ + ∫(ξ_AL*(∇⋅u)*(∇⋅v))dΩ + ∫(m*(u⋅nᵤ)*q)dΓᵤ
    l((v,q))        = ∫(v⋅f + q*g)dΩ - ∫((v⋅nₚ)*p_Γₚ)dΓₚ + ∫(γ*(v⋅nᵤ)*(u_Γᵤ⋅nᵤ))dΓᵤ - ξ_AL*∫((g)*(∇⋅v))dΩ + ∫(m*q*(u_Γᵤ⋅nᵤ))dΓᵤ
    wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
    vec_wr=Gridap.FESpaces.collect_cell_vector(Y,l(dy))

    # Add stabilization terms
    if !isapprox(ξ,0.0)
        
        # Retrieve additional stabilization parameters
        ξᵤ = ξ*stabilization.ξᵤ
        ξₚ = ξ*stabilization.ξₚ
        ξₛ = ξ*stabilization.ξₛ 
        ξₘ = ξ*stabilization.ξₘ 

        # Selecting relevant global dofs ids of aggregate cells (from background mesh)
        Ωbg_agg_cell_dof_ids   = get_cell_dof_ids(X,Ωbg_agg_cells)
        U_Ωbg_agg_cell_dof_ids = _restrict_to_block(Ωbg_agg_cell_dof_ids, 1) 
        P_Ωbg_agg_cell_dof_ids = _restrict_to_block(Ωbg_agg_cell_dof_ids, 2)

        # Computing local (per aggregate) dof ids 
        aggregate_to_local_cells =setup_aggregate_to_local_cells(aggregate_to_cells)
        U_agg_cells_local_dof_ids=compute_agg_cells_local_dof_ids(U_Ωbg_agg_cell_dof_ids, aggregate_to_local_cells)
        P_agg_cells_local_dof_ids=compute_agg_cells_local_dof_ids(P_Ωbg_agg_cell_dof_ids, aggregate_to_local_cells)

        # Compute global dofs ids per aggregate and reindex these 
        U_aggregate_dof_ids=compute_aggregate_dof_ids(U_Ωbg_agg_cell_dof_ids,aggregate_to_cells)
        U_agg_cells_to_aggregate_dof_ids=lazy_map(Reindex(U_aggregate_dof_ids),agg_cells_to_aggregate)
        P_aggregate_dof_ids=compute_aggregate_dof_ids(P_Ωbg_agg_cell_dof_ids,aggregate_to_cells)
        P_agg_cells_to_aggregate_dof_ids=lazy_map(Reindex(P_aggregate_dof_ids),agg_cells_to_aggregate)

        # Setup cut cells and triangulation
        aggregate_to_cut_cells = restrict_aggregate_to_cells(cutgeo,aggregate_to_cells,GridapEmbedded.Interfaces.CUT)
        cut_cells = flatten(aggregate_to_cut_cells)
        Ωbg_cut_cells   = view(Ωbg,cut_cells)
        dΩbg_cut_cells  = Measure(Ωbg_cut_cells,degree)

        # Selecting relevant global dofs ids of cut cells (from background mesh)
        Ωbg_cut_cell_dof_ids   = get_cell_dof_ids(X,Ωbg_cut_cells)
        U_Ωbg_cut_cell_dof_ids = _restrict_to_block(Ωbg_cut_cell_dof_ids, 1) 
        P_Ωbg_cut_cell_dof_ids = _restrict_to_block(Ωbg_cut_cell_dof_ids, 2)

        # Compute global dofs ids per aggregate and reindex these 
        cut_cells_to_aggregate = setup_cells_to_aggregate(aggregate_to_cut_cells)
        U_cut_cells_to_aggregate_dof_ids=lazy_map(Reindex(U_aggregate_dof_ids),cut_cells_to_aggregate)
        P_cut_cells_to_aggregate_dof_ids=lazy_map(Reindex(P_aggregate_dof_ids),cut_cells_to_aggregate)

        # Setup projections
        du_proj_Vbb, dv_proj_Vbb = setup_L2_proj_in_bb_space(
            dΩbg_agg_cells,             # measure of aggregated cells in background domain
            ref_agg_cell_to_ref_bb_map, # map
            agg_cells_to_aggregate,     # 
            aggregate_to_local_cells,   # 
            du,                         # Trial basis (to project) 
            dv,                         # Test basis 
            ubb,                        # Trial basis of bounding box space Vbb
            vbb,                        # Test basis of bounding box space Vbb
            identity,                   # operation to be applied to u and v
            U_agg_cells_local_dof_ids)  # aggregates local dof ids for space U   

        dp_proj_Qbb, dq_proj_Qbb = setup_L2_proj_in_bb_space(
            dΩbg_agg_cells,             # measure of aggregated cells in background domain
            ref_agg_cell_to_ref_bb_map, # map
            agg_cells_to_aggregate,     # 
            aggregate_to_local_cells,   # 
            dp,                         # Trial basis (to project) 
            dq,                         # Test basis 
            pbb,                        # Trial basis of bounding box space Qbb
            qbb,                        # Test basis of bounding box space Qbb
            identity,                   # operation to be applied to u and v
            P_agg_cells_local_dof_ids)  # aggregates local dof ids for space P   

        div_du_proj_Qbb, div_dv_proj_Qbb = setup_L2_proj_in_bb_space(
            dΩbg_agg_cells,             # measure of aggregated cells in background domain
            ref_agg_cell_to_ref_bb_map, # map
            agg_cells_to_aggregate,     # 
            aggregate_to_local_cells,   # 
            du,                         # Trial basis (to project) 
            dv,                         # Test basis 
            pbb,                        # Trial basis of bounding box space Vbb
            qbb,                        # Test basis of bounding box space Vbb
            divergence,                 # operation to be applied to u and v
            U_agg_cells_local_dof_ids)  # aggregates local dof ids for space U   

        if !isapprox(ξₘ,0.0) 
            g_proj_Qbb = setup_L2_proj_in_bb_space(dΩbg_agg_cells,
                ref_agg_cell_to_ref_bb_map, 
                agg_cells_to_aggregate,      
                aggregate_to_local_cells,    
                g,                   
                pbb,                        
                qbb) 
        end

        # Apply BlockMap
        if (PaperUnfittedFEMDarcy.GridapEmbedded._is_multifield_fe_basis_component(dv))
            nfields=PaperUnfittedFEMDarcy.GridapEmbedded._nfields(dv)
            fieldid=PaperUnfittedFEMDarcy.GridapEmbedded._fieldid(dv)
            U_Ωbg_cut_cell_dof_ids=lazy_map(PaperUnfittedFEMDarcy.Gridap.Fields.BlockMap(nfields,fieldid),U_Ωbg_cut_cell_dof_ids)
        end

        if (PaperUnfittedFEMDarcy.GridapEmbedded._is_multifield_fe_basis_component(du))
            nfields=PaperUnfittedFEMDarcy.GridapEmbedded._nfields(du)
            fieldid=PaperUnfittedFEMDarcy.GridapEmbedded._fieldid(du)
            U_cut_cells_to_aggregate_dof_ids=
            lazy_map(PaperUnfittedFEMDarcy.Gridap.Fields.BlockMap(nfields,fieldid),U_cut_cells_to_aggregate_dof_ids)
        end

        if (PaperUnfittedFEMDarcy.GridapEmbedded._is_multifield_fe_basis_component(dq))
            nfields=PaperUnfittedFEMDarcy.GridapEmbedded._nfields(dq)
            fieldid=PaperUnfittedFEMDarcy.GridapEmbedded._fieldid(dq)
            P_Ωbg_cut_cell_dof_ids=lazy_map(PaperUnfittedFEMDarcy.Gridap.Fields.BlockMap(nfields,fieldid),P_Ωbg_cut_cell_dof_ids)
        end

        if (PaperUnfittedFEMDarcy.GridapEmbedded._is_multifield_fe_basis_component(dp))
            nfields=PaperUnfittedFEMDarcy.GridapEmbedded._nfields(dp)
            fieldid=PaperUnfittedFEMDarcy.GridapEmbedded._fieldid(dp)
            P_cut_cells_to_aggregate_dof_ids=
            lazy_map(PaperUnfittedFEMDarcy.Gridap.Fields.BlockMap(nfields,fieldid),P_cut_cells_to_aggregate_dof_ids)
        end

        # Determine domain of intergration for stabilization terms
        dD = dΩbg_cut_cells

        # Compute stabilization terms for u
        if !isapprox(ξᵤ,0.0)
            if stabilization.ζ # "U_DIFF"
                wu,ru,cu=bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dD,
                            Ωbg_agg_cells,
                            ξᵤ,
                            du,
                            dv,
                            du_proj_Vbb,
                            dv_proj_Vbb,
                            U_Ωbg_cut_cell_dof_ids,
                            U_Ωbg_cut_cell_dof_ids,
                            U_cut_cells_to_aggregate_dof_ids,
                            U_cut_cells_to_aggregate_dof_ids,
                            identity,
                            identity)
            else # "U_FULL"
                wu,ru,cu=bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dD,
                            Ωbg_agg_cells,
                            ξᵤ,
                            du,
                            dv,
                            du_proj_Vbb,
                            U_Ωbg_cut_cell_dof_ids,
                            U_Ωbg_cut_cell_dof_ids,
                            U_cut_cells_to_aggregate_dof_ids,
                            identity,
                            identity)
            end
            push!(wrc[1], wu...)
            push!(wrc[2], ru...)
            push!(wrc[3], cu...)
        end

        # Compute stabilization terms for p
        if !isapprox(ξₚ,0.0)
            if stabilization.ζ # "P_DIFF"
                wp,rp,cp=bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dD,
                            Ωbg_agg_cells,
                            ξₚ,
                            dp,
                            dq,
                            dp_proj_Qbb,
                            dq_proj_Qbb,
                            P_Ωbg_cut_cell_dof_ids,
                            P_Ωbg_cut_cell_dof_ids,
                            P_cut_cells_to_aggregate_dof_ids,
                            P_cut_cells_to_aggregate_dof_ids,
                            identity,
                            identity)
            else # "P_FULL"
                wp,rp,cp=bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dD,
                            Ωbg_agg_cells,
                            ξₚ,
                            dp,
                            dq,
                            dp_proj_Qbb,
                            P_Ωbg_cut_cell_dof_ids,
                            P_Ωbg_cut_cell_dof_ids,
                            P_cut_cells_to_aggregate_dof_ids,
                            identity,
                            identity)
            end
            push!(wrc[1], wp...)
            push!(wrc[2], rp...)
            push!(wrc[3], cp...)
        end

        # Compute stabilization terms for divu
        if !isapprox(ξₛ,0.0)
            if stabilization.ζ # "DIVU_DIFF"
                wdiv, rdiv, cdiv = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dD,
                                    Ωbg_agg_cells,
                                    ξₛ,
                                    du,
                                    dv,
                                    div_du_proj_Qbb,
                                    div_dv_proj_Qbb,
                                    U_Ωbg_cut_cell_dof_ids,
                                    U_Ωbg_cut_cell_dof_ids,
                                    U_cut_cells_to_aggregate_dof_ids,
                                    U_cut_cells_to_aggregate_dof_ids,
                                    divergence,
                                    divergence)
            else # "DIVU_FULL"
                wdiv, rdiv, cdiv = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dD,
                                    Ωbg_agg_cells,
                                    ξₛ,
                                    du,
                                    dv,
                                    div_du_proj_Qbb,
                                    U_Ωbg_cut_cell_dof_ids,
                                    U_Ωbg_cut_cell_dof_ids,
                                    U_cut_cells_to_aggregate_dof_ids,
                                    divergence,
                                    divergence)
            end
            push!(wrc[1], wdiv...)
            push!(wrc[2], rdiv...)
            push!(wrc[3], cdiv...)
        end

        # Compute mixed stabilization terms
        if !isapprox(ξₘ,0.0)
            if stabilization.ζ 
                # "DIVUQ_DIFF"
                wdivuq, rdivuq, cdivuq = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dD,
                    Ωbg_agg_cells,
                    -ξₘ,
                    du,
                    dq,
                    div_du_proj_Qbb,
                    dq_proj_Qbb,
                    U_Ωbg_cut_cell_dof_ids,
                    P_Ωbg_cut_cell_dof_ids,
                    U_cut_cells_to_aggregate_dof_ids,
                    P_cut_cells_to_aggregate_dof_ids,
                    divergence,
                    identity)
                # "DIVP_DIFF"
                wpdivv, rpdivv, cpdivv = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dD,
                    Ωbg_agg_cells,
                    -ξₘ,
                    dp,
                    dv,
                    dp_proj_Qbb,
                    div_dv_proj_Qbb,
                    P_Ωbg_cut_cell_dof_ids,
                    U_Ωbg_cut_cell_dof_ids,
                    P_cut_cells_to_aggregate_dof_ids,
                    U_cut_cells_to_aggregate_dof_ids,
                    identity,
                    divergence)
                # "GQ_DIFF"
                wgq, rgq = bulk_ghost_penalty_stabilization_collect_cell_vector_on_D(dD,
                    ξₘ,  
                    dq,
                    dq_proj_Qbb,
                    P_Ωbg_cut_cell_dof_ids,
                    P_cut_cells_to_aggregate_dof_ids,
                    identity,
                    g,
                    g_proj_Qbb)
            else
                # "DIVUQ_FULL"
                wdivuq, rdivuq, cdivuq = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dD,
                    Ωbg_agg_cells,
                    -ξₘ,
                    du,
                    dq,
                    div_du_proj_Qbb,
                    U_Ωbg_cut_cell_dof_ids,
                    P_Ωbg_cut_cell_dof_ids,
                    U_cut_cells_to_aggregate_dof_ids,
                    divergence,
                    identity)
                # "DIVP_FULL"
                wpdivv, rpdivv, cpdivv = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dD,
                    Ωbg_agg_cells,
                    -ξₘ,
                    dp,
                    dv,
                    dp_proj_Qbb,
                    P_Ωbg_cut_cell_dof_ids,
                    U_Ωbg_cut_cell_dof_ids,
                    P_cut_cells_to_aggregate_dof_ids,
                    identity,
                    divergence)
                # "GQ_FULL"
                wgq, rgq = bulk_ghost_penalty_stabilization_collect_cell_vector_on_D(dD,
                    ξₘ, 
                    dq, 
                    P_Ωbg_cut_cell_dof_ids,
                    identity,
                    g,
                    g_proj_Qbb)
            end
            push!(wrc[1], wdivuq...)
            push!(wrc[2], rdivuq...)
            push!(wrc[3], cdivuq...)
            push!(wrc[1], wpdivv...)
            push!(wrc[2], rpdivv...)
            push!(wrc[3], cpdivv...)
            push!(vec_wr[1],wgq...)
            push!(vec_wr[2],rgq...)
        end
    end

    # Assembly
    assem   = SparseMatrixAssembler(X,Y)
    A       = assemble_matrix(assem, wrc)
    b       = assemble_vector(assem, vec_wr)
    @assert num_free_dofs(X) == size(A,2) "Incompatible trial space and matrix"
    @assert num_free_dofs(Y) == size(A,1) "Incompatible test space and matrix"
    op = PaperUnfittedFEMDarcy.Gridap.Algebra.AffineOperator(A,b)
    operator = AffineFEOperator(X,Y,op)
    if method.PETSc
        options = """
            -ksp_type preonly -ksp_error_if_not_converged true
            -pc_type lu -pc_factor_mat_solver_type mumps
            -mat_mumps_icntl_1 4
            -mat_mumps_icntl_4 0 
            -mat_mumps_icntl_7 0 
            -mat_mumps_icntl_14 100
            -mat_mumps_icntl_28 1
            -mat_mumps_icntl_29 2
            -mat_mumps_cntl_3 1.0e-6
            """
        GridapPETSc.Init(args=split(options))
        ls = GridapPETSc.PETScLinearSolver()
        solver = LinearFESolver(ls)
        xh = solve(solver,operator)
    else
        xh = solve(operator)
    end
    uh, ph = xh

    # Zero-mean pressure correction
    if isempty(problem.bnd_p)
        area  = sum(∫(1.0)dΩ)
        mean_p = sum(∫(p)dΩ)/area # mean presure exact sol
        ph  = ph + mean_p
    end

    # Export results as vtu file
    if analysis.vtu
        writevtk(Ω,"$(problem.name)-$(method.name)-n$(problem.n)-k$(method.k)",cellfields=["uh"=>uh,"ph"=>ph,"divuh"=>(∇⋅uh), "u-uh"=>(u-uh),"p-ph"=>(p-ph),"g+divuh"=>(g+(∇⋅uh))])
    end

    # Compute (error) norms
    if analysis.results
        results = setup_results(problem,method,Ω,Γᵤ,uh,ph,γ,h,t)
    else
        results               = Dict{Symbol,Any}()
    end

    # Compute linf of mass conservation error
    #TODO: check this
    if analysis.linf_mass
        results[:linf_mass] = compute_linf_mass_err(Ωact,dΩ,∇⋅(uh)+g,refFEₚ)
    end

    # Compute linf of the residual
    if analysis.linf_res
        if isempty(problem.bnd_p)
            p_zeromean(x) = p(x) - mean_p
            results[:linf_res],results[:linf_res_ex] = compute_linf_res(operator,u,p_zeromean,X,xh)
        else
            results[:linf_res],results[:linf_res_ex] = compute_linf_res(operator,u,p,X,xh)
        end
    end

    # Compute condition nr
    if analysis.cond2_A || analysis.cond1_A
        array_A = Array(A)
        if analysis.cond2_A
            results[:cond2_A] = cond(array_A,2)
        else 
            results[:cond2_A] = NaN
        end
        if analysis.cond1_A
            results[:cond1_A] = cond(array_A,1)
        else 
            results[:cond1_A] = NaN
        end
    end
    results
end