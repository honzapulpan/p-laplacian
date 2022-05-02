using DrWatson
@quickactivate "Eigenvalue bounds of p-Laplacian"

using Gridap, GridapGmsh, GridapMakie
using GLMakie
using FileIO
import Random
using LinearAlgebra: norm
using LineSearches: BackTracking
using Gridap.Geometry
using Gridap.ReferenceFEs
using Dates


function pl_num_eigenpair(p, mesh, order, degree)
    # start time measure
    startt = Dates.now()

    # load mesh from file
    model = GmshDiscreteModel(srcdir("mesh",mesh))

    # define test and trial spaces
    reffe = ReferenceFE(lagrangian,Float64,order)
    V0 = TestFESpace(model,reffe,conformity=:H1, dirichlet_tags="bnd") # dirichlet všude na hranici. Ta je v modelu označena bnd
    g = 0 # okrajové podmínky
    Ug = TrialFESpace(V0,g)

    # model triangulation 
    Ω = Triangulation(model)
    # Lebegues measure definition
    dΩ = Measure(Ω,degree)

    # residual and jacobian
    flux(∇u) = norm(∇u)^(p-2) * ∇u
    f(x) = 1
    res0(u,v) = ∫( ∇(v)⊙(flux∘∇(u)) - v*f )*dΩ

    dflux(∇du,∇u) = (p-2)*norm(∇u)^(p-4)*(∇u⊙∇du)*∇u+norm(∇u)^(p-2)*∇du
    jac(u,du,v) = ∫( ∇(v)⊙(dflux∘(∇(du),∇(u))) )*dΩ

    op = FEOperator(res0,jac,Ug,V0)

    # solver definition
    nls = NLSolver(
      show_trace=false, method=:newton, linesearch=BackTracking()) # show_trace --> vypíše průběh řešiče
    solver = FESolver(nls)

    # first solution
    Random.seed!(1234)
    x = rand(Float64,num_free_dofs(Ug))
    uh0 = FEFunction(Ug,x)
    uh, = solve!(uh0,solver,op)

    # algorithm
    q = p-.01
    pts = get_node_coordinates(Ω) # potřebujeme k evaluaci řešení v uzlových bodech

    for m in 0:8
        uhsup = maximum(abs.(lazy_map(uh,pts)))
        raiseprw(g1,g2) = g1*(g2/uhsup)^(q-1)
        res(u,v) = ∫( ∇(v)⊙(flux∘∇(u)) - raiseprw∘(v,uh) )*dΩ
        op = FEOperator(res,jac,Ug,V0)

        Random.seed!(1234)
        x = rand(Float64,num_free_dofs(Ug))
        uh0 = FEFunction(Ug,x)
        uh, = solve!(uh0,solver,op);

    end

    #uhsup = maximum(abs.(lazy_map(uh,pts)))
    λ₁ = (1/(maximum(abs.(lazy_map(uh,pts))))^(p-1))

    stopt = Dates.now()
    runtime = Dates.canonicalize(Dates.CompoundPeriod(stopt-startt))
    return λ₁, runtime
end
