using DrWatson
@quickactivate "Eigenvalue bounds of p-Laplacian"

include(srcdir("numerics.jl"))

allparams = Dict(
    "p" => [2],
    "mesh" => ["disk7.msh"],#"disk3.msh", "disk4.msh", "disk5.msh", "disk6.msh"],
    "order" => [1,2,3], #,4,5],
    "degree" => [2,3],
    )
#p = 2
#mesh = "disk4.msh"
#order = 1
#degree = 2

dicts = dict_list(allparams)


function makesim(d::Dict)
    @unpack p, mesh, order, degree = d
    λ₁, runtime = pl_num_eigenpair(p, mesh, order, degree)
    fulld = copy(d)
    fulld["λ₁"] = λ₁
    fulld["runtime"] = runtime
    return fulld
end

for (i,d) in enumerate(dicts)
    f = makesim(d)
    @tagsave(datadir("simulations", savename(d, "jld2")), f)
end
# println("První vlastní číslo pro p=$(p): λ₁ ≈ $λ₁ (runtime=$runtime)")
