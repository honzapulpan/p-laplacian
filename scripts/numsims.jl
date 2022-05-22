using DrWatson
@quickactivate "Eigenvalue bounds of p-Laplacian"

include(srcdir("numerics.jl"))

# first batch
#allparams = Dict(
#    "p" => [2],
#    "mesh" => ["disk1.msh","disk2.msh", "disk3.msh"],
#    "order" => [1,],
#    "degree" => [1,2,],
#    "itnum" => [12,],
#    )

# second batch
#allparams = Dict(
#    "p" => [2],
#    "mesh" => ["disk3.msh", "disk4.msh", "disk5.msh", "disk6.msh"],
#    "order" => [1,2,3,], #[1,2,3],
#    "degree" => [2,3,], #[2,3],
#    "itnum" => [12,],
#    )

# third batch
allparams = Dict(
    "p" => [2],
    "mesh" => ["disk6.msh", "disk7.msh"],
    "order" => [3,], #[1,2,3],
    "degree" => [4,5,], #[2,3],
    "itnum" => [12,],
    )




dicts = dict_list(allparams)

function makesim(d::Dict)
    @unpack p, mesh, order, degree, itnum = d
    λ₁, runtime, λs = pl_num_eigenpair(p, mesh, order, degree, itnum)
    fulld = copy(d)
    fulld["λ₁"] = λ₁
    fulld["runtime"] = runtime
    fulld["λs"] = λs
    return fulld
end

for (i,d) in enumerate(dicts)
    f = makesim(d)
    @tagsave(datadir("simulations", savename(d, "jld2")), f)
end
# println("První vlastní číslo pro p=$(p): λ₁ ≈ $λ₁ (runtime=$runtime)")
