using DrWatson
@quickactivate "Eigenvalue bounds of p-Laplacian"


using DifferentialEquations
using IntervalArithmetic


"""
    plaplace_solve(λi, p, n; u₂0=1.0, dom=(0.0, 1.0))

Numericaly solves p-Laplace equation using shooting method.

Arguments:
λᵢₙᵢₜ ... initial interval for λ₁
p    ... p of p-Laplacian
n    ... number of solution points
u₂0  ... initial condition for u₂
dom  ... domain

Return:
t    ... division points where the solution is interpolated,
            double precision value and its interval representation
tᴵ   ... homogenous mesh of domain interval
U₁, U₁ᴵ  ... numerical approximation of u₁ at t, 
            double precision value and its interval representation
U₂, U₂ᴵ  ... numerical approximation of u₂ at t,
            double precision value and its interval representation
Λ₁       ... numerical approximation of first eigenvalue λ₁
"""
function plaplace_solve(λᵢₙᵢₜ, p, n; u₂0=1.0, dom=(0.0, 1.0))
   
    function sl(du,u,P,t)
        λ, p = P
        du[1] = abs(u[2])^(1/(p-1)) * sign(u[2])
        du[2] = -λ * abs(u[1])^(p-1)*sign(u[1]) 
    end
    
    tl, tr = dom
    
    u0 = [0.0; u₂0;] # initial condition
    a, b = λᵢₙᵢₜ
    Λ₁ = (a + b)/2
    Δt = (tr-tl)/(n-1)
    e = 1e-12 # stop condition

    while (b-a) >= e
        prob = ODEProblem(sl, u0, dom, (Λ₁, p))
        sol = solve(prob, saveat=Δt, abstol=1e-8,reltol=1e-8)
        if sol(tr)[1] == 0
            break
        else
            probA = ODEProblem(sl, u0, dom, (a, p))
            solA = solve(probA, saveat=Δt, abstol=1e-8,reltol=1e-8)
            probS = ODEProblem(sl, u0, dom, (Λ₁, p))
            solS = solve(probS, saveat=Δt, abstol=1e-8,reltol=1e-8)
            if solA(tr)[1] * solS(tr)[1] < 0
                b = Λ₁
            else
                a = Λ₁
            end
            Λ₁ = (a+b)/2
        end
    end

    prob = ODEProblem(sl, u0, dom, (Λ₁, p))
    sol = solve(prob, saveat=Δt, abstol=1e-8,reltol=1e-8)
    
    t = collect(LinRange(0,1,n-1))
    #tᴵ = [@interval(i) for i in t]
    tᴵ = []
    for i in 1:length(t)-1
        append!(tᴵ, t[i]..t[i+1])  
    end
    
    
    U₁ = [u[1] for u in sol(t).u]
    U₁ᴵ = [@interval(u[1]) for u in sol(t).u]
    U₂ = [u[2] for u in sol(t).u]
    U₂ᴵ = [@interval(u[2]) for u in sol(t).u]

    # první derivace použitá pro polynomy
    Ud = []
    for u2 in U₂
        if u2 >= 0
           append!(Ud, u2^(1/(p-1))) 
        else
           append!(Ud, -(-u2)^(1/(p-1)))  
        end
    end
    
    Udᴵ = []
    for u2I in U₂ᴵ
        if u2I >= 0
           append!(Udᴵ, u2I^(1/(p-1))) 
        else
           append!(Udᴵ, -(-u2I)^(1/(p-1)))  
        end
    end

    return t, tᴵ, U₁, U₁ᴵ, U₂, U₂ᴵ, Ud, Udᴵ, Λ₁
end;