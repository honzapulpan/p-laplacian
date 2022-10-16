using DrWatson
@quickactivate "Eigenvalue bounds of p-Laplacian"


using IntervalArithmetic


"""
    get_v1(p, V₂, t)

Rebuilds interval expression of `V₁` by integrating `V₂`.

Arguments:
p             ... p of p-Laplacian
V₂            ... intervals values of v₂
t             ... division points
"""
function get_v1(p, V₂, t)
    
    f(x) = abs(x)^(1/(p-1))*sign(x)
    ni = mince(0..1,length(V₂))

    V₁_tmp = Interval[0..0]
    for i in 1:length(V₂)
        append!(V₁_tmp, V₁_tmp[end] + f(V₂[i]) * diam(ni[i]))
    end
    
    V₁ = Interval[]
    for i in 1:length(V₁_tmp)-1
        append!(V₁, V₁_tmp[i] ∪ V₁_tmp[i+1])
    end
    
    V₁ = V₁ .- inf(minimum(V₁))
    
    return V₁
end;


"""
    der_cubic_spline(csc, t; ns=1)

Computes the first derivative `V` of a given interval cubic spline. 

Arguments:
csc           ... cubic spline coefficients interval representation
t             ... division points
ns            ... number of division points for each single piece of spline
"""
function der_cubic_spline(csc, t; ns=1)
    
    V_tmp = Interval[]
    csc_Vder = [ [@interval(3) * c[1], @interval(2) * c[2], c[3]] for c in csc ] 
    for i in 1:length(t)-1
        x_dom = t[i]..t[i+1] 
        x_int = mince(x_dom,ns)
        f(x) = csc_Vder[i][3] + (x-t[i])*(csc_Vder[i][2] + (x-t[i])*csc_Vder[i][1])
        append!(V_tmp, f.(x_int))
    end
    
    V = Interval[]
    for i in 1:length(V_tmp)-1
        append!(V, V_tmp[i] ∪ V_tmp[i+1])
    end
    append!(V, V_tmp[end])
        
    return V
end;



"""
    lower_estimate(V₂_der, V₁, p)

Returns guaranteed lower estimate of first eigenvalue `λ₁ˡᵒʷ` and interval values of -u'₂/u₁⁽ᵖ⁻¹⁾ in `Fˡᵒʷ`.

Arguments:
V₂_der  ... interval values of v₂'
V₁      ... interval values of v₁
"""
function lower_estimate(V₂_der, V₁, p)
    f(x,y) = -x / y^(p-1)
    λ₁_tmp = f.(V₂_der, V₁)

    Fˡᵒʷ = Interval[]
    for i in 1:length(λ₁_tmp)-1
        append!(Fˡᵒʷ, λ₁_tmp[i] ∪ λ₁_tmp[i+1])
    end
    append!(Fˡᵒʷ,λ₁_tmp[end])
    
    λ₁ˡᵒʷ = inf(minimum(Fˡᵒʷ))
    
    return λ₁ˡᵒʷ, Fˡᵒʷ
end;


"""
    upper_estimate(V₁, V₁_der, p)

Returns guaranteed upper estimate of the first eigenvalue `λ₁ᵘᵖ`.

Arguments:
p             ... p of p-Laplacian
V₁           ... interval values of v₁
V₁_der       ... interval values of v₁'
"""
function upper_estimate(V₁, V₁_der, p)
    
    f(x) = abs(x)^(p)
    ni = mince(0..1,length(V₁))

    numerator = 0..0
    for i in 1:length(V₁_der)
        numerator = numerator + f(V₁_der[i]) * diam(ni[i])
    end
    
    denominator = 0..0
    for i in 1:length(V₁)
        denominator = denominator + f(V₁[i]) * diam(ni[i])
    end

    λ₁ᵘᵖ = sup(numerator/denominator)

    return λ₁ᵘᵖ 
end;