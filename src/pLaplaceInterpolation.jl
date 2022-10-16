using DrWatson
@quickactivate "Eigenvalue bounds of p-Laplacian"


#using DifferentialEquations
using IntervalArithmetic
using LinearAlgebra

"""
   linear_interpolation(t, U₁ᴵ)

IA interpolation of given points by 3rd degree polynomials . 
Returns coefficients `csc_V` as well as interval values `V` 
of the polynomials.

Arguments:
t            ... division points
Uᴵ           ... interval values to be interpolated
"""

function linear_interpolation(t, U₁ᴵ)
   
    csc_V = Vector[]
 
    for i in 1:length(U₁ᴵ)-1
        a = b = c = d = 0
        c = (U₁ᴵ[i]-U₁ᴵ[i+1]) / (t[1]-t[2]) 
        d = (t[1]*U₁ᴵ[i+1]-t[2]*U₁ᴵ[i]) / (t[1]-t[2])
        
        append!(csc_V, [Interval{Float64}[a,b,c,d]])
    end
    
    # napočítat boxy výsledné fce
    V = Interval{Float64}[] 
    for i in 1:length(U₁ᴵ)-1
        x_int = t[i]..t[i+1]
        f(x) = csc_V[i][4] + (x-t[i])*csc_V[i][3]
        append!(V, f.(x_int))
    end
        
    return csc_V, V
end;


"""
   polynomial_interpolation(t, U₁ᴵ, Udᴵ)

IA interpolation of given points by 3rd degree polynomials . 
Returns coefficients `csc_V` as well as interval values `V` 
of the polynomials.

Arguments:
t            ... division points
Uᴵ           ... interval values to be interpolated
Udᴵ          ... interval first derivative values at division points
"""

function polynomial_interpolation(t, U₁ᴵ, Udᴵ)
   
    csc_V = Vector[]
 
    for i in 1:length(Udᴵ)-1
        a = b = c = d = 0
        
        a = (2*(U₁ᴵ[i+1]-U₁ᴵ[i]) + (Udᴵ[i]+Udᴵ[i+1])*(t[1]-t[2]))/(t[1]-t[2])^3
        b = (Udᴵ[i+1] - Udᴵ[i] + 3*a*(t[1]^2-t[2]^2))/(2*(t[2]-t[1]))
        c = Udᴵ[i] - t[1]*(2*b+t[1]*3*a)
        d = U₁ᴵ[i] - t[1]*(c+t[1]*(b+t[1]*a))
        
        append!(csc_V, [Interval{Float64}[a,b,c,d]])
    end
    
    # napočítat boxy výsledné fce
    V = Interval{Float64}[] 
    for i in 1:length(Udᴵ)-1
        x_int = t[i]..t[i+1]
        f(x) = csc_V[i][4] + (x-t[i])*(csc_V[i][3] + 
            (x-t[i])*(csc_V[i][2] + csc_V[i][1]*(x-t[i])))
        append!(V, f.(x_int))
    end
        
    return csc_V, V
end;



"""
   cubic_natural_spline(t, U, Uᴵ, Uₗd2, Uᵣd2; ns=10)

IA interpolation of given points by natural cubic spline. 
Returns spline coefficients `csc_V` as well as interval values `V` 
of the spline function.

Arguments:
t            ... division points
U, Uᴵ        ... values to be interpolated
Uₗd2, Uᵣd2    ... left and right boundary values of second derivative
ns=10        ... number of division points for each single piece of spline
"""
function cubic_natural_spline(t, U, Uᴵ, Uₗd2, Uᵣd2; ns=10)
    # A matrix
    n=length(Uᴵ)
    dv = [4..4 for i in 1:n-2]
    ev = [1..1 for i in 1:n-3]
    A = Array(SymTridiagonal(dv,ev))
    A⁻¹ = inv(A)

    # right-hand side
    h = 1.0/(n-1)
    rhs = []
    for i in 3:length(Uᴵ)
        append!(rhs, 6/h^2 * (Uᴵ[i] - 2 * Uᴵ[i-1] + Uᴵ[i-2]))
    end
    
    rhs[1] = rhs[1]-Uₗd2
    rhs[end] = rhs[end]-Uᵣd2    

    # second derivatives vector 
    Ud2 = []
    append!(Ud2, @interval(Uₗd2))
    append!(Ud2, A⁻¹*rhs)
    append!(Ud2, @interval(Uᵣd2))

    # spline coefficients
    csc_V = Vector[]
    for i in 1:length(Uᴵ)-1
        a = b = c = d = 0
        a = (Ud2[i+1]-Ud2[i])/(6*h)
        b = Ud2[i]/2 
        c = (Uᴵ[i+1] - Uᴵ[i])/h - h*(2*Ud2[i]+Ud2[i+1])/6
        d = Uᴵ[i]
        append!(csc_V, [Interval{Float64}[a,b,c,d]])
    end 
    
    V = Interval{Float64}[] 
    for i in 1:length(Uᴵ)-1
        x_dom = t[i]..t[i+1]
        x_int = mince(x_dom,ns)
        f(x) = csc_V[i][4] + (x-t[i])*(csc_V[i][3] + 
            (x-t[i])*(csc_V[i][2] + csc_V[i][1]*(x-t[i])))
        append!(V, f.(x_int))
    end
        
    return csc_V, V
end;



"""
    cubic_end_slope_spline(t, U, Uᴵ, Uₗd1, Uᵣd1; ns=10)

IA interpolation of given points with end slope cubic spline. 
Returns spline coefficients `csc_V`as well as interval values `V`.

Arguments:
t            ... division points
U, Uᴵ        ... values to be interpolated
Uₗd1, Uᵣd1    ... left and right boundary values of first derivative
ns=10        ... number of division points for each single piece of spline
"""

function cubic_end_slope_spline(t, U, Uᴵ, Uₗd1, Uᵣd1; ns=1)
    # A matrix
    n=length(Uᴵ)
    dv = [4..4 for i in 1:n-2]
    ev = [1..1 for i in 1:n-3]
    A = Array(SymTridiagonal(dv,ev))
    A[1,1] = 3.5..3.5
    A[end,end] = 3.5..3.5
    A⁻¹ = inv(A)

    # right-hand side
    h = 1.0/(n-1)
    rhs = []
    for i in 3:length(Uᴵ)
        append!(rhs, 6/h^2 * (Uᴵ[i] - 2 * Uᴵ[i-1] + Uᴵ[i-2]))
    end
    
    rhs[1] = rhs[1] - 3/h * ( (Uᴵ[2]-Uᴵ[1])/h - Uₗd1)
    rhs[end] = rhs[end] - 3/h * (Uᵣd1 - (Uᴵ[end]-Uᴵ[end-1])/h)    

    # second derivatives vector 
    sol = A⁻¹*rhs
    Ud2 = []
    
    σ₀ = 3/h * ( (Uᴵ[2]-Uᴵ[1])/h - Uₗd1) - sol[1]/2
    σ₁ = 3/h * (Uᵣd1 - (Uᴵ[end]-Uᴵ[end-1])/h)   - sol[end]/2
    append!(Ud2, @interval(σ₀))
    append!(Ud2, sol)
    append!(Ud2, @interval(σ₁))

    # spline coefficients
    csc_V = Vector[]
    for i in 1:length(Uᴵ)-1
        a=b=c=d=0
        a = (Ud2[i+1]-Ud2[i])/(6*h)
        b = Ud2[i]/2 
        c = (Uᴵ[i+1] - Uᴵ[i])/h - h*(2*Ud2[i]+Ud2[i+1])/6
        d = Uᴵ[i]
        append!(csc_V, Vector[Interval{Float64}[a,b,c,d]])
    end 
    
    V = Interval{Float64}[]
    for i in 1:length(Uᴵ)-1
        x_dom = t[i]..t[i+1]
        x_int = mince(x_dom,ns)
        f(x) = csc_V[i][4] + (x-t[i])*(csc_V[i][3] + 
            (x-t[i])*(csc_V[i][2] + csc_V[i][1]*(x-t[i])))
        append!(V, f.(x_int))
    end
        
    return csc_V, V
end;