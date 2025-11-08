using LinearAlgebra

include("../Part_1/fredholm_solver.jl")

function solve_fredholm_first_kind(K::Function, f::Function, a::Float64, b::Float64, h::Float64, alpha=0.1, n_points=20)
    if alpha <= 0
        error("Параметр регуляризации alpha должен быть положительным")
    end
    
    K_regularized(t, x) = - (1/alpha) * K(t, x)
    f_regularized(t) = (1/alpha) * f(t)
    
    return fredholm_method(K_regularized, f_regularized, a, b, h, n_points)
end

export solve_fredholm_first_kind