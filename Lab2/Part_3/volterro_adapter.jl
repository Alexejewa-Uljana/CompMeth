using LinearAlgebra
using Plots

include("../Part_1/fredholm_solver.jl")

graphics_dir = "Lab2/Part_3/Part3_graphics"
if !isdir(graphics_dir)
    mkdir(graphics_dir)
end

function solve_volterra_simple(K::Function, f::Function, a::Float64, b::Float64, h::Float64, n_points=20)
    K_adapted(t, s) = (s <= t) ? K(t, s) : 0.0
    return fredholm_method(K_adapted, f, a, b, h, n_points)
end