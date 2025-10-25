include("main.jl")
using Distributed
using Base.Threads

function parallel_variation_diminishing_approximation(f::Function, a::Float64, b::Float64, n_segments::Int)
    knots, n, h = uniform_quadratic_knots(a, b, n_segments)
    t_star = knot_averages(knots, 2)
    control_points = zeros(n)
    @threads for j in 1:n
        control_points[j] = f(t_star[j])
    end
    return QuadraticBSpline(knots, control_points, a, b, h)
end

function parallel_three_point_approximation(f::Function, a::Float64, b::Float64, n_segments::Int)
    knots, n, h = uniform_quadratic_knots(a, b, n_segments)
    t_star = knot_averages(knots, 2)
    A = zeros(n, n)
    b_vec = zeros(n)
    @threads for i in 1:n
        x = t_star[i]
        m = find_interval(knots, x, 2)
        if m != nothing
            for j in (m-2):m
                if j >= 1 && j <= n
                    A[i, j] = bspline_basis(j, 2, knots, x)
                end
            end
            b_vec[i] = f(x)
        end
    end
    control_points = pinv(A) * b_vec
    return QuadraticBSpline(knots, control_points, a, b, h)
end