using LinearAlgebra
using Plots
using Statistics

struct QuadraticBSpline
    knots::Vector{Float64}
    control_points::Vector{Float64}
    a::Float64
    b::Float64
    h::Float64
end

function uniform_quadratic_knots(a::Float64, b::Float64, n_segments::Int)
    h = (b - a) / n_segments
    n = n_segments + 2
    knots = Vector{Float64}(undef, n + 3)
    knots[1] = a
    knots[2] = a
    knots[3] = a
    for i in 1:n_segments
        knots[3 + i] = a + i * h
    end
    knots[end-2] = b
    knots[end-1] = b
    knots[end] = b
    return knots, n, h
end

function bspline_basis(i::Int, d::Int, knots::Vector{Float64}, x::Float64)
    if d == 0
        if knots[i] <= x && x < knots[i+1]
            return 1.0
        else
            return 0.0
        end
    else
        term1 = 0.0
        term2 = 0.0
        denom1 = knots[i+d] - knots[i]
        denom2 = knots[i+d+1] - knots[i+1]
        if denom1 > 0
            term1 = (x - knots[i]) / denom1 * bspline_basis(i, d-1, knots, x)
        end
        if denom2 > 0
            term2 = (knots[i+d+1] - x) / denom2 * bspline_basis(i+1, d-1, knots, x)
        end
        return term1 + term2
    end
end

function knot_averages(knots::Vector{Float64}, d::Int)
    n = length(knots) - d - 1
    t_star = Vector{Float64}(undef, n)
    for j in 1:n
        t_star[j] = (knots[j+1] + knots[j+2]) / 2.0
    end
    return t_star
end

function evaluate_spline(spline::QuadraticBSpline, x::Float64)
    knots = spline.knots
    coeffs = spline.control_points
    d = 2
    m = find_interval(knots, x, d)
    if m == nothing
        return 0.0
    end
    result = 0.0
    for j in (m-d):m
        if j >= 1 && j <= length(coeffs)
            basis_val = bspline_basis(j, d, knots, x)
            result += coeffs[j] * basis_val
        end
    end
    return result
end

function find_interval(knots::Vector{Float64}, x::Float64, d::Int)
    for m in (d+1):(length(knots)-1)
        if knots[m] <= x && x <= knots[m+1]
            return m
        end
    end
    return nothing
end

function variation_diminishing_approximation(f::Function, a::Float64, b::Float64, n_segments::Int)
    knots, n, h = uniform_quadratic_knots(a, b, n_segments)
    t_star = knot_averages(knots, 2)
    control_points = zeros(n)
    for j in 1:n
        control_points[j] = f(t_star[j])
    end
    return QuadraticBSpline(knots, control_points, a, b, h)
end

function three_point_approximation(f::Function, a::Float64, b::Float64, n_segments::Int)
    knots, n, h = uniform_quadratic_knots(a, b, n_segments)
    t_star = knot_averages(knots, 2)
    A = zeros(n, n)
    b_vec = zeros(n)
    for i in 1:n
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