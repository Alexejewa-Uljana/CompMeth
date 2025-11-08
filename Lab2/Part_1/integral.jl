using LinearAlgebra

include("utils.jl")

function integrate_simpson(f, a, b, n_points=50)
    h = (b - a) / n_points
    result = f(a) + f(b)
    for i in 1:n_points-1
        x = a + i * h
        if i % 2 == 0
            result += 2 * f(x)
        else
            result += 4 * f(x)
        end
    end
    return result * h / 3
end

function compute_spline_product_integral(i::Int, j::Int, knots::Vector{Float64}, a::Float64, b::Float64, n_points)
    function temp_function(x)
        return quad_bspline(x, i, knots) * quad_bspline(x, j, knots)
    end
    integral = temp_function
    support_start = max(knots[i], knots[j])
    support_end = min(knots[i+3], knots[j+3])
    integration_start = max(a, support_start)
    integration_end = min(b, support_end)
    if integration_start >= integration_end
        return 0.0
    end
    return integrate_simpson(integral, integration_start, integration_end, n_points)
end

function compute_omega_matrix(knots::Vector{Float64}, a::Float64, b::Float64, n_points)
    total_splines = length(knots) - 2 - 1
    gram_matrix = zeros(total_splines, total_splines)
    for spline_i in 1:total_splines
        for spline_j in 1:total_splines
            integral_value = compute_spline_product_integral(spline_i, spline_j, knots, a, b, n_points)
            gram_matrix[spline_i, spline_j] = integral_value
        end
    end
    return gram_matrix
end