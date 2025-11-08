using LinearAlgebra

include("integral.jl")

function fredholm_method(K::Function, f::Function, a::Float64, b::Float64, h::Float64, n_points=20)
    knots = get_points(a, b, h)
    n = length(knots) - 3
    rhs_vector = [get_coef(f, knots, idx) for idx in 1:n]
    gram_mat = compute_omega_matrix(knots, a, b, n_points)
    op_matrix = zeros(n, n)
    for col in 1:n
        transformed_func = create_kernel_transform(K, gram_mat, knots, col)
        for row in 1:n
            op_matrix[row, col] = get_coef(transformed_func, knots, row)
        end
    end
    
    system_matrix = Matrix(I, n, n) - op_matrix
    coeffs = system_matrix \ rhs_vector
    
    return x -> sum(coeffs[j] * quad_bspline(x, j, knots) for j in 1:n)
end

function create_kernel_transform(K::Function, gram_mat::Matrix{Float64}, knots::Vector{Float64}, j::Int)
    function kernel_transform(t::Float64)
        total = 0.0
        for k in 1:size(gram_mat, 1)
            coef_func = x -> K(t, x)
            basis_coef = get_coef(coef_func, knots, k)
            total += basis_coef * gram_mat[k, j]
        end
        return total
    end
    return kernel_transform
end