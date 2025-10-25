include("main.jl")
using Printf


function compute_errors(f::Function, a::Float64, b::Float64, n_segments::Int)
    spline_vd = variation_diminishing_approximation(f, a, b, n_segments)
    spline_3pt = three_point_approximation(f, a, b, n_segments)
    n_test_points = 10 * n_segments + 1
    x_test = range(a, b, length=n_test_points)
    max_error_vd = maximum(abs(f(x) - evaluate_spline(spline_vd, x)) for x in x_test)
    max_error_3pt = maximum(abs(f(x) - evaluate_spline(spline_3pt, x)) for x in x_test)
    mse_vd = sqrt(mean((f(x) - evaluate_spline(spline_vd, x))^2 for x in x_test))
    mse_3pt = sqrt(mean((f(x) - evaluate_spline(spline_3pt, x))^2 for x in x_test))
    return (max_error_vd, max_error_3pt, mse_vd, mse_3pt)
end

test_functions = [
    (x -> x^2, "x^2", "Гладкая квадратичная функция"),
    (x -> sin(2*pi*x), "sin(2pi x)", "Гладкая периодическая функция"), 
    (x -> abs(x - 0.5), "|x-0.5|", "Функция с изломом"),
    (x -> x < 0.5 ? 0.0 : 1.0, "H(x-0.5)", "Ступенчатая функция"),
    (x -> exp(-50*(x-0.5)^2), "exp(-50(x-0.5)^2)", "Резкий пик")
]
a, b = 0.0, 1.0
n_segments_list = [4, 8, 16, 32]
println("ЭКСПЕРИМЕНТ: Сравнение методов аппроксимации")
println("="^60)
for (f, f_name, f_description) in test_functions
    println("\nФункция: $f_name")
    println("Описание: $f_description")
    println("-"^40)
    for n_segments in n_segments_list
        max_err_vd, max_err_3pt, mse_vd, mse_3pt = compute_errors(f, a, b, n_segments)  
        println("Сегментов: $n_segments")
        println("  Уменьшение вариации:  max_error=$(@sprintf "%.2e" max_err_vd), MSE=$(@sprintf "%.2e" mse_vd)")
        println("  Трехточечный метод:   max_error=$(@sprintf "%.2e" max_err_3pt), MSE=$(@sprintf "%.2e" mse_3pt)")
        println()
    end
end