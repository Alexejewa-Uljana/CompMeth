include("fredholm_solver.jl")

graphics_dir = "Lab2/Part_1/Part1_graphics"
if !isdir(graphics_dir)
    mkdir(graphics_dir)
    println("Создана папка: $graphics_dir")
end

function test1()
    println("="^50)
    println("ТЕСТ 1: Уравнение с экспоненциальным ядром")
    
    lambda = 0.5
    K(t, x) = lambda * exp(t - x)
    f(t) = (1 - lambda) * exp(t)
    exact(t) = exp(t)
    a, b = 0.0, 1.0
    
    solutions = []
    errors = []
    
    for h in [0.2, 0.1, 0.05]
        try
            u_approx = fredholm_method(K, f, a, b, h, 20)
            push!(solutions, (h, u_approx))
            
            t_test = 0.5
            exact_val = exact(t_test)
            approx_val = u_approx(t_test)
            error = abs(exact_val - approx_val)
            push!(errors, error)
            
            println("h = $h, u(0.5): точное = $exact_val, приближенное = $approx_val, ошибка = $error")
        catch e
            println("Ошибка при h = $h: $e")
        end
    end
    
    if !isempty(solutions)
        plot_test1(solutions, exact)
    else
        println("Нет успешных решений для построения графиков")
    end
    
    return solutions, errors
end

function test2()
    println("="^50)
    println("ТЕСТ 2: Уравнение с полиномиальным ядром")
    
    K(t, x) = t * x
    f(t) = (2/3) * t
    exact(t) = t
    a, b = 0.0, 1.0
    
    solutions = []
    errors = []
    
    for h in [0.2, 0.1, 0.05]
        try
            u_approx = fredholm_method(K, f, a, b, h, 20)
            push!(solutions, (h, u_approx))
            
            t_test = 0.7
            exact_val = exact(t_test)
            approx_val = u_approx(t_test)
            error = abs(exact_val - approx_val)
            push!(errors, error)
            
            println("h = $h, u(0.7): точное = $exact_val, приближенное = $approx_val, ошибка = $error")
        catch e
            println("Ошибка при h = $h: $e")
        end
    end
    
    if !isempty(solutions)
        plot_test2(solutions, exact)
    else
        println("Нет успешных решений для построения графиков")
    end
    
    return solutions, errors
end

function plot_test1(solutions, exact)
    p = plot(legend=:topleft, title="ТЕСТ 1: Экспоненциальное ядро\nu(t) = e^t", 
             xlabel="t", ylabel="u(t)", grid=true, size=(800, 400))
    
    x_plot = range(0.1, 0.9, length=100)
    
    plot!(p, x_plot, exact.(x_plot), label="Точное решение", 
          linewidth=3, color=:black, linestyle=:dash)
    
    colors = [:blue, :red, :green]
    for (i, (h, u_approx)) in enumerate(solutions)
        y_approx = [u_approx(x) for x in x_plot]
        plot!(p, x_plot, y_approx, label="h = $h", 
              linewidth=2, color=colors[i])
    end
    
    savefig(p, "$graphics_dir/test1_exponential.png")
    println("График для теста 1 сохранен: $graphics_dir/test1_exponential.png")
    
    return p
end

function plot_test2(solutions, exact)
    p = plot(legend=:topleft, title="ТЕСТ 2: Полиномиальное ядро\nu(t) = t", 
             xlabel="t", ylabel="u(t)", grid=true, size=(800, 400))
    
    x_plot = range(0.1, 0.9, length=100)
    
    plot!(p, x_plot, exact.(x_plot), label="Точное решение", 
          linewidth=3, color=:black, linestyle=:dash)
    
    colors = [:blue, :red, :green]
    for (i, (h, u_approx)) in enumerate(solutions)
        y_approx = [u_approx(x) for x in x_plot]
        plot!(p, x_plot, y_approx, label="h = $h", 
              linewidth=2, color=colors[i])
    end
    savefig(p, "$graphics_dir/test2_polynomial.png")
    println("График для теста 2 сохранен: $graphics_dir/test2_polynomial.png")
    return p
end

function plot_convergence(solutions1, errors1, solutions2, errors2)
    p = plot(layout=(1,2), size=(1000, 400), legend=:topright)
    h_values1 = [s[1] for s in solutions1]
    plot!(p[1], h_values1, errors1, 
          label="Экспоненциальное ядро", marker=:circle, linewidth=2, color=:blue,
          title="Сходимость: тест 1 (экспоненциальное ядро)",
          xlabel="Шаг сетки h", ylabel="Ошибка в точке t=0.5", 
          grid=true, xscale=:log10, yscale=:log10)
    h_values2 = [s[1] for s in solutions2]
    plot!(p[2], h_values2, errors2, 
          label="Полиномиальное ядро", marker=:square, linewidth=2, color=:red,
          title="Сходимость: тест 2 (полиномиальное ядро)",
          xlabel="Шаг сетки h", ylabel="Ошибка в точке t=0.7", 
          grid=true, xscale=:log10, yscale=:log10)
    savefig(p, "$graphics_dir/convergence_analysis.png")
    println("График сходимости сохранен: $graphics_dir/convergence_analysis.png")
    return p
end
function plot_combined(solutions1, exact1, solutions2, exact2)
    p = plot(layout=(2,1), size=(800, 600), legend=:topleft)
    x_plot = range(0.1, 0.9, length=100)
    plot!(p[1], x_plot, exact1.(x_plot), label="Точное решение", 
          linewidth=3, color=:black, linestyle=:dash,
          title="ТЕСТ 1: Экспоненциальное ядро - u(t) = e^t",
          xlabel="", ylabel="u(t)", grid=true)
    colors = [:blue, :red, :green]
    for (i, (h, u_approx)) in enumerate(solutions1)
        y_approx = [u_approx(x) for x in x_plot]
        plot!(p[1], x_plot, y_approx, label="h = $h", 
              linewidth=2, color=colors[i])
    end
    plot!(p[2], x_plot, exact2.(x_plot), label="Точное решение", 
          linewidth=3, color=:black, linestyle=:dash,
          title="ТЕСТ 2: Полиномиальное ядро - u(t) = t",
          xlabel="t", ylabel="u(t)", grid=true)
    for (i, (h, u_approx)) in enumerate(solutions2)
        y_approx = [u_approx(x) for x in x_plot]
        plot!(p[2], x_plot, y_approx, label="h = $h", 
              linewidth=2, color=colors[i])
    end
    savefig(p, "$graphics_dir/combined_solutions.png")
    println("Комбинированный график сохранен: $graphics_dir/combined_solutions.png")
    return p
end

function run_all_tests()
    println("ЗАПУСК ТЕСТОВ МЕТОДА СПЛАЙН-КОЛЛОКАЦИЙ")
    println("="^60)
    solutions1, errors1 = test1() 
    solutions2, errors2 = test2()
    if !isempty(solutions1) && !isempty(solutions2)
        plot_convergence(solutions1, errors1, solutions2, errors2)
        plot_combined(solutions1, t -> exp(t), solutions2, t -> t)
        println("="^60)
        println("ГРАФИКИ СОХРАНЕНЫ В ПАПКЕ: $graphics_dir/")
        println("test1_exponential.png - тест 1 (экспоненциальное ядро)")
        println("test2_polynomial.png - тест 2 (полиномиальное ядро)") 
        println("convergence_analysis.png - анализ сходимости")
        println("combined_solutions.png - комбинированный график")
    else
        println("Невозможно построить графики - недостаточно успешных решений")
    end
    println("="^60)
    println("ТЕСТИРОВАНИЕ ЗАВЕРШЕНО")
end

run_all_tests()