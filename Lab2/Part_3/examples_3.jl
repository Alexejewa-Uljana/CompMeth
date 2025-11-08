include("volterro_adapter.jl")

function analyze_volterra_performance()
    println("АНАЛИЗ ПРОИЗВОДИТЕЛЬНОСТИ АДАПТАЦИИ")
    println("="^50)
    K(t, s) = exp(t - s)
    f(t) = exp(t)
    exact(t) = 1.0
    a, b = 0.0, 1.0
    println("\nЗависимость времени от размера сетки:")
    println("N \t\t Время (сек) \t Ошибка")
    println("-"^45)
    grid_sizes = [10, 20, 30, 40, 50]
    times = []
    errors = []
    for n in grid_sizes
        h = (b - a) / n
        start_time = time()
        u_approx = solve_volterra_simple(K, f, a, b, h, 20)
        end_time = time()
        elapsed_time = end_time - start_time
        test_points = range(0.2, 0.8, length=10)
        max_error = maximum([abs(u_approx(t) - exact(t)) for t in test_points])
        push!(times, elapsed_time)
        push!(errors, max_error)
        println("$n \t $(round(elapsed_time, digits=4)) \t $(round(max_error, digits=6))")
    end
    return grid_sizes, times, errors
end


function test_volterra_adaptation()
    println("ТЕСТИРОВАНИЕ АДАПТАЦИИ ДЛЯ УРАВНЕНИЙ ВОЛЬТЕРРА")
    println("="^50)
    println("\nПРИМЕР 1: u(t) = e^t - ∫₀ᵗ e^(t-s) u(s) ds")
    println("Точное решение: u(t) = 1")
    K1(t, s) = exp(t - s)
    f1(t) = exp(t)
    exact1(t) = 1.0
    a, b = 0.0, 1.0
    solutions_ex1 = []
    for h in [0.2, 0.1, 0.05]
        try
            u_approx = solve_volterra_simple(K1, f1, a, b, h, 30)
            push!(solutions_ex1, (h, u_approx))
            println("h = $h: решение получено")
        catch e
            println("h = $h: ошибка")
        end
    end
    println("\nПРИМЕР 2: u(t) = t² + ∫₀ᵗ (t-s)u(s) ds")
    println("Точное решение: u(t) = 2(cosh(t) - 1)")
    K2(t, s) = t - s
    f2(t) = t^2
    exact2(t) = 2*(cosh(t) - 1)
    u_approx2 = solve_volterra_simple(K2, f2, a, b, 0.1, 30)
    solutions_ex2 = [(0.1, u_approx2)]
    return solutions_ex1, solutions_ex2
end

function create_combined_plot(solutions_ex1, solutions_ex2, grid_sizes, times, errors)
    println("\nСОЗДАНИЕ КОМБИНИРОВАННОГО ГРАФИКА")
    println("="^50)
    p = plot(layout=(2,2), size=(1000, 800), legend=:topleft)
    x_plot = range(0.1, 0.9, length=100)
    plot!(p[1], x_plot, ones(length(x_plot)), 
          label="Точное решение", linewidth=3, color=:black, linestyle=:dash,
          title="(A) Пример 1: u(t) = 1",
          xlabel="t", ylabel="u(t)", grid=true)
    colors = [:blue, :red, :green]
    for (i, (h, u_approx)) in enumerate(solutions_ex1)
        if i <= length(colors)
            y_approx = [u_approx(x) for x in x_plot]
            plot!(p[1], x_plot, y_approx, label="h = $h", 
                  linewidth=2, color=colors[i])
        end
    end
    x_plot2 = range(0.0, 1.0, length=100)
    exact_vals = [2*(cosh(t) - 1) for t in x_plot2]
    plot!(p[2], x_plot2, exact_vals,
          label="Точное решение", linewidth=3, color=:black, linestyle=:dash,
          title="(B) Пример 2: u(t) = 2(cosh(t)-1)", 
          xlabel="t", ylabel="u(t)", grid=true)
    for (h, u_approx) in solutions_ex2
        y_approx = [u_approx(x) for x in x_plot2]
        plot!(p[2], x_plot2, y_approx, label="h = $h", 
              linewidth=2, color=:red)
    end
    plot!(p[3], grid_sizes, times,
          label="Время решения", marker=:circle, linewidth=2, color=:blue,
          title="(C) Зависимость времени от размера сетки",
          xlabel="Размер сетки N", ylabel="Время (сек)", grid=true)
    plot!(p[4], grid_sizes, errors,
          label="Максимальная ошибка", marker=:square, linewidth=2, color=:red,
          title="(D) Зависимость ошибки от размера сетки",
          xlabel="Размер сетки N", ylabel="Максимальная ошибка", grid=true)
    savefig(p, "$graphics_dir/volterra_combined_analysis.png")
    println("Комбинированный график сохранен: $graphics_dir/volterra_combined_analysis.png")
    return p
end


function run_volterra_analysis()
    println("ЧАСТЬ 3: АНАЛИЗ АДАПТАЦИИ ДЛЯ УРАВНЕНИЙ ВОЛЬТЕРРА")
    println("="^70)
    solutions_ex1, solutions_ex2 = test_volterra_adaptation()
    grid_sizes, times, errors = analyze_volterra_performance()
    create_combined_plot(solutions_ex1, solutions_ex2, grid_sizes, times, errors)
    println("\n")
    println("="^70)
    println("ВЫВОДЫ ПО АДАПТАЦИИ:")
    println("="^70)
    println("Адаптация ВОЗМОЖНА через модификацию ядра")
    println("Метод работает для небольших задач")
    println("Вычислительно неэффективен для больших систем") 
    
    println("\n Комбинированный график сохранен: $graphics_dir/volterra_combined_analysis.png")
    println("   (A) Пример 1 - простое уравнение")
    println("   (B) Пример 2 - уравнение с аналитическим решением")
    println("   (C) Зависимость времени от размера сетки")
    println("   (D) Зависимость ошибки от размера сетки")
    
    println("\n АНАЛИЗ ЧАСТИ 3 ЗАВЕРШЕН")
end

export solve_volterra_simple, run_volterra_analysis

run_volterra_analysis()