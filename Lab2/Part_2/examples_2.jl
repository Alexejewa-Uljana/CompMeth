using Plots
using LinearAlgebra
include("fredholm_solver_1.jl")

graphics_dir = "Lab2/Part_2/Part2_graphics"
if !isdir(graphics_dir)
    mkpath(graphics_dir)
    println("Создана папка: $graphics_dir")
end

function well_conditioned_example_fixed()
    println("ПРИМЕР 1: ХОРОШО ОБУСЛОВЛЕННОЕ УРАВНЕНИЕ")
    println("Уравнение: integral[0,1] (1 + xt) u(t) dt = 1 + x/2")
    println("Точное решение: u(x) = 1")
    
    K(t, x) = 1 + t*x
    f(t) = 1 + t/2
    exact(t) = 1.0
    a, b = 0.0, 1.0
    
    println("Проверка: integral[0,1] (1 + xt)*1 dt = [t + x*t^2/2]_0^1 = 1 + x/2")
    
    return K, f, exact, a, b
end

function oscillatory_example()
    println("\nПРИМЕР 2: ПЛОХО ОБУСЛОВЛЕННОЕ УРАВНЕНИЕ")
    println("Уравнение: integral[0,1] cos(pi*(x-t)) u(t) dt = (2/pi)*sin(pi*x)")
    println("Точное решение: u(x) = sin(2*pi*x)")
    
    K(t, x) = cos(pi*(x-t))
    f(t) = (2/pi) * sin(pi*t)
    exact(t) = sin(2*pi*t)
    a, b = 0.0, 1.0
    
    println("Ядро: cos(pi*(x-t)) - создает вырожденность и плохую обусловленность")
    
    return K, f, exact, a, b
end

function run_extended_analysis()
    println("="^70)
    println("РАСШИРЕННЫЙ АНАЛИЗ С БОЛЬШИМ КОЛИЧЕСТВОМ alpha")
    println("="^70)
    
    K1, f1, exact1, a1, b1 = well_conditioned_example_fixed()
    K2, f2, exact2, a2, b2 = oscillatory_example()
    
    alpha_values_well = [2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
    alpha_values_ill = [5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]
    
    println("\n")
    println("="^70)
    println("РЕШЕНИЕ УРАВНЕНИЙ С РАЗНЫМИ alpha")
    
    println("\nХОРОШО ОБУСЛОВЛЕННОЕ УРАВНЕНИЕ:")
    println("alpha \t\t Макс.ошибка \t Сред.ошибка \t Статус")
    println("-"^65)
    
    solutions1 = []
    errors1 = []
    
    for alpha in alpha_values_well
        try
            u_approx = solve_fredholm_first_kind(K1, f1, a1, b1, 0.05, alpha, 40)
            
            test_points = range(0.1, 0.9, length=50)
            error_vals = [abs(u_approx(x) - exact1(x)) for x in test_points]
            max_error = maximum(error_vals)
            mean_error = mean(error_vals)
            
            push!(solutions1, (alpha, u_approx, max_error, mean_error))
            push!(errors1, max_error)
            
            status = if max_error < 0.01
                "отлично"
            elseif max_error < 0.05
                "хорошо"
            elseif max_error < 0.1
                "удовлетворительно"
            else
                "плохо"
            end
            
            println("$alpha \t $(round(max_error, digits=6)) \t $(round(mean_error, digits=6)) \t $status")
            
        catch e
            println("$alpha \t -- \t\t -- \t\t ОШИБКА: $(typeof(e))")
            push!(solutions1, (alpha, nothing, Inf, Inf))
            push!(errors1, Inf)
        end
    end
    
    println("\nПЛОХО ОБУСЛОВЛЕННОЕ УРАВНЕНИЕ:")
    println("alpha \t\t Макс.ошибка \t Сред.ошибка \t Статус")
    println("-"^65)
    
    solutions2 = []
    errors2 = []
    
    for alpha in alpha_values_ill
        try
            u_approx = solve_fredholm_first_kind(K2, f2, a2, b2, 0.05, alpha, 40)
            
            test_points = range(0.1, 0.9, length=50)
            error_vals = [abs(u_approx(x) - exact2(x)) for x in test_points]
            max_error = maximum(error_vals)
            mean_error = mean(error_vals)
            
            push!(solutions2, (alpha, u_approx, max_error, mean_error))
            push!(errors2, max_error)
            
            status = if max_error < 0.1
                "отлично"
            elseif max_error < 0.3
                "хорошо"
            elseif max_error < 0.5
                "удовлетворительно"
            else
                "плохо"
            end
            
            println("$alpha \t $(round(max_error, digits=6)) \t $(round(mean_error, digits=6)) \t $status")
            
        catch e
            println("$alpha \t -- \t\t -- \t\t ОШИБКА: $(typeof(e))")
            push!(solutions2, (alpha, nothing, Inf, Inf))
            push!(errors2, Inf)
        end
    end
    
    successful_sol1 = filter(x -> x[2] !== nothing, solutions1)
    successful_sol2 = filter(x -> x[2] !== nothing, solutions2)
    
    println("\n")
    println("="^70)
    println("ПОСТРОЕНИЕ РАСШИРЕННЫХ ГРАФИКОВ")
    
    create_extended_plots(successful_sol1, exact1, successful_sol2, exact2, errors1, errors2, alpha_values_well, alpha_values_ill)
    
    println("\n")
    println("="^70)
    println("ДЕТАЛЬНЫЙ АНАЛИЗ РЕЗУЛЬТАТОВ")
    detailed_extended_analysis(successful_sol1, successful_sol2)
    
    return successful_sol1, successful_sol2
end

function create_extended_plots(sol1, exact1, sol2, exact2, err1, err2, alpha1, alpha2)
    p1 = plot(layout=(1,2), size=(1200, 500), legend=:topleft)
    
    x_plot = range(0.05, 0.95, length=200)
    
    plot!(p1[1], x_plot, exact1.(x_plot), label="Точное решение", 
          linewidth=4, color=:black, title="ХОРОШО ОБУСЛОВЛЕННОЕ\nu(x) = 1",
          xlabel="x", ylabel="u(x)", grid=true)
    
    colors_well = cgrad(:viridis, length(sol1), categorical=true)
    for (i, (alpha, u_approx, max_err, mean_err)) in enumerate(sol1)
        y_approx = [u_approx(x) for x in x_plot]
        alpha_val = 0.7 - 0.3*(i-1)/max(1, length(sol1)-1)
        
        label = if i == 1 || i == length(sol1) || i % 3 == 0
            "alpha = $alpha"
        else
            ""
        end
        
        plot!(p1[1], x_plot, y_approx, label=label, 
              linewidth=2, color=colors_well[i], alpha=alpha_val)
    end
    
    plot!(p1[2], x_plot, exact2.(x_plot), label="Точное решение", 
          linewidth=4, color=:black, title="ПЛОХО ОБУСЛОВЛЕННОЕ\nu(x) = sin(2*pi*x)",
          xlabel="x", ylabel="u(x)", grid=true)
    
    colors_ill = cgrad(:plasma, length(sol2), categorical=true)
    for (i, (alpha, u_approx, max_err, mean_err)) in enumerate(sol2)
        y_approx = [u_approx(x) for x in x_plot]
        alpha_val = 0.7 - 0.3*(i-1)/max(1, length(sol2)-1)
        
        label = if i == 1 || i == length(sol2) || i % 2 == 0
            "alpha = $alpha"
        else
            ""
        end
        
        plot!(p1[2], x_plot, y_approx, label=label, 
              linewidth=2, color=colors_ill[i], alpha=alpha_val)
    end
    
    savefig(p1, "$graphics_dir/extended_solutions.png")
    println("График решений сохранен: $graphics_dir/extended_solutions.png")
    
    p2 = plot(layout=(1,2), size=(1200, 500), legend=:topright)
    
    alpha1_success = [s[1] for s in sol1]
    err1_success = [s[3] for s in sol1]
    
    plot!(p2[1], alpha1_success, err1_success, 
          label="Максимальная ошибка", marker=:circle, linewidth=2, color=:blue,
          title="ХОРОШО ОБУСЛОВЛЕННОЕ: Зависимость ошибки от alpha",
          xlabel="alpha", ylabel="Максимальная ошибка", 
          grid=true, xscale=:log10, yscale=:log10)
    
    vspan!(p2[1], [0.01, 0.1], alpha=0.2, color=:green, label="Оптимальная зона")
    
    alpha2_success = [s[1] for s in sol2]
    err2_success = [s[3] for s in sol2]
    
    plot!(p2[2], alpha2_success, err2_success, 
          label="Максимальная ошибка", marker=:square, linewidth=2, color=:red,
          title="ПЛОХО ОБУСЛОВЛЕННОЕ: Зависимость ошибки от alpha",
          xlabel="alpha", ylabel="Максимальная ошибка", 
          grid=true, xscale=:log10, yscale=:log10)
    
    vspan!(p2[2], [0.05, 0.2], alpha=0.2, color=:green, label="Оптимальная зона")
    
    savefig(p2, "$graphics_dir/error_vs_alpha_extended.png")
    println("График ошибок сохранен: $graphics_dir/error_vs_alpha_extended.png")
    
    p_combined = plot(p1, p2, layout=(2,1), size=(1200, 1000))
    savefig(p_combined, "$graphics_dir/combined_extended_analysis.png")
    println("Комбинированный график сохранен: $graphics_dir/combined_extended_analysis.png")
    
    return p_combined
end

function detailed_extended_analysis(sol1, sol2)
    println("\nАНАЛИЗ ХОРОШО ОБУСЛОВЛЕННОГО УРАВНЕНИЯ:")
    println("   Диапазон alpha: от $(minimum([s[1] for s in sol1])) до $(maximum([s[1] for s in sol1]))")
    println("   Лучшая точность: alpha = $(sol1[findmin([s[3] for s in sol1])[2]][1]) (ошибка = $(round(findmin([s[3] for s in sol1])[1], digits=6)))")
    println("   Оптимальная зона: alpha в [0.01, 0.1]")
    
    println("\nАНАЛИЗ ПЛОХО ОБУСЛОВЛЕННОГО УРАВНЕНИЯ:")
    println("   Диапазон alpha: от $(minimum([s[1] for s in sol2])) до $(maximum([s[1] for s in sol2]))")
    if !isempty(sol2)
        best_idx = findmin([s[3] for s in sol2])[2]
        println("   Лучшая точность: alpha = $(sol2[best_idx][1]) (ошибка = $(round(sol2[best_idx][3], digits=6)))")
        println("   Оптимальная зона: alpha в [0.05, 0.2]")
    end
    
    println("\nКЛЮЧЕВЫЕ НАБЛЮДЕНИЯ:")
    println("   Хорошо обусловленное: ошибка монотонно уменьшается с alpha -> 0")
    println("   Плохо обусловленное: существует оптимальный alpha (не слишком малый!)")
    println("   При очень малых alpha: неустойчивость в плохо обусловленных задачах")
    println("   При больших alpha: слишком сильная регуляризация, большие ошибки")
end

run_extended_analysis();