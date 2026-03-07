using DrWatson
@quickactivate "project"

using DifferentialEquations
using DataFrames
using StatsPlots
using LaTeXStrings
using Plots
using Statistics
using FFTW

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

# Описание модели Лотки-Вольтерры
"""
Модель Лотки-Вольтерры (хищник-жертва)
Система уравнений:
dx/dt = αx - βxy  Изменение популяции жертв
dy/dt = δxy - γy  Изменение популяции хищников
Где:
x - популяция жертв (например, зайцы)
y - популяция хищников (например, лисы)
α - естественный прирост жертв (в отсутствие хищников)
β - коэффициент поедания жертв хищниками
δ - коэффициент прироста хищников за счет поедания жертв
γ - естественная смертность хищников (в отсутствие жертв)
"""

function lotka_volterra!(du, u, p, t)
    x, y = u 
    α, β, δ, γ = p 
    @inbounds begin
        du[1] = α*x - β*x*y 
        du[2] = δ*x*y - γ*y 
    end
    nothing
end

# Начальные условия и временные параметры
u0_lv = [40.0, 9.0] # начальная популяция [жертвы, хищники]
tspan_lv = (0.0, 200.0) # длительность симуляции
dt_lv = 0.01 # шаг интегрирования

# СОЗДАЕМ НАБОР ПАРАМЕТРОВ ДЛЯ ЦИКЛИЧЕСКОГО ПЕРЕБОРА
# Каждый набор параметров: [α, β, δ, γ]
parameter_sets = [
    [0.20, 0.06, 0.04, 0.10],  # Классический набор (базовый)
    [0.30, 0.06, 0.04, 0.10],  # Увеличенный рост жертв
    [0.10, 0.06, 0.04, 0.10],  # Уменьшенный рост жертв
    [0.20, 0.08, 0.04, 0.10],  # Увеличенное поедание жертв
    [0.20, 0.04, 0.04, 0.10],  # Уменьшенное поедание жертв
    [0.20, 0.06, 0.06, 0.10],  # Увеличенная конверсия в хищников
    [0.20, 0.06, 0.02, 0.10],  # Уменьшенная конверсия в хищников
    [0.20, 0.06, 0.04, 0.15],  # Увеличенная смертность хищников
    [0.20, 0.06, 0.04, 0.05],  # Уменьшенная смертность хищников
]

# Сохраняем все результаты
results = []

# Функция для нахождения первого пика
function find_first_peak(signal, time)
    for i in 2:length(signal)-1
        if signal[i] > signal[i-1] && signal[i] > signal[i+1]
            return time[i], signal[i]
        end
    end
    return NaN, NaN
end

# Функция для вычисления спектра
function compute_fft(signal, dt)
    n = length(signal)
    spectrum = abs.(rfft(signal))
    freq = rfftfreq(n, 1/dt)
    return freq, spectrum
end

# Проходим циклом по каждому набору параметров
for (idx, p) in enumerate(parameter_sets)
    println("\n" * "="^60)
    println("НАБОР ПАРАМЕТРОВ $idx")
    println("="^60)
    println("Параметры модели:")
    println("α (скорость размножения жертв) = ", p[1])
    println("β (скорость поедания жертв) = ", p[2])
    println("δ (коэффициент конверсии) = ", p[3])
    println("γ (смертность хищников) = ", p[4])
    x_star = p[4] / p[3] # стационарная точка для жертв
    y_star = p[1] / p[2] # стационарная точка для хищников
    println("\nСтационарные точки (положения равновесия):")
    println("x* = γ/δ = ", round(x_star, digits=3))
    println("y* = α/β = ", round(y_star, digits=3))
    prob_lv = ODEProblem(lotka_volterra!, u0_lv, tspan_lv, p)
    sol_lv = solve(prob_lv,
                   dt = dt_lv,
                   Tsit5(),
                   reltol=1e-8,
                   abstol=1e-10,
                   saveat=0.1,
                   dense=true)
    df_lv = DataFrame()
    df_lv[!, :t] = sol_lv.t
    df_lv[!, :prey] = [u[1] for u in sol_lv.u] 
    df_lv[!, :predator] = [u[2] for u in sol_lv.u] 
    df_lv[!, :набор] .= "Набор $idx"
    df_lv[!, :dprey_dt] = p[1] .* df_lv.prey .- p[2] .* df_lv.prey .* df_lv.predator
    df_lv[!, :dpredator_dt] = p[3] .* df_lv.prey .* df_lv.predator .- p[4] .* df_lv.predator
    df_lv[!, :prey_pct_change] = df_lv.dprey_dt ./ df_lv.prey .* 100
    df_lv[!, :predator_pct_change] = df_lv.dpredator_dt ./ df_lv.predator .* 100
    df_lv[!, :x_star] .= x_star
    df_lv[!, :y_star] .= y_star
    push!(results, df_lv)
    println("\nОсновные статистики:")
    println("Жертвы: min = ", round(minimum(df_lv.prey), digits=2),
            ", max = ", round(maximum(df_lv.prey), digits=2),
            ", mean = ", round(mean(df_lv.prey), digits=2))
    println("Хищники: min = ", round(minimum(df_lv.predator), digits=2),
            ", max = ", round(maximum(df_lv.predator), digits=2),
            ", mean = ", round(mean(df_lv.predator), digits=2))
    peak_time_prey, peak_value_prey = find_first_peak(df_lv.prey, df_lv.t)
    peak_time_predator, peak_value_predator = find_first_peak(df_lv.predator, df_lv.t)
    
    if !isnan(peak_time_prey) && !isnan(peak_time_predator)
        phase_shift = peak_time_predator - peak_time_prey
        println("\nАнализ колебаний:")
        println("Первый пик жертв: время = ", round(peak_time_prey, digits=2), 
                ", значение = ", round(peak_value_prey, digits=2))
        println("Первый пик хищников: время = ", round(peak_time_predator, digits=2), 
                ", значение = ", round(peak_value_predator, digits=2))
        println("Сдвиг фаз (хищники отстают): ", round(phase_shift, digits=2))
    end
    freq_prey, spectrum_prey = compute_fft(df_lv.prey .- mean(df_lv.prey), dt_lv)
    if length(spectrum_prey) > 0
        idx_prey = argmax(spectrum_prey[2:end]) + 1
        dominant_freq_prey = freq_prey[idx_prey]
        period_prey = 1/dominant_freq_prey
        println("Доминирующий период колебаний жертв: ", round(period_prey, digits=2), " единиц времени")
    end
end

# Сохраняем все результаты
using Serialization
serialize(datadir(script_name, "lv_all_results.jls"), results)

# Объединяем результаты для сравнительных графиков
combined_df = vcat(results...)

# СОЗДАЕМ СРАВНИТЕЛЬНЫЕ ГРАФИКИ

# График 1: Сравнение динамики жертв для всех наборов
plt_prey_comparison = plot(xlabel="Время", ylabel="Популяция жертв", 
                           title="Сравнение динамики жертв для разных параметров", 
                           size=(900, 500), legend=:topright)
for (idx, df) in enumerate(results)
    plot!(plt_prey_comparison, df.t, df.prey, 
          label="Набор $idx (α=$(round(parameter_sets[idx][1], digits=2)), γ=$(round(parameter_sets[idx][4], digits=2)))", 
          linewidth=1.5)
end
savefig(plt_prey_comparison, plotsdir(script_name, "lv_prey_comparison.png"))

# График 2: Сравнение динамики хищников для всех наборов
plt_predator_comparison = plot(xlabel="Время", ylabel="Популяция хищников", 
                               title="Сравнение динамики хищников для разных параметров", 
                               size=(900, 500), legend=:topright)
for (idx, df) in enumerate(results)
    plot!(plt_predator_comparison, df.t, df.predator, 
          label="Набор $idx (β=$(round(parameter_sets[idx][2], digits=2)), δ=$(round(parameter_sets[idx][3], digits=2)))", 
          linewidth=1.5)
end
savefig(plt_predator_comparison, plotsdir(script_name, "lv_predator_comparison.png"))

# График 3: Сравнение фазовых портретов
plt_phase_comparison = plot(xlabel="Популяция жертв (x)", ylabel="Популяция хищников (y)", 
                            title="Сравнение фазовых портретов", 
                            size=(800, 600), legend=:topright)
for (idx, df) in enumerate(results)
    plot!(plt_phase_comparison, df.prey, df.predator, 
          label="Набор $idx", linewidth=1.5, alpha=0.7)
    scatter!(plt_phase_comparison, [df.x_star[1]], [df.y_star[1]], 
             markersize=5, label="Равновесие $idx", alpha=0.7)
end
savefig(plt_phase_comparison, plotsdir(script_name, "lv_phase_comparison.png"))

# График 4: Сравнение периодов колебаний
periods = []
for (idx, df) in enumerate(results)
    freq_prey, spectrum_prey = compute_fft(df.prey .- mean(df.prey), dt_lv)
    if length(spectrum_prey) > 0
        idx_prey = argmax(spectrum_prey[2:end]) + 1
        period = 1/freq_prey[idx_prey]
        push!(periods, (набор=idx, период=period))
    end
end

period_df = DataFrame(periods)
plt_periods = bar(period_df.набор, period_df.период,
                  xlabel="Набор параметров",
                  ylabel="Период колебаний",
                  title="Периоды колебаний для разных наборов",
                  label="Период",
                  color=:steelblue,
                  size=(800, 400))
savefig(plt_periods, plotsdir(script_name, "lv_periods.png"))

# График 5: Сравнение амплитуд колебаний
amplitudes = []
for (idx, df) in enumerate(results)
    prey_amp = (maximum(df.prey) - minimum(df.prey))/2
    predator_amp = (maximum(df.predator) - minimum(df.predator))/2
    push!(amplitudes, (набор=idx, prey_amp=prey_amp, predator_amp=predator_amp))
end

amp_df = DataFrame(amplitudes)
plt_amplitudes = groupedbar(amp_df.набор, [amp_df.prey_amp amp_df.predator_amp],
                            bar_position=:dodge,
                            label=["Амплитуда жертв" "Амплитуда хищников"],
                            xlabel="Набор параметров",
                            ylabel="Амплитуда",
                            title="Сравнение амплитуд колебаний",
                            size=(900, 500))
savefig(plt_amplitudes, plotsdir(script_name, "lv_amplitudes.png"))

# СОЗДАЕМ СВОДНУЮ СТАТИСТИКУ
summary_stats = DataFrame(
    набор = 1:length(parameter_sets),
    α = [p[1] for p in parameter_sets],
    β = [p[2] for p in parameter_sets],
    δ = [p[3] for p in parameter_sets],
    γ = [p[4] for p in parameter_sets],
    x_стационарное = [p[4]/p[3] for p in parameter_sets],
    y_стационарное = [p[1]/p[2] for p in parameter_sets],
    средние_жертвы = [mean(df.prey) for df in results],
    средние_хищники = [mean(df.predator) for df in results],
    макс_жертвы = [maximum(df.prey) for df in results],
    макс_хищники = [maximum(df.predator) for df in results],
    мин_жертвы = [minimum(df.prey) for df in results],
    мин_хищники = [minimum(df.predator) for df in results]
)

# Добавляем периоды в сводную статистику
period_dict = Dict(row.набор => row.период for row in eachrow(period_df))
summary_stats[!, :период] = [get(period_dict, i, NaN) for i in 1:nrow(summary_stats)]

# Сохраняем сводную статистику
using CSV
CSV.write(datadir(script_name, "lv_summary_stats.csv"), summary_stats)

println("\n" * "="^60)
println("СВОДНАЯ СТАТИСТИКА")
println("="^60)
println(summary_stats)

# График 6: Тепловая карта зависимости от α и γ
plt_heatmap_αγ = scatter([p[1] for p in parameter_sets], 
                         [p[4] for p in parameter_sets],
                         zcolor=[maximum(df.prey) for df in results],
                         xlabel="α (скорость размножения жертв)",
                         ylabel="γ (смертность хищников)",
                         title="Максимальная популяция жертв",
                         label="Max жертвы",
                         size=(600, 400))
savefig(plt_heatmap_αγ, plotsdir(script_name, "lv_heatmap_alphagamma.png"))

# График 7: Для базового набора создаем детальную панель
base_idx = 1  # Классический набор
df_base = results[base_idx]

plt_panel = plot(layout=(3, 2), size=(1200, 900))
plot!(plt_panel[1], df_base.t, df_base.prey, label=L"x(t)", color=:green, linewidth=2, 
      title="Популяция жертв (Набор $base_idx)", grid=true)
plot!(plt_panel[2], df_base.t, df_base.predator, label=L"y(t)", color=:red, linewidth=2, 
      title="Популяция хищников (Набор $base_idx)", grid=true)
plot!(plt_panel[3], df_base.prey, df_base.predator, label=false, color=:blue, linewidth=1.5, 
      title="Фазовый портрет", xlabel=L"x", ylabel=L"y", grid=true)
scatter!(plt_panel[3], [df_base.x_star[1]], [df_base.y_star[1]], 
         color=:black, markersize=5, label="Равновесие")
plot!(plt_panel[4], df_base.t, [df_base.dprey_dt df_base.dpredator_dt], 
      label=[L"dx/dt" L"dy/dt"], color=[:green :red], linewidth=1.5, 
      title="Скорости изменения", grid=true, legend=:topright)
freq_base, spec_base = compute_fft(df_base.prey .- mean(df_base.prey), dt_lv)
plot!(plt_panel[5], freq_base, spec_base, label=L"Спектр x", color=:green, linewidth=1.5, 
      title="Спектр жертв", xscale=:log10, yscale=:log10, grid=true)
plot!(plt_panel[6], df_base.t, [df_base.prey_pct_change df_base.predator_pct_change], 
      label=[L"dx/x" L"dy/y"], color=[:green :red], linewidth=1.5, 
      title="Относительные изменения", grid=true, legend=:topright)
savefig(plt_panel, plotsdir(script_name, "lv_base_panel.png"))

println("\n" * "="^60)
println("Моделирование завершено успешно!")

