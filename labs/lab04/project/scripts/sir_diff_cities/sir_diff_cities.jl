using DrWatson
@quickactivate "project"
using Agents, DataFrames, Plots, CSV, Statistics

include(srcdir("sir_model.jl"))

println(repeat("=", 60))
println("ЗАДАНИЕ 3: Эффект гетерогенности")
println(repeat("=", 60))

function get_infections_by_city(model, Ns)
    city_infections = zeros(Int, length(Ns))
    for agent in allagents(model)
        if agent.status == :I
            if agent.id <= Ns[1]
                city_infections[1] += 1
            elseif agent.id <= Ns[1] + Ns[2]
                city_infections[2] += 1
            else
                city_infections[3] += 1
            end
        end
    end
    return city_infections
end

function run_simulation(β_und_values; seed=42, n_steps=100)
    println("Запуск симуляции с β = $(β_und_values)...")

    model = initialize_sir(;
        Ns = [1000, 1000, 1000],
        β_und = β_und_values,
        β_det = β_und_values ./ 10,
        infection_period = 14,
        detection_time = 7,
        death_rate = 0.02,
        reinfection_probability = 0.1,
        Is = [0, 0, 1],
        seed = seed,
        n_steps = n_steps,
    )

    city_infections = [Int[] for _ in 1:3]
    times = 1:n_steps

    for step in times
        Agents.step!(model, 1)
        city_inf = get_infections_by_city(model, [1000, 1000, 1000])
        for city in 1:3
            push!(city_infections[city], city_inf[city])
        end
    end

    return times, city_infections
end

println("\n1. Запуск однородного сценария (β = 0.5 для всех городов)...")
β_homog = [0.5, 0.5, 0.5]
times, infections_homog = run_simulation(β_homog, seed=42)

println("\n2. Запуск гетерогенного сценария (β = [0.03, 0.6, 0.4])...")
β_heterog = [0.03, 0.6, 0.4]
times, infections_heterog = run_simulation(β_heterog, seed=42)

total_pop = 3000

peak_homog = maximum([maximum(infections_homog[i]) for i in 1:3])
final_infected_homog = sum([infections_homog[i][end] for i in 1:3])

peak_heterog = maximum([maximum(infections_heterog[i]) for i in 1:3])
final_infected_heterog = sum([infections_heterog[i][end] for i in 1:3])

println("\n=== Результаты ===")
println("\nОднородный сценарий (β = [0.5, 0.5, 0.5]):")
println("  Пик заболеваемости: $(peak_homog) ($(round(peak_homog/total_pop*100, digits=1))%)")
println("  Конечная доля инфицированных: $(final_infected_homog) ($(round(final_infected_homog/total_pop*100, digits=1))%)")

println("\nГетерогенный сценарий (β = [0.03, 0.6, 0.4]):")
println("  Пик заболеваемости: $(peak_heterog) ($(round(peak_heterog/total_pop*100, digits=1))%)")
println("  Конечная доля инфицированных: $(final_infected_heterog) ($(round(final_infected_heterog/total_pop*100, digits=1))%)")
println("  Город A, пик: $(maximum(infections_heterog[1])) ($(round(maximum(infections_heterog[1])/1000*100, digits=1))%)")
println("  Город B, пик: $(maximum(infections_heterog[2])) ($(round(maximum(infections_heterog[2])/1000*100, digits=1))%)")
println("  Город C, пик: $(maximum(infections_heterog[3])) ($(round(maximum(infections_heterog[3])/1000*100, digits=1))%)")

p1 = plot(times, infections_homog[1] ./ 1000, label="Город A (β=0.5)", linewidth=2,
          title="Гомогенная передача - Город A", xlabel="Дни", ylabel="Доля инфицированных")

p2 = plot(times, infections_homog[2] ./ 1000, label="Город B (β=0.5)", linewidth=2,
          title="Гомогенная передача - Город B", xlabel="Дни", ylabel="Доля инфицированных")

p3 = plot(times, infections_homog[3] ./ 1000, label="Город C (β=0.5)", linewidth=2,
          title="Гомогенная передача - Город C", xlabel="Дни", ylabel="Доля инфицированных")

p4 = plot(times, infections_heterog[1] ./ 1000, label="Город A (β=0.03, высокий риск)",
          linewidth=2, title="Гетерогенная передача - Город A", xlabel="Дни", ylabel="Доля инфицированных")

p5 = plot(times, infections_heterog[2] ./ 1000, label="Город B (β=0.6)", linewidth=2,
          title="Гетерогенная передача - Город B", xlabel="Дни", ylabel="Доля инфицированных")

p6 = plot(times, infections_heterog[3] ./ 1000, label="Город C (β=0.4)", linewidth=2,
          title="Гетерогенная передача - Город C", xlabel="Дни", ylabel="Доля инфицированных")

plot(p1, p2, p3, layout=(3, 1), size=(800, 900))
savefig(plotsdir("task3_homogeneity.png"))

plot(p4, p5, p6, layout=(3, 1), size=(800, 900))
savefig(plotsdir("task3_heterogeneous.png"))

println("РЕЗУЛЬТАТЫ СОХРАНЕНЫ")
