using DrWatson
@quickactivate "project"
using Agents, DataFrames, Plots, CSV, Random

include(srcdir("sir_model.jl"))

function create_migration_matrix(C, intensity)
	M = ones(C, C) .* intensity ./ (C-1)
	for i = 1:C
		M[i, i] = 1 - intensity
	end
	return M
end

function peak_time(p)
	migration_rates = create_migration_matrix(p[:C], p[:migration_intensity])
	model = initialize_sir(;
		Ns = p[:Ns],
		β_und = p[:β_und],
		β_det = p[:β_det],
		infection_period = p[:infection_period],
		detection_time = p[:detection_time],
		death_rate = p[:death_rate],
		reinfection_probability = p[:reinfection_probability],
		Is = p[:Is],
		seed = p[:seed],
		migration_rates = migration_rates,
	)

	infected_frac(model) = count(a.status == :I for a in allagents(model)) / nagents(model)
	peak = 0.0
	peak_step = 0

	for step = 1:p[:n_steps]
		agent_ids = collect(allids(model))
		for id in agent_ids
			agent = try
				model[id]
			catch
				nothing
			end
			if agent !== nothing
				sir_agent_step!(agent, model)
			end
		end
		frac = infected_frac(model)
		if frac > peak
			peak = frac
			peak_step = step
		end
	end
	return (peak_time=peak_step, peak_value=peak)
end

migration_intensities = 0.0:0.1:0.8
seeds = [42, 43, 44]
params_list = []
for mig in migration_intensities
	for s in seeds
	push!(
		params_list,
		Dict(
			:migration_intensity => mig, # скаляр
			:C => 3,
			:Ns => [1000, 1000, 1000],
			:β_und => [0.5, 0.5, 0.5],
			:β_det => [0.05, 0.05, 0.05],
			:infection_period => 14,
			:detection_time => 7,
			:death_rate => 0.02,
			:reinfection_probability => 0.1,
			:Is => [1, 0, 0],
			:seed => s,
			:n_steps => 150,
			),
		)
	end
end

results = []
for params in params_list
	data = peak_time(params)
	push!(results, merge(params, Dict(pairs(data))))
	println("Завершён эксперимент с migration_intensity = $(params[:migration_intensity]), seed = $(params[:seed])", )
end



df = DataFrame(results)
CSV.write(datadir("migration_scan_all2.csv"), df)

using Statistics
grouped = combine(
	groupby(df, [:migration_intensity]),
	:peak_time => mean => :mean_peak_time,
	:peak_value => mean => :mean_peak_value,
)

plot(
	grouped.migration_intensity,
	grouped.mean_peak_time,
	marker = :circle,
	xlabel = "Интенсивность миграции",
	ylabel = "Время до пика (дни)",
	label = "Время пика",
	)
plot!(
	grouped.migration_intensity,
	grouped.mean_peak_value .* 3000,
	marker = :square,
	xlabel = "Интенсивность миграции",
	ylabel = "Численность в пике",
	label = "Пиковая заболеваемость",
	)
savefig(plotsdir("migration_effect2.png"))
println("Результаты сохранены в data/migration_scan_all2.csv и plots/migration_effect2.png")


function run_quarantine_scenario(quarantine_threshold)
    M = create_migration_matrix(3, 0.3)

    model = initialize_sir(;
        Ns=[1000, 1000, 1000],
        β_und=[0.5, 0.5, 0.5],
        β_det=[0.05, 0.05, 0.05],
        infection_period=14,
        detection_time=7,
        death_rate=0.02,
        reinfection_probability=0.1,
        Is=[1, 0, 0],  # Инфекция начинается в городе 1
        seed=42,
        migration_rates=M,
        n_steps=150
    )

    peak_infected = 0.0
    quarantined = false

    for step in 1:150

        if !quarantined && quarantine_threshold > 0

            city1_infected = count(a -> a.status == :I && a.pos == 1, allagents(model))
            city1_pop = 1000
            if city1_infected / city1_pop > quarantine_threshold

                M[1, :] .= 0
                M[:, 1] .= 0
                model.migration_rates = M
                quarantined = true
                println("  Город 1 закрыт на карантин на шаге $step")
            end
        end

        Agents.step!(model, 1)
        current_peak = infected_count(model) / 3000
        if current_peak > peak_infected
            peak_infected = current_peak
        end
    end

    total_deaths = 3000 - nagents(model)
    return (peak=peak_infected, deaths=total_deaths)
end

thresholds = [0.0, 0.05, 0.10, 0.20]  # 0 = без карантина
scenario_names = ["Без карантина", "Ранний (5%)", "Умеренный (10%)", "Поздний (20%)"]

results_q = []
for (i, thresh) in enumerate(thresholds)
    println("\nТестирование: $(scenario_names[i])")
    result = run_quarantine_scenario(thresh)
    push!(results_q, (scenario=scenario_names[i],
                      peak=round(result.peak*100, digits=1),
                      deaths=result.deaths))
end

df_q = DataFrame(results_q)
println("\n=== Эффективность карантина ===")
println(df_q)

baseline_deaths = df_q.deaths[1]
for row in eachrow(df_q)
    if row.scenario != "Без карантина"
        reduction = (baseline_deaths - row.deaths) / baseline_deaths * 100
        println("$(row.scenario): $(round(reduction, digits=1))% снижение смертности")
    end
end

plot(thresholds[2:end], df_q.peak[2:end], marker=:circle,
     label="Пик заболеваемости", xlabel="Порог карантина",
     ylabel="Пик заболеваемости (%)", linewidth=2)
plot!(twinx(), thresholds[2:end], df_q.deaths[2:end], marker=:square,
      label="Общее число умерших", ylabel="Число умерших", legend=:topright)
title!("Эффективность карантина для города 1")
savefig(plotsdir("quarantine_effect.png"))
println("\nГрафик карантина сохранен в plots/quarantine_effect.png")
