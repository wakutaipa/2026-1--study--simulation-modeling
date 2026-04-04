using DrWatson
@quickactivate "project"

using BlackBoxOptim, Random, Statistics

include(srcdir("sir_model.jl"))

function cost_multi(x)
	model = initialize_sir(;
		Ns = [1000, 1000, 1000],
		β_und = fill(x[1], 3),
		β_det = fill(x[1]/10, 3),
		infection_period = 14,
		detection_time = round(Int, x[2]),
		death_rate = x[3],
		reinfection_probability = 0.1,
		Is = [0, 0, 1],
		seed = 42,
		n_steps = 100,
	)
	infected_frac(model) = count(a.status == :I for a in allagents(model)) / nagents(model)
	dead_count(model) = 3000 - nagents(model)
	peak_infected = 0.0
	replicates = 5
	peak_vals = Float64[]
	dead_vals = Int[]
	for rep = 1:replicates
		model = initialize_sir(;
		Ns = [1000, 1000, 1000],
		β_und = fill(x[1], 3),
		β_det = fill(x[1]/10, 3),
		infection_period = 14,
		detection_time = round(Int, x[2]),
		death_rate = x[3],
		reinfection_probability = 0.1,
		Is = [0, 0, 1],
		seed = 42 + rep,
		n_steps = 100,
		)
		for step = 1:100
			Agents.step!(model, 1)
			frac = infected_frac(model)
			if frac > peak_infected
				peak_infected = frac
			end
		end
		push!(peak_vals, peak_infected)
		push!(dead_vals, dead_count(model))
	end
	return (mean(peak_vals), mean(dead_vals) / 3000)
end

# Запуск оптимизации
result = bboptimize(
	cost_multi,
	Method = :borg_moea,
	FitnessScheme = ParetoFitnessScheme{2}(is_minimizing = true),
	SearchRange = [
		(0.1, 1.0),
		(3.0, 14.0),
		(0.01, 0.1),
	],
	NumDimensions = 3,
	MaxTime = 120, # 2 минуты
	TraceMode = :compact,
	)
	
best = best_candidate(result)
fitness = best_fitness(result)
println("Оптимальные параметры:")
println("β_und = $(best[1])")
println("Время выявления = $(round(Int, best[2])) дней")
println("Смертность = $(best[3])")
println("Достигнутые показатели:")
println("Пик заболеваемости: $(fitness[1])")

println("Доля умерших: $(fitness[2])")
save(datadir("optimization_result.jld2"), Dict("best" => best, "fitness" => fitness))

# Add this to sir_optimize_parameters.jl

# Task 6: Constrained optimization (peak < 30%)
function cost_multi_constrained(x)
    β_und = x[1]
    detection_time = round(Int, x[2])
    death_rate = x[3]
    
    replicates = 5
    peak_vals = Float64[]
    dead_vals = Int[]
    
    for rep in 1:replicates
        model = initialize_sir(;
            Ns = [1000, 1000, 1000],
            β_und = fill(β_und, 3),
            β_det = fill(β_und/10, 3),
            infection_period = 14,
            detection_time = detection_time,
            death_rate = death_rate,
            reinfection_probability = 0.1,
            Is = [0, 0, 1],
            seed = 42 + rep,
            n_steps = 100,
        )
        
        peak = 0.0
        for step in 1:100
            Agents.step!(model, 1)
            frac = infected_count(model) / 3000
            if frac > peak
                peak = frac
            end
        end
        
        deaths = 3000 - nagents(model)
        
        # Penalty if peak exceeds 30%
        if peak > 0.3
            peak = peak + 10.0  # Heavy penalty
        end
        
        push!(peak_vals, peak)
        push!(dead_vals, deaths)
    end
    
    return (mean(dead_vals), mean(peak_vals))
end

result_constrained = bboptimize(
    cost_multi_constrained,
    Method = :borg_moea,
    FitnessScheme = ParetoFitnessScheme{2}(is_minimizing = true),
    SearchRange = [
        (0.1, 1.0),     
        (3.0, 14.0),   
        (0.01, 0.1),    
    ],
    NumDimensions = 3,
    MaxTime = 120,
    TraceMode = :compact,
)

best_constrained = best_candidate(result_constrained)
fitness_constrained = best_fitness(result_constrained)

println("\nОптимальные параметры (с ограничением пик < 30%):")
println("β_und = $(round(best_constrained[1], digits=3))")
println("Время выявления = $(round(Int, best_constrained[2])) дней")
println("Коэффициент смертности = $(round(best_constrained[3], digits=4))")
println("\nДостигнутые показатели:")
println("Среднее число умерших = $(round(fitness_constrained[1]))")
println("Средний пик заболеваемости = $(round(fitness_constrained[2]*100, digits=1))%")

# Save results
save(datadir("optimization_constrained.jld2"), 
     Dict("best" => best_constrained, "fitness" => fitness_constrained))

