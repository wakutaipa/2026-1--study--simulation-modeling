using DrWatson
using Random
using DataFrames
using Plots
using OrdinaryDiffEq
@quickactivate "Project"
include(srcdir("SIRPetri.jl"))
using .SIRPetri
using DataFrames, CSV, Plots

β = 0.3
γ = 0.1
tmax = 100.0
step = 0.2  

net, u0, states = build_sir_network(β, γ)

df = simulate_deterministic(net, u0, (0.0, tmax), saveat = step, rates = [β, γ])

@info "Creating animation... This may take a moment."
anim = @animate for t in eachrow(df)
    S_val = t.S
    I_val = t.I
    R_val = t.R
    time_val = round(t.time, digits=1)
    
    bar(
        ["S", "I", "R"],
        [S_val, I_val, R_val],
        ylabel = "Population",
        title = "SIR Dynamics at time = $(time_val)",
        color = [:green, :red, :blue],
        legend = false,
        ylims = (0, sum(u0)),
        bar_width = 0.6,
        label = ["S" "I" "R"]
    )
end

gif_path = plotsdir("sir_animation.gif")
gif(anim, gif_path, fps = 15)
println("Animation saved to: $(gif_path)")β = 0.3
γ = 0.1
tmax = 100.0
step = 0.2  

net, u0, states = build_sir_network(β, γ)

df = simulate_deterministic(net, u0, (0.0, tmax), saveat = step, rates = [β, γ])

@info "Creating animation... This may take a moment."
anim = @animate for t in eachrow(df)
    S_val = t.S
    I_val = t.I
    R_val = t.R
    time_val = round(t.time, digits=1)
    
    bar(
        ["S", "I", "R"],
        [S_val, I_val, R_val],
        ylabel = "Population",
        title = "SIR Dynamics at time = $(time_val)",
        color = [:green, :red, :blue],
        legend = false,
        ylims = (0, sum(u0)),
        bar_width = 0.6,
        label = ["S" "I" "R"]
    )
end

gif_path = plotsdir("sir_animation.gif")
gif(anim, gif_path, fps = 15)
println("Animation saved to: $(gif_path)")
