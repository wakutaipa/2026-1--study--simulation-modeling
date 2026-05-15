using DrWatson
@quickactivate "project"
using StableRNGs
using Distributions
using ConcurrentSim
using ResumableFunctions
using Plots

rng = StableRNG(123)
num_customers = 10  # total number of customers generated

num_servers = 2  # number of servers
mu = 1.0 / 2  # service rate
lam = 0.9  # arrival rate
arrival_dist = Exponential(1 / lam)  # interarrival time distribution
service_dist = Exponential(1 / mu)  # service time distribution

mutable struct SimData
    arrival_times::Vector{Float64}
    service_start_times::Vector{Float64}
    service_end_times::Vector{Float64}
    customer_ids::Vector{Int}
end

sim_data = SimData([], [], [], [])

@resumable function customer(
    env::Environment,
    server::Resource,
    id::Integer,
    t_a::Float64,
    d_s::Distribution,
    data::SimData
)
    @yield timeout(env, t_a)  # customer arrives
    println("Customer $id arrived: ", now(env))
    push!(data.arrival_times, now(env))  # ← Record arrival HERE
    push!(data.customer_ids, id)         # ← Move ID recording here too

    @yield request(server)  # customer starts service
    println("Customer $id entered service: ", now(env))
    push!(data.service_start_times, now(env))  # ← Record start HERE

    @yield timeout(env, rand(rng, d_s))  # server is busy
    @yield release(server)  # customer exits service
    println("Customer $id exited service: ", now(env))
    push!(data.service_end_times, now(env))  # ← Record end HERE
end

function setup_and_run()
    sim = Simulation()  # initialize simulation environment
    server = Resource(sim, num_servers)  # initialize servers
    arrival_time = 0.0

    for i in 1:num_customers  # initialize customers
        arrival_time += rand(rng, arrival_dist)
        @process customer(sim, server, i, arrival_time, service_dist, sim_data)
    end

    run(sim)  # run simulation
end

setup_and_run()

function make_graphs(data::SimData)

    p1 = plot(data.arrival_times, data.customer_ids, label="Arrivals", color=:red, marker=:circle)
    xlabel!(p1, "Time")
    ylabel!(p1, "Customer ID")
    title!(p1, "Customer Arrival Times")
    savefig(p1, plotsdir("customer_arrival_times.png"))
    display(p1)

    p2 = plot(data.service_start_times, data.customer_ids, label="Service Start", color=:blue, marker=:square)
    xlabel!(p2, "Time")
    ylabel!(p2, "Customer ID")
    title!(p2, "Customer Service Start Times")
    savefig(p2, plotsdir("customer_service_start_times.png"))
    display(p2)

    p3 = plot(data.service_end_times, data.customer_ids, label="Service End", color=:green, marker=:diamond)
    xlabel!(p3, "Time")
    ylabel!(p3, "Customer ID")
    title!(p3, "Customer Service End Times")
    savefig(p3, plotsdir("customer_service_end_times.png"))
    display(p3)

    wait = data.service_end_times .- data.service_start_times
    p4 = plot(wait, data.customer_ids, label="Service End", color=:green, marker=:diamond)
    xlabel!(p4, "wait")
    ylabel!(p4, "Customer ID")
    title!(p4, "Customer Wait Times")
    savefig(p4, plotsdir("customer_wait_times.png"))
    display(p4)

    println("All plots saved in plots directory.")
end

make_graphs(sim_data)
