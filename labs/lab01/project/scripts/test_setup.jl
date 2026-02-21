##!/usr/bin/env julia
## test_setup.jl

using DrWatson
@quickactivate "project"

println(" Project activated: ", projectdir())

packages = [
	"DrWatson",
	"DifferentialEquations",
	"Plots",
	"DataFrames",
	"CSV",
	"JLD2",
	"Literate",
	"IJulia",
	"BenchmarkTools",
	"Quarto"
]

println("\nTesting packages: ")
for pkg in packages
	try
		eval(Meta.parse("using $pkg"))
		println("  ✓ $pkg")
	catch e
		println("  ✗ $pkg: error loading")
		showerror(stdout, e)
	end
end

println("\nProject structure:")
println(" Root:             ", projectdir())
println(" Data:             ", datadir())
println(" Scripts:          ", srcdir())
println(" Graphs:           ", plotsdir())
