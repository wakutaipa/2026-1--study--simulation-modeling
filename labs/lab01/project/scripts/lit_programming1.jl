using Literate

input_file = joinpath("scripts", "01_exponential_growth.jl")

Literate.notebook(input_file, "notebooks")
Literate.markdown(input_file, "papers")
