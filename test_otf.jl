#  Test file for the OTF decoding process. We attempt to load a pcm and then run OTF to it. 

using CodingTheory
using CodingTheory: CTMatrixTypes
using CodingTheory: AbstractClassicalNoiseChannel
using Distributions
# Import the function from the specified path
include("codes/julia_parser.jl")
include("src/Quantum/decoders/OTF.jl")

# Now you can call the function directly if it is in the global scope
# Assuming the function `read_csv_as_matrix` is defined at the top level in `julia_parser.jl`
d = 3
p = .1
# Call the function with a CSV file path
file_path = "codes/pcm_matrices/distance_$(d)_surface_code.csv"
H = read_csv_as_matrix(file_path)
channel_probs = fill(p, size(H,2))



# TODO, generate depolarizing channel
error = rand(Bernoulli(p), size(H,2)) .== 1



syndrome = (H * error).%2
# syndrome = rand(Bool, num_columns)

recovered_error = ordered_Tanner_forest(H, syndrome, channel_probs, :SP)

# Print or use the matrix

println(error)
println(syndrome)