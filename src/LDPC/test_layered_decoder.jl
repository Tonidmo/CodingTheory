# test_layered_decoder.jl

# Load the package (assumes CodingTheory is a module with a Project.toml)
using Pkg
Pkg.activate(".")  # activate local environment
Pkg.instantiate()  # make sure dependencies are installed

# Load the module
using CodingTheory

# Now you can call the function (we'll just show how to set up dummy data for now)
H = [1 0 1; 0 1 1]              # Small example parity check matrix
syndrome = [0, 1]               # Dummy syndrome
chn_inits = [0.9, -1.2, 0.3]    # Dummy LLRs
current_bits = [0, 0, 0]
totals = zeros(Float64, 3)
syn = zeros(Int, 2)

# Dummy graph structure (usually constructed from H)
var_adj_list = [[1,3], [2,3], [1,2]]
check_adj_list = [[1,3], [2,3]]

# Dummy message arrays
check_to_var_messages = zeros(Float64, 2, 3, 2)
var_to_check_messages = zeros(Float64, 3, 2, 2)

# Layers (here we assume 2 layers of check nodes)
layers = [[1], [2]]

# Dummy message function
c_to_v_mess(c, v, curr_iter, check_adj_list, var_to_check_messages, attenuation) = 0.5

# Call the function
converged, decoded, iterations, totals = CodingTheory._message_passing_layered(
    H, syndrome, chn_inits, c_to_v_mess, var_adj_list, check_adj_list,
    10, :layered, current_bits, totals, syn, check_to_var_messages,
    var_to_check_messages, 0.75, layers
)

println("Converged: $converged")
println("Decoded bits: $decoded")
