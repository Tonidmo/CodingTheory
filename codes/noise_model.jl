using Random
"""
In this function we generate an error in the binary symplectic form (z|x) for convenience, so we can do the pcm, error product as a binary symplectic product.
"""
function Pauli_noise_channel(matrix::AbstractMatrix, p0::Float64, px::Float64, py::Float64, pz::Float64)
    # Check the conditions on the float values
    # @assert p0 > (px + py + pz) "v1 must be larger than the sum of v2, v3, and v4"
    @assert p0+px+py+pz == 1.0 "The sum of v1, v2, v3, and v4 must be equal to 1"
    
    # Generate a binary array of size equal to the number of columns in the matrix
    num_cols = size(matrix, 2)
    binary_array = falses(num_cols) .|> Int
    
    # Number of qubits (half the number of columns)
    num_qubits = div(num_cols, 2)
    
    # Iterate through each element to assign events
    for i in 1:num_qubits
        # Choose an event based on probabilities p0, p1, p2, p3
        event = rand(Float64)
        
        if event < p0
            # Event 0: Nothing happens
            continue
        elseif event < p0 + px
            # Event 1: X ERROR, Set binary_array[i+ num_qubits] to true
            binary_array[i + num_qubits] = true
        elseif event < p0 + px + py
            # Event 2: Y ERROR, Set both binary_array[i] and binary_array[i + num_realizations] to true
            binary_array[i] = true
            binary_array[i + num_qubits] = true
        else
            # Event 3: Z ERROR Set binary_array[i] to true
            binary_array[i] = true
        end
    end
    return binary_array
end


"""
In this function we generate two matrices, one containing a number of errors and another containing a number of their correspondent syndromes for a given QECC.
The arguments are the number of iterations, how many errors and syndromes are generated, the parity check matrix of the code and the Pauli probabilities.
"""
function noiseless_syndrome_extractions(
    number_of_iterations::Int64, 
    matrix::AbstractMatrix, 
    p0::Float64, 
    px::Float64, 
    py::Float64, 
    pz::Float64)

    # We allocate memory for the desired matrices
    error_matrix = falses( number_of_iterations, size(matrix,2))
    syndrome_matrix = falses( number_of_iterations, size(matrix,1))

    # For every iteration we generate an error and compute its syndrome.
    for i in 1:number_of_iterations
        error =  Pauli_noise_channel(matrix, p0, px, py, pz)
        syndrome = matrix * error .%2
        error_matrix[i,:] = error
        syndrome_matrix[i,:] = syndrome
    end

    return error_matrix, syndrome_matrix

end



function assign_noiseless_weights(
    matrix::AbstractMatrix, 
    px::Float64, 
    py::Float64, 
    pz::Float64)

    priors = zeros(size(matrix, 2))


    num_qubits = div(size(matrix, 2), 2)

    first_columns = pz + py - (pz*py)
    second_columns = px + py - (px*py)

    for i in 1:num_qubits

        priors[i] = log((1-first_columns)/first_columns)
        priors[i+num_qubits] = log((1-second_columns)/second_columns)

    end

    return priors
end


function adapt_pcm(pcm::BitMatrix)
    new_row = falses(1, size(pcm, 2))
    add_row = false
    for col in 1:size(pcm, 2)
        # Check if there is exactly one `true` element in the column
        if count(pcm[:, col]) == 1
            new_row[1, col] = true
            add_row = true
        end
    end
    pcm_otf = vcat(pcm, new_row)
    return pcm_otf
end