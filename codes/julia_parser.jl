using CSV
using DataFrames


function read_csv_as_matrix(file_path::String)
    # Read the CSV file into a DataFrame
    df = CSV.read(file_path, DataFrame)
    
    # Convert the DataFrame to a matrix
    matrix = Matrix(df)
    
    # Return the matrix
    return matrix
end