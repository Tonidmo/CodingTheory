function _message_passing_filling(
    H::Matrix{T}, 
    syndrome::Union{Missing, Vector{T}},
    chn_inits::Vector{Float64}, 
    c_to_v_mess::Function, 
    var_adj_list::Vector{Vector{Int}},
    check_adj_list::Vector{Vector{Int}}, 
    max_iter::Int,
    current_bits::Vector{Int}, 
    totals::Array{Float64, 2}, 
    syn::Vector{Int},
    check_to_var_messages::Array{Float64, 3}, 
    var_to_check_messages::Array{Float64, 3},
    attenuation::Float64, 
    layers::Vector{Vector{Int}},
    antilayers::Vector{Vector{Int}},
    type_B::Vector{Bool}
    ) where T <: Integer

    # first iteration for variable nodes - set to channel initialization
    num_check, num_var = size(H)


    @inbounds for v in 1:num_var
        @simd for c in var_adj_list[v]
            var_to_check_messages[v, c, 2] = chn_inits[v]
        end
    end


    iter::Int = 1
    curr_iter::Int = 2
    prev_iter::Int = 1
    


    # FILLING PHASE, each layer contains a number of variable nodes that we will
    #= 
        CUP BEING FILLED
        |
        |     ← Syndromes stream in
        v
    
    |       |
    |       |
    |       |     ← Syndromes rising
    |       |
    \_______/
        =#
    for (_ , layer) in enumerate(layers)
        for _ in 1:1
            @simd for c in layer
                for v in check_adj_list[c]
                    if !ismissing(syndrome)
                        check_to_var_messages[c, v, curr_iter] = (-1)^syndrome[c] * c_to_v_mess(c, v,
                        curr_iter, check_adj_list, var_to_check_messages, attenuation)
                    else
                        check_to_var_messages[c, v, curr_iter] = c_to_v_mess(c, v,
                        curr_iter, check_adj_list, var_to_check_messages, attenuation)
                    end
                end
            end

            # TODO the only values that should be changing here are the ones connected to the check nodes in the layer
            @simd for v in 1:num_var
                totals[v, curr_iter] = chn_inits[v]
                for c in var_adj_list[v]
                    totals[v, curr_iter] += check_to_var_messages[c, v, curr_iter]
                end
                if sign(totals[v,curr_iter]) != sign(totals[v,prev_iter]) && type_B[v]
                    # if the sign changed, we need to update the previous value
                    totals[v, curr_iter] += totals[v, prev_iter]
                end 
                current_bits[v] = totals[v,curr_iter] >= 0 ? 0 : 1
            end


            @simd for v in 1:num_var
                for c in var_adj_list[v]
                    # this includes the channel inputs in total
                    # TODO: here prev was changed to curr
                    var_to_check_messages[v, c, curr_iter] = totals[v,curr_iter] -
                        check_to_var_messages[c, v, curr_iter]
                end
            end

            # switch these current values to previous values so the next layer can use them
            #TODO Check what is going on here, is this wrong?
            # @inbounds @simd for c in layer
            #     for v in check_adj_list[c]
            #         check_to_var_messages[c, v, prev_iter] = check_to_var_messages[c, v, curr_iter]
            #     end
            # end


            temp = curr_iter
            curr_iter = prev_iter
            prev_iter = temp
        end
    end


        

    # Tanner graph has been filled, now we need to check the syndrome

        #= 
            FILLED CUP
            
                 ← Syndromes stop streaming in
            
        
        |~~~~~~~|
        |       |
        |       |     ← d syndromes
        |       |
        \_______/
        =#


    LinearAlgebra.mul!(syn, H, current_bits)
    if !ismissing(syndrome)
        all(syn[i] % 2 == syndrome[i] for i in 1:num_check) && return true, current_bits, iter, totals[:, prev_iter]
    else
        all(iszero(syn[i] % 2) for i in 1:num_check) && return true, current_bits, iter, totals[:,  prev_iter]
    end

    # Even though the Tanner graph has been filled, bp does not converge, we need to proceed with parallel iterations.

    for _ in 1:max_iter
        @simd for c in 1:num_check
            for v in check_adj_list[c]
                if !ismissing(syndrome)
                    check_to_var_messages[c, v, curr_iter] = (-1)^syndrome[c] * c_to_v_mess(c, v,
                    curr_iter, check_adj_list, var_to_check_messages, attenuation)
                else
                    check_to_var_messages[c, v, curr_iter] = c_to_v_mess(c, v,
                    curr_iter, check_adj_list, var_to_check_messages, attenuation)
                end
            end
        end

        @simd for v in 1:num_var
            totals[v, curr_iter] = chn_inits[v]
            for c in var_adj_list[v]
                totals[v,curr_iter] += check_to_var_messages[c, v, curr_iter]
            end
            if sign(totals[v,curr_iter]) != sign(totals[v,prev_iter]) && type_B[v]
                # if the sign changed, we need to update the previous value
                totals[v, curr_iter] += totals[v, prev_iter]
            end 
            current_bits[v] = totals[v,curr_iter] >= 0 ? 0 : 1
        end


        @simd for v in 1:num_var
            for c in var_adj_list[v]
                # this includes the channel inputs in total
                # TODO: here prev was changed to curr
                var_to_check_messages[v, c, curr_iter] = totals[v,curr_iter] -
                    check_to_var_messages[c, v, curr_iter]
            end
        end

        temp = curr_iter
        curr_iter = prev_iter
        prev_iter = temp
        iter += 1

        LinearAlgebra.mul!(syn, H, current_bits)
        if !ismissing(syndrome)
            all(syn[i] % 2 == syndrome[i] for i in 1:num_check) && return true, current_bits, iter, totals[:, prev_iter]
        else
            all(iszero(syn[i] % 2) for i in 1:num_check) && return true, current_bits, iter, totals[:, prev_iter]
        end
    end


    # totals[:, curr_iter] = copy(totals[:, prev_iter])

    # Emptying Phase, as the process continues, syndromes are continuously being processed.

        #= 
            EMPTYING CUP
            
        
        |       |
        |       |
        |~~~~~~~|     ← < d syndromes
        |       |
        \_______/

            |
            |       Processed syndromes stream out
            |
            v
        =#

    #TODO How do we store the 
    println("Onto emptying stage")
    for (_ , layer) in enumerate(antilayers)
        for _ in 1:1
            @simd for c in layer
                for v in check_adj_list[c]
                    if !ismissing(syndrome)
                        check_to_var_messages[c, v, curr_iter] = (-1)^syndrome[c] * c_to_v_mess(c, v,
                        curr_iter, check_adj_list, var_to_check_messages, attenuation)
                    else
                        check_to_var_messages[c, v, curr_iter] = c_to_v_mess(c, v,
                        curr_iter, check_adj_list, var_to_check_messages, attenuation)
                    end
                end
            end

            @simd for v in 1:num_var
                # TODO for variable nodes not experiencing additional messages, the following line should be:
                # totals[v, curr_iter] += check_to_var_messages[c, v, prev_iter]
                # I do not know how to implement that efficiently
                for c in var_adj_list[v]
                    totals[v, curr_iter] += check_to_var_messages[c, v, curr_iter]
                end

                if sign(totals[v,curr_iter]) != sign(totals[v,prev_iter]) && type_B[v]
                    # if the sign changed, we need to update the previous value
                    totals[v, curr_iter] += totals[v, prev_iter]
                end
                current_bits[v] = totals[v,curr_iter] >= 0 ? 0 : 1
            end


            @simd for v in 1:num_var
                for c in var_adj_list[v]
                    # this includes the channel inputs in total
                    # TODO: here prev was changed to curr
                    var_to_check_messages[v, c, curr_iter] = totals[v,curr_iter] -
                        check_to_var_messages[c, v, curr_iter]
                end
            end

            # switch these current values to previous values so the next layer can use them
            # @inbounds @simd for c in layer
            #     for v in check_adj_list[c]
            #         check_to_var_messages[c, v, prev_iter] = check_to_var_messages[c, v, curr_iter]
            #     end
            # end


            temp = curr_iter
            curr_iter = prev_iter
            prev_iter = temp

        end
    end

    # Tanner graph has been emptied again, let's check if syndrome converged
    LinearAlgebra.mul!(syn, H, current_bits)
    if !ismissing(syndrome)
        all(syn[i] % 2 == syndrome[i] for i in 1:num_check) && return true, current_bits, iter, totals[:, prev_iter]
    else
        all(iszero(syn[i] % 2) for i in 1:num_check) && return true, current_bits, iter, totals[:,  prev_iter]
    end

    log_likelihood_ratios = copy(totals[:, prev_iter])

    return false, current_bits, iter, log_likelihood_ratios
end