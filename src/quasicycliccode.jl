# Copyright (c) 2022, 2023 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    QuasiCyclicCode(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}}, parity::Bool=false)

Return the quasi-cycle code specified by the matrix `A` of polynomial circulant generators. If the
optional paramater `parity` is set to `true`, the input is used to construct the parity-check matrix.
"""
function QuasiCyclicCode(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}}, parity::Bool=false)
    R = parent(A[1, 1])
    S = base_ring(R)
    F = base_ring(S)
    g = modulus(R)
    m = degree(g)
    g == gen(S)^m - 1 || throw(ArgumentError("Residue ring not of the form x^m - 1."))
    l = ncols(A)
    if parity
        Atype = 'H'
        H = lift(A)
        k, _ = right_kernel(H)
        W = weightmatrix(A)
        return QuasiCyclicCode(F, R, ncols(H), k, missing, 1, ncols(H), missing, missing, missing, missing, missing, missing, l, m, A, Atype, W, maximum(W))
    else
        Atype = 'G'
        G = lift(A)
        k = rank(G)
        W = weightmatrix(A)
        return QuasiCyclicCode(F, R, ncols(G), k, missing, 1, ncols(G), missing, missing, missing, missing, missing, missing, l, m, A, Atype, W, maximum(W))
    end
end

"""
    QuasiCyclicCode(v::Vector{fq_nmod_mat}, l::Int, circgens::Bool, parity::Bool=false)

Return the quasi-cyclic code of index `l` generated by right-bit shifts of size `l` of the
generator vectors `v`. If `circgens` is `true`, the vectors are taken to be (column) generators
for the circulant matrices instead of generator vectors for the code. If the optional paramater
`parity` is set to `true`, the input is used to construct the parity-check matrix.

Note:
- If `circgens` is `false`, then the length of the code is `ncols(v[1])` and must be divisible by `l`.
- If `circgens` is `true`, then the length of the code is `ncols(v[1]) * l`. Circulant matrices are
  stacked in rows of length `l`, so `l` must divide `length(v)`.
"""
function QuasiCyclicCode(v::Vector{fq_nmod_mat}, l::Int, circgens::Bool, parity::Bool=false)
    F = base_ring(v[1])
    lenv = length(v)
    if circgens
        lenv >= 2 || throw(ArgumentError("Length of input vector must be at least two."))
        lenv % l == 0 || throw(ArgumentError("The length of the input vector must be divisible by l."))
        nr = div(lenv, l)
        r, m = size(v[1])
        (r != 1 && m != 1) && throw(ArgumentError("The input matrices must be vectors."))
        m == 1 && (v[1] = transpose(v[1]); (r, m = size(v[1]));)
        for i in 2:lenv
            F == base_ring(v[i]) || throw(ArgumentError("All inputs must be over the same base ring."))
            r2, m2 = size(v[i])
            (r2 != 1 && m2 != 1) && throw(ArgumentError("The input matrices must be vectors."))
            m2 == 1 && (v[i] = transpose(v[i]); (r2, m2 = size(v[i]));)
            m == m2 || throw(ArgumentError("The input vectors must all be the same length."))
        end

        S, x = PolynomialRing(F, "x")
        R = ResidueRing(S, x^m - 1)
        A = zero_matrix(R, nr, l)
        for r in 1:nr
            for c in 1:l
                temp = [v[(r - 1) * l + c][i] for i in 1:m]
                A[r, c] = R(S(temp))
            end
        end
        return QuasiCyclicCode(A, parity)
    else
        r, n = size(v[1])
        (r != 1 && n != 1) && throw(ArgumentError("The input matrices must be vectors."))
        n == 1 && (v[1] = transpose(v[1]); (r, n = size(v[1]));)
        n % l == 0 || throw(ArgumentError("Parameter l must divide the length of the vector."))
        m = div(n, l)
        for i in 2:lenv
            F == base_ring(v[i]) || throw(ArgumentError("All vectors must be over the same base ring."))
            r2, n2 = size(v[i])
            (r2 != 1 && n2 != 1) && throw(ArgumentError("The input matrices must be vectors."))
            n2 == 1 && (v[i] = transpose(v[i]); (r2, n2 = size(v[i]));)
            n == n2 || throw(ArgumentError("The input vectors must all be the same length."))
        end
        
        S, x = PolynomialRing(F, "x")
        R = ResidueRing(S, x^m - 1)
        A = zero_matrix(R, lenv, l)
        for k in 1:lenv
            for i in 1:l
                row = v[k]
                topcircrow = zero_matrix(F, m, 1)
                for j in 1:m
                    topcircrow[j, 1] = row[i + (j - 1) * l]
                end
                # transpose to get first circulant column
                topcircrow[2:end, :] = topcircrow[end:-1:2, :]
                A[k, i] = R(S([topcircrow[i, 1] for i in 1:m]))
            end
        end
        return QuasiCyclicCode(A, parity)
    end
end

"""
    QuasiCyclicCode(v::fq_nmod_mat, l::Int, parity::Bool=false)

Return the quasi-cyclic code of index `l` generated by right-bit shifts of size `l` of the
generator vector `v`. If the optional paramater `parity` is set to `true`, the input is used
to construct the parity check matrix.
"""
QuasiCyclicCode(v::fq_nmod_mat, l::Int, parity::Bool=false) = QuasiCyclicCode([v], l, false, parity)

"""
    QuasiCyclicCode(v::Vector{fq_nmod_poly}, n::Int, l::Int, parity::Bool=false)

Return the quasi-cyclic code of index `l` whose circulants are defined by the generator
polynomials `v`. If the optional paramater `parity` is set to `true`, the input is used
to construct the parity check matrix.
"""
function QuasiCyclicCode(v::Vector{fq_nmod_poly}, n::Int, l::Int, parity::Bool=false)
    # if g = x^10 + α^2*x^9 + x^8 + α*x^7 + x^3 + α^2*x^2 + x + α
    # g.coeffs = [α  1  α^2  1  0  0  0  α  1  α^2  1]
    genvecs = Vector{fq_nmod_mat}()
    F = base_ring(v[1])
    for g in v
        temp = zero_matrix(F, 1, n)
        coeffs = collect(coefficients(g))
        temp[1, 1:length(coeffs)] = coeffs
        push!(genvecs, temp)
    end
    return QuasiCyclicCode(genvecs, l, true, parity)
end

"""
    QuasiCyclicCode(v::Vector{AbstractCyclicCode}, l::Int, parity::Bool=false)

Return the quasi-cyclic code of index `l` whose circulants are determined by the cyclic
codes in `v`. If the optional paramater `parity` is set to `true`, the input is used to
construct the parity check matrix.
"""
function QuasiCyclicCode(v::Vector{AbstractCyclicCode}, l::Int, parity::Bool=false)
    genvecs = Vector{fq_nmod_mat}()
    for C in v
        push!(genvecs, C.G[1, :])
    end
    return QuasiCyclicCode(gences, l, true, parity)
end

#############################
      # getter functions
#############################

"""
    index(C::AbstractQuasiCyclicCode)

Return the index of the quasi-cyclic code.
"""
index(C::AbstractQuasiCyclicCode) = C.l

"""
    expansionfactor(C::AbstractQuasiCyclicCode)

Return the expansion factor of the quasi-cycle code `C`.
"""
expansionfactor(C::AbstractQuasiCyclicCode) = C.m

"""
    type(C::AbstractQuasiCyclicCode)

Return the type of the quasi-cycle code `C`.
"""
type(C::AbstractQuasiCyclicCode) = C.type

# need these because C::AbstractLinearCode requires it
originalgeneratormatrix(C::AbstractQuasiCyclicCode) = missing

originalparitycheckmatrix(C::AbstractQuasiCyclicCode) = missing

"""
    polynomialmatrix(C::AbstractQuasiCyclicCode)

Return the polynomial matrix used to define the code.

Use `polynomialmatrixtype` to determine if specifies the generator or parity-check matrix.
"""
polynomialmatrix(C::AbstractQuasiCyclicCode) = C.A

"""
    polynomialmatrixtype(C::AbstractQuasiCyclicCode)

Return `'G'` if the polynomial matrix of `C` specifies the generator or parity-check matrix.
"""
polynomialmatrixtype(C::AbstractQuasiCyclicCode) = C.Atype

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

"""
    basematrix(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}})
    protographmatrix(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}})
    weightmatrix(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}})

Return the base/protograph/weight matrix of `A`.
"""
function weightmatrix(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}})
    nr, nc = size(A)
    W = zeros(Int, nr, nc)
    for c in 1:nc
        for r in 1:nr
            W[r, c] = wt(Nemo.lift(A[r, c]))
        end
    end
    return W
end
basematrix(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}}) = weightmatrix(A)
protographmatrix(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}}) = weightmatrix(A)

"""
    issinglegenerator(C::AbstractQuasiCyclicCode)

Return `true` if `C` is a single-generator quasi-cyclic code.
"""
issinglegenerator(C::AbstractQuasiCyclicCode) = (nrows(C.A) == 1;)

function generatormatrix(C::AbstractQuasiCyclicCode, standform::Bool=false)
    if ismissing(C.G)
        if C.Atype == 'G'
            G = lift(C.A);
        else
            _, G = right_kernel(lift(C.A))
            G = transpose(G)
        end
        C.G = G
    end
    standform ? (return _standardform(C.G);) : (return C.G;)
end

function paritycheckmatrix(C::AbstractQuasiCyclicCode, standform::Bool=false)
    if ismissing(C.H)
        if C.Atype == 'H'
            H = lift(C.A)
        else
            _, H = right_kernel(lift(C.A))
            H = transpose(H)
        end
        C.H = H
    end
    standform ? (return _standardform(C.H);) : (return C.H;)
end

"""
    noncirculantgeneratormatrix(C::AbstractQuasiCyclicCode)

Return the non-circulant form of the generator matrix for the quasi-cyclic code `C` if the
polynomial matrix specifies the generator matrix. Otherwise, return `missing`.
"""
function noncirculantgeneratormatrix(C::AbstractQuasiCyclicCode)
    if C.Atype == 'G'
        flag = true
        for r in 1:C.m
            Ginner = zero_matrix(C.F, C.l, C.l * C.n)
            for col in 1:C.l
                C = lift(C.A[r, col])
                Ginner = zero(C)
                for i in 1:m
                    c = k % l
                    c == 0 && (c = l;)
                    Ginner[:, (i - 1) * l + c] = C[:, i]
                end
            end
            flag ? (G = Ginner;) : (G = vcat(G, Ginner);)
        end
        return G
    else
        return missing
    end
end

"""
    noncirculantparitycheckmatrix(C::AbstractQuasiCyclicCode)

Return the non-circulant form of the parity-check matrix for the quasi-cyclic code `C`
if the polynomial matrix specifies the parity-check matrix. Otherwise, return `missing`.
"""
function noncirculantparitycheckmatrix(C::AbstractQuasiCyclicCode)
    if C.Atype == 'H'
        flag = true
        for r in 1:C.m
            Hinner = zero_matrix(C.F, C.l, C.l * C.n)
            for col in 1:C.l
                C = lift(C.A[r, col])
                Hinner = zero(C)
                for i in 1:m
                    c = k % l
                    c == 0 && (c = l;)
                    Hinner[:, (i - 1) * l + c] = C[:, i]
                end
            end
            flag ? (H = Hinner;) : (H = vcat(H, Hinner);)
        end
        return H
    else
        return missing
    end
end

"""
    generators(C::AbstractQuasiCyclicCode)

Return the generators of the quasi-cyclic code.
"""
function generators(C::AbstractQuasiCyclicCode)
    G = noncirculantgeneratormatrix(C)
    genvecs = Vector{fq_nmod_mat}()
    nr = nrows(G)
    for i in 1:nr
        if i % l == 1
            push!(genvecs, G[i, :])
        end
    end
    return genvecs
end

"""
    circulants(C::AbstractQuasiCyclicCode)

Return the circulant matrices of the quasi-cyclic code.
"""
function circulants(C::AbstractQuasiCyclicCode)
    circulants = Vector{fq_nmod_mat}()
    nr, nc = size(C.A)
    # want stored in row order
    for r in 1:nr
        for c in 1:nc
            push!(circulants, lift(C.A[r, c]))
        end
    end
    return circulants
end
