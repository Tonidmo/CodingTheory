# Copyright (c) 2021, 2023 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    ReedMullergeneratormatrix(r::Int, m::Int, alt::Bool=false)

Return the recursive form of the generator matrix for the ``\\mathcal{RM}(r, m)`` Reed-Muller code.

# Notes
* If `alt` is `true`, the identity is used for the generator matrix for ``\\mathcal{RM}(1, 1)``, as in common in some sources.
  Otherwise, `[1 1; 0 1]` is used, as is common in other sources.
"""
function ReedMullergeneratormatrix(r::Int, m::Int, alt::Bool=false)
    (0 ≤ r ≤ m) || throw(DomainError("Reed-Muller codes require 0 ≤ r ≤ m, received r = $r and m = $m."))

    F = GF(2)
    if r == 1 && m == 1 && !alt
        return matrix(F, 2, 2, [1, 1, 0, 1])
    elseif r == m
        return identity_matrix(F, 2^m)
    elseif r == 0
        return matrix(F, ones(Int, 1, 2^m))
    else
        Grm1 = ReedMullergeneratormatrix(r, m - 1, alt)
        Gr1m1 = ReedMullergeneratormatrix(r - 1, m - 1, alt)
        return vcat(hcat(Grm1, Grm1), hcat(zero_matrix(F, nrows(Gr1m1), ncols(Gr1m1)), Gr1m1))
    end
end

"""
    ReedMullerCode(r::Int, m::Int, alt::Bool=false)

Return the ``\\mathcal{RM}(r, m)`` Reed-Muller code.

# Notes
* If `alt` is `true`, the identity is used for the generator matrix for ``\\mathcal{RM}(1, 1)``, as in common in some sources.
  Otherwise, `[1 1; 0 1]` is used, as is common in other sources.
"""
function ReedMullerCode(r::Int, m::Int, alt::Bool=false)
    (0 ≤ r < m) || throw(DomainError("Reed-Muller codes require 0 ≤ r < m, received r = $r and m = $m."))
    m < 64 || throw(DomainError("This Reed-Muller code requires the implmentation of BigInts. Change if necessary."))

    G = ReedMullergeneratormatrix(r, m, alt)
    H = ReedMullergeneratormatrix(m - r - 1, m, alt)
    Gstand, Hstand, P, rnk = _standardform(G)

    # verify
    ncols(G) == 2^m || error("Generator matrix computed in ReedMuller has the wrong number of columns; received: $(ncols(G)), expected: $(BigInt(2)^m).")
    k = sum([binomial(m, i) for i in 0:r])
    nrows(G) == k || error("Generator matrix computed in ReedMuller has the wrong number of rows; received: $(nrows(G)), expected: $k.")
    size(H) == (2^m - k, k) && (H = transpose(H);)
    d = 2^(m - r)

    if r == 1
        R, vars = PolynomialRing(Nemo.ZZ, 2)
        poly = vars[1]^(2^m) + (2^(m + 1) - 2) * vars[1]^(2^m - 2^(m - 1))*vars[2]^(2^(m - 1)) + vars[2]^(2^m)
        return ReedMullerCode(base_ring(G), ncols(G), nrows(G), d, d, d, r, m, G,
            H, Gstand, Hstand, P, WeightEnumerator(poly, :complete))
    end

    return ReedMullerCode(base_ring(G), ncols(G), nrows(G), d, d, d, r, m, G,
        H, Gstand, Hstand, P, missing)
end

#############################
      # getter functions
#############################

"""
    order(C::ReedMullerCode)
    RMr(C::ReedMullerCode)

Return the order, `r`, of the ``\\mathcal{RM}(r, m)`` Reed-Muller code.
"""
order(C::ReedMullerCode) = C.r
RMr(C::ReedMullerCode) = order(C)

"""
    numberofvariables(C::ReedMullerCode)
    RMm(C::ReedMullerCode)

Return the number of variables, `m`, of the ``\\mathcal{RM}(r, m)`` Reed-Muller code.
"""
numberofvariables(C::ReedMullerCode) = C.m
RMm(C::ReedMullerCode) = C.m

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

"""
    dual(C::ReedMullerCode)

Return the dual of the Reed-Muller code `C`.
"""
function dual(C::ReedMullerCode)
    d = 2^(C.r + 1)
    return ReedMullerCode(C.F, C.n, C.n - C.k, d, d, d, C.m - C.r - 1, C.m, C.H, C.G,
        C.Hstand, C.Gstand, C.Pstand, missing)
end

# "On products and powers of linear codes under componentwise multiplication"
# Hugues Randriambololona
"""
    entrywiseproductcode(C::ReedMullerCode, D::ReedMullerCode)
    *(C::ReedMullerCode, D::ReedMullerCode)
    Schurproductcode(C::ReedMullerCode, D::ReedMullerCode)
    Hadamardproductcode(C::ReedMullerCode, D::ReedMullerCode)
    componentwiseproductcode(C::ReedMullerCode, D::ReedMullerCode)

Return the entrywise product of `C` and `D`.

# Notes
* This is known to often be the full ambient space.
"""
function entrywiseproductcode(C::ReedMullerCode, D::ReedMullerCode)
    C.F == D.F || throw(ArgumentError("Codes must be over the same field in the Schur product."))
    C.n == D.n || throw(ArgumentError("Codes must have the same length in the Schur product."))

    r = C.r + D.r
    if r <= C.n
        return ReedMullerCode(r, C.m)
    else
        return ReedMullerCode(C.m, C.m)
    end
end
*(C::ReedMullerCode, D::ReedMullerCode) = entrywiseproductcode(C, D)
Schurproductcode(C::ReedMullerCode, D::ReedMullerCode) = entrywiseproductcode(C, D)
Hadamardproductcode(C::ReedMullerCode, D::ReedMullerCode) = entrywiseproductcode(C, D)
componentwiseproductcode(C::ReedMullerCode, D::ReedMullerCode) = entrywiseproductcode(C, D)
