abstract type AbstractBracketAlgebra <: Nemo.Ring end

mutable struct BracketAlgebra <: AbstractBracketAlgebra
    d::Int
    n::Int
    R::QQMPolyRing
    variables::Dict{Vector{Int},Nemo.QQMPolyRingElem}
    ordering::DegRevLex{QQMPolyRingElem}
    groebner_basis::Union{Nothing,Vector{Nemo.QQMPolyRingElem}}

    function BracketAlgebra(n, d)
        vars = Nemo.AbstractAlgebra.variable_names(:x => combinations(1:n, d + 1))
        R, x = polynomial_ring(QQ, vars; internal_ordering=:degrevlex)
        variable_dict = Dict{Vector{Int},typeof(x[1])}()

        for (i, bracket) in enumerate(combinations(1:n, d + 1))
            variable_dict[bracket] = x[i]
        end

        ordering = Groebner.DegRevLex(x)

        return new(d, n, R, variable_dict, ordering, nothing)
    end
end

function BracketAlgebra(g::Graphs.AbstractSimpleGraph, d::Integer=2)
    return BracketAlgebra(Graphs.nv(g), d)
end

function BracketAlgebra(poly::AbstractEmbOrCombPolyhedron, d::Integer=3)
    return BracketAlgebra(Graphs.SimpleGraph(poly), d)
end

function sizyges(B::BracketAlgebra)
    # Sturmfels: Algorithms in invariant theory p.81 & p.84 exercise 3
    # in Sturmfels d is length of brackets, for us it is dimension. So every d in Sturmfels needs to be substituted by d+1

    n = B.n
    d = B.d
    R = B.R
    variable_dict = B.variables

    function sign(λ::AbstractVector{<:Integer}, k)
        # Sturmfels p.79
        λ_ast = setdiff(collect(1:k), λ)
        perm = Perm(vcat(λ, λ_ast))^(-1)
        return Nemo.AbstractAlgebra.sign(perm)
    end

    function summand(α, β, γ, τ)
        τ_ast = setdiff(collect(1:d+2), τ)
        if length(unique(vcat(α, β[τ_ast]))) < length(vcat(α, β[τ_ast]))
            return 0
        elseif length(unique(vcat(β[τ], γ))) < length(vcat(β[τ], γ))
            return 0
        end

        return sign(τ, d + 2) * variable_dict[sort(vcat(α, β[τ_ast]))] * variable_dict[sort(vcat(β[τ], γ))]
    end

    if d <= 2
        s = 1

        return (sum(τ -> summand(α, β, γ, τ), combinations(1:d+2, s)) for α in combinations(1:n, s - 1), β in combinations(1:n, d + 2), γ in combinations(1:n, d + 1 - s) if (α <= β && (length(α) > 0 ? α[end] < β[s+1] : true) && (β[s] < γ[1])))
    end

    return (sum(τ -> summand(α, β, γ, τ), combinations(1:d+2, s)) for s in 1:d for α in combinations(1:n, s - 1), β in combinations(1:n, d + 2), γ in combinations(1:n, d + 1 - s) if (α <= β && (length(α) > 0 ? α[end] < β[s+1] : true) && (β[s] < γ[1])))
end

function sizyges_vector(B)
    vec = collect(sizyges(B))
    filter!(x -> x != 0, vec)
    vec = leading_coefficient.(vec) .^ (-1) .* vec
    unique!(vec)
    return vec
end

mutable struct Tabloid{T<:Integer}
    rows::Vector{Vector{T}}

    function Tabloid(rows::AbstractVector{<:AbstractVector{T}}) where {T<:Integer}
        @assert all(length(rows[1]) == length(row) for row in rows)
        return new{T}(sort(sort.(rows)))
    end
end

function Matrix(t::Tabloid)
    return transpose(hcat(t.rows...))
end

# function Base.show(io::IO, t::Tabloid)
#     println(io, "Tabloid with $(length(t.rows)) rows and matrix")
#     show(io, Matrix(t))
# end

# function Base.show(io::IO, ::MIME"text/plain", t::Tabloid)
#     println(io, "Tabloid with $(length(t.rows)) rows and matrix")
#     display(Matrix(t))
# end

function Base.vcat(t::Tabloid...)
    return Tabloid(vcat([t[i].rows for i in 1:length(t)]...))
end

function is_standard(t::Tabloid)
    # tabloid is standard, if the entries in each column are increasing
    mat = Matrix(t)

    for i in 1:size(mat, 2)
        col = mat[:, i]
        if any(col[1:end-1] .> col[2:end])
            return false
        end
    end

    return true
end

function bracket_monomial(t::Tabloid, B::BracketAlgebra)
    return prod(B.variables[row] for row in t.rows)
end

function nonstandard_tabloids(B::BracketAlgebra)
    d = B.d
    n = B.n
    return (Tabloid([row1, row2]) for row1 in combinations(1:n, d + 1), row2 in combinations(1:n, d + 1) if (row1 < row2 && !is_standard(Tabloid([row1, row2]))))
end

function reduced_groebner_basis!(B::BracketAlgebra)
    if !isnothing(B.groebner_basis)
        return B.groebner_basis
    end

    basis = sizyges_vector(B)
    tobereduced = [bracket_monomial(t, B) for t in collect(nonstandard_tabloids(B))]

    if length(tobereduced) == 0
        return basis
    end

    reduced = tobereduced .- Groebner.normalform(basis, tobereduced, ordering=B.ordering)
    filter!(b -> b != 0, reduced)

    B.groebner_basis = reduced
    return reduced
end

abstract type AbstractBracketAlgebraElem end

mutable struct BracketAlgebraElem <: AbstractBracketAlgebraElem
    parent::BracketAlgebra
    polynomial::Nemo.QQMPolyRingElem
end

# function Base.display(b::BracketAlgebraElem)
#     exponents = Nemo.exponent_vectors(b)
#     str = ""

#     for (i, exp) in enumerate(exponents)
#         coeff = Nemo.coeff(b, i)

#         if Nemo.sign(coeff) == -1
#             str = str * " - "
#         elseif i > 1
#             str = str * " + "
#         end

#         if !(coeff in [Nemo.one(Nemo.base_ring(b)), -Nemo.one(Nemo.base_ring(b))])
#             str = str * "$coeff"
#         end

#         for (j, val) in enumerate(exp)
#             if val == 0
#                 continue
#             elseif val == 1
#                 str = str * "$(collect(keys(parent(b).variables))[j])"
#             else
#                 str = str * "$(collect(keys(parent(b).variables))[j])" * "^$val"
#             end
#         end
#     end

#     display(str)
# end

Base.display(b::BracketAlgebraElem) = display(b.polynomial)

Base.parent(b::BracketAlgebraElem) = b.parent
Nemo.elem_type(::BracketAlgebra) = BracketAlgebraElem
Nemo.parent_type(::BracketAlgebraElem) = BracketAlgebra
Nemo.base_ring(b::BracketAlgebraElem) = Nemo.base_ring(Base.parent(b).R)
Nemo.base_ring(B::BracketAlgebra) = Nemo.base_ring(B.R)

Base.one(B::BracketAlgebra) = BracketAlgebraElem(B, Base.one(B.R))
Base.zero(B::BracketAlgebra) = BracketAlgebraElem(B, Base.zero(B.R))
(B::BracketAlgebra)(A::Vector{T}, m::Vector{Vector{Int}}) where {T<:Nemo.RingElem} = BracketAlgebraElem(B, B.R(A, m))
(B::BracketAlgebra)(p::Nemo.MPolyRingElem) = BracketAlgebraElem(B, p)

# Bracket expression from array: B([1,2,3,4]) = [1,2,3,4]
(B::BracketAlgebra)(bracket::Vector{<:Integer}) = length(unique(bracket)) == length(bracket) ? BracketAlgebraElem(B, B.variables[sort(bracket)]) : zero(B)

# Bracket polynomial from array of array of arrays. They encode the bracket polynomial as a sum of monomials. B([[[1,2], [3,4]], [2,3]]) = [1,2]*[3,4] + [2,3]
(B::BracketAlgebra)(A::Vector{<:Vector{<:Vector{<:Integer}}}) = sum(prod(B(bracket) for bracket in monomial) for monomial in A)

Nemo.length(b::BracketAlgebraElem) = Nemo.length(b.polynomial)
Nemo.degrees(b::BracketAlgebraElem) = Nemo.degrees(b.polynomial)
Nemo.total_degree(b::BracketAlgebraElem) = Nemo.total_degree(b.polynomial)
Nemo.coefficients(b::BracketAlgebraElem) = Nemo.coefficients(b.polynomial)
Nemo.monomials(b::BracketAlgebraElem) = (parent(b)(p) for p in Nemo.monomials(b.polynomial))
Nemo.terms(b::BracketAlgebraElem) = (parent(b)(p) for p in Nemo.terms(b.polynomial))
Nemo.exponent_vectors(b::BracketAlgebraElem) = Nemo.exponent_vectors(b.polynomial)
Nemo.coeff(b::BracketAlgebraElem, n::Int) = Nemo.coeff(b.polynomial, n)
Nemo.coeff(b::BracketAlgebraElem, exps::Vector{Int}) = Nemo.coeff(b.polynomial, exps)
Nemo.monomial(b::BracketAlgebraElem, n::Int) = parent(b)(Nemo.monomial(b.polynomial, n))
Nemo.term(b::BracketAlgebraElem, n::Int) = parent(b)(Nemo.term(b.polynomial, n))

# return all brackets that appear in b as arrays
brackets(b::BracketAlgebraElem) = collect(keys(parent(b).variables))[sum(Nemo.exponent_vectors(b)).>0]


Base.:*(a::BracketAlgebraElem, b::BracketAlgebraElem) = BracketAlgebraElem(a.parent, a.polynomial * b.polynomial)
Base.:+(a::BracketAlgebraElem, b::BracketAlgebraElem) = BracketAlgebraElem(a.parent, a.polynomial + b.polynomial)
Base.:-(a::BracketAlgebraElem, b::BracketAlgebraElem) = BracketAlgebraElem(a.parent, a.polynomial - b.polynomial)
Base.:-(b::BracketAlgebraElem) = BracketAlgebraElem(b.parent, -b.polynomial)
Base.:^(b::BracketAlgebraElem, n::Int) = BracketAlgebraElem(b.parent, b.polynomial^n)


# function evaluate(b::BracketAlgebraElem, coordinization::Vector{<:Nemo.RingElem})

# end