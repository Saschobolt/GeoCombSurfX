abstract type AbstractBracketAlgebra <: Nemo.Ring end

mutable struct BracketAlgebra{T<:Union{Nemo.RingElem,Number}} <: AbstractBracketAlgebra
    d::Int
    n::Int
    R::Nemo.MPolyRing{T}
    variables::Dict{Vector{Int},<:Nemo.MPolyRingElem{T}}
    # ordering::DegRevLex{<:Nemo.MPolyRingElem{T}}
    # groebner_basis::Union{Nothing,Vector{<:Nemo.MPolyRingElem{T}}}

    function BracketAlgebra(n, d, T::Type=Nemo.ZZRingElem)
        S = Nemo.parent_type(T)
        brackets = reverse(collect(combinations(1:n, d + 1)))
        vars = Nemo.AbstractAlgebra.variable_names("x#" => brackets)
        R, x = Nemo.polynomial_ring(S(), vars; internal_ordering=:lex)
        variable_dict = Dict{Vector{Int},typeof(x[1])}()

        for (i, bracket) in enumerate(brackets)
            variable_dict[bracket] = x[i]
        end

        # ordering = Groebner.DegRevLex(x)

        return new{Nemo.elem_type(S)}(d, n, R, variable_dict)
    end
end

function BracketAlgebra(g::Graphs.AbstractSimpleGraph, d::Integer=2, T::Type=Nemo.ZZRingElem)
    return BracketAlgebra(Graphs.nv(g), d, T)
end

function BracketAlgebra(poly::AbstractEmbOrCombPolyhedron, d::Integer=3, T::Type=Nemo.ZZRingElem)
    return BracketAlgebra(Graphs.SimpleGraph(poly), d, T)
end

# function sizyges(B::BracketAlgebra)
#     # Sturmfels: Algorithms in invariant theory p.81 & p.84 exercise 3
#     # in Sturmfels d is length of brackets, for us it is dimension. So every d in Sturmfels needs to be substituted by d+1

#     n = B.n
#     d = B.d
#     R = B.R
#     variable_dict = B.variables

#     function sign(λ::AbstractVector{<:Integer}, k)
#         # Sturmfels p.79
#         λ_ast = setdiff(collect(1:k), λ)
#         perm = Perm(vcat(λ, λ_ast))^(-1)
#         return Nemo.AbstractAlgebra.sign(perm)
#     end

#     function summand(α, β, γ, τ)
#         τ_ast = setdiff(collect(1:d+2), τ)
#         if length(unique(vcat(α, β[τ_ast]))) < length(vcat(α, β[τ_ast]))
#             return 0
#         elseif length(unique(vcat(β[τ], γ))) < length(vcat(β[τ], γ))
#             return 0
#         end

#         return sign(τ, d + 2) * variable_dict[sort(vcat(α, β[τ_ast]))] * variable_dict[sort(vcat(β[τ], γ))]
#     end

#     if d <= 2
#         s = 1

#         return (sum(τ -> summand(α, β, γ, τ), combinations(1:d+2, s)) for α in combinations(1:n, s - 1), β in combinations(1:n, d + 2), γ in combinations(1:n, d + 1 - s) if (α <= β && (length(α) > 0 ? α[end] < β[s+1] : true) && (β[s] < γ[1])))
#     end

#     return (sum(τ -> summand(α, β, γ, τ), combinations(1:d+2, s)) for s in 1:d for α in combinations(1:n, s - 1), β in combinations(1:n, d + 2), γ in combinations(1:n, d + 1 - s) if (α <= β && (length(α) > 0 ? α[end] < β[s+1] : true) && (β[s] < γ[1])))
# end

# function sizyges_vector(B)
#     vec = collect(sizyges(B))
#     filter!(x -> x != 0, vec)
#     vec = leading_coefficient.(vec) .^ (-1) .* vec
#     unique!(vec)
#     return vec
# end

# mutable struct Tabloid{T<:Integer}
#     rows::Vector{Vector{T}}

#     function Tabloid(rows::AbstractVector{<:AbstractVector{T}}) where {T<:Integer}
#         @assert all(length(rows[1]) == length(row) for row in rows)
#         return new{T}(sort(sort.(rows)))
#     end
# end

# function Matrix(t::Tabloid)
#     return transpose(hcat(t.rows...))
# end

# # function Base.show(io::IO, t::Tabloid)
# #     println(io, "Tabloid with $(length(t.rows)) rows and matrix")
# #     show(io, Matrix(t))
# # end

# # function Base.show(io::IO, ::MIME"text/plain", t::Tabloid)
# #     println(io, "Tabloid with $(length(t.rows)) rows and matrix")
# #     display(Matrix(t))
# # end

# function Base.vcat(t::Tabloid...)
#     return Tabloid(vcat([t[i].rows for i in 1:length(t)]...))
# end

# function is_standard(t::Tabloid)
#     # tabloid is standard, if the entries in each column are increasing
#     mat = Matrix(t)

#     for i in 1:size(mat, 2)
#         col = mat[:, i]
#         if any(col[1:end-1] .> col[2:end])
#             return false
#         end
#     end

#     return true
# end

# function bracket_monomial(t::Tabloid, B::BracketAlgebra)
#     return prod(B.variables[row] for row in t.rows)
# end

# function nonstandard_tabloids(B::BracketAlgebra)
#     d = B.d
#     n = B.n
#     return (Tabloid([row1, row2]) for row1 in combinations(1:n, d + 1), row2 in combinations(1:n, d + 1) if (row1 < row2 && !is_standard(Tabloid([row1, row2]))))
# end

# function reduced_groebner_basis!(B::BracketAlgebra)
#     if !isnothing(B.groebner_basis)
#         return B.groebner_basis
#     end

#     basis = sizyges_vector(B)
#     tobereduced = [bracket_monomial(t, B) for t in collect(nonstandard_tabloids(B))]

#     if length(tobereduced) == 0
#         return basis
#     end

#     reduced = tobereduced .- Groebner.normalform(basis, tobereduced, ordering=B.ordering)
#     filter!(b -> b != 0, reduced)

#     B.groebner_basis = reduced
#     return reduced
# end

abstract type AbstractBracketAlgebraElem end

mutable struct BracketAlgebraElem{T<:Union{Nemo.RingElem,Number}} <: AbstractBracketAlgebraElem
    parent::BracketAlgebra{T}
    polynomial::Nemo.MPolyRingElem{T}
end

function display(b::BracketAlgebraElem)
    exponents = Nemo.exponent_vectors(b)
    str = ""

    for (i, exp) in enumerate(exponents)
        coeff = Nemo.coeff(b, i)

        if Nemo.sign(coeff) == -1
            str = str * " - "
        elseif i > 1
            str = str * " + "
        end

        if !(coeff in [Nemo.one(Nemo.base_ring(b)), -Nemo.one(Nemo.base_ring(b))])
            str = str * "$coeff"
        end

        for (j, val) in enumerate(exp)
            if val == 0
                continue
            elseif val == 1
                str = str * "$(sort(collect(keys(parent(b).variables)), rev = true)[j])"
            else
                str = str * "$(sort(collect(keys(parent(b).variables)), rev = true)[j])" * "^$val"
            end
        end
    end

    println(str)
end

# display(b::BracketAlgebraElem) = display(b.polynomial)

Base.parent(b::BracketAlgebraElem) = b.parent
Nemo.elem_type(::BracketAlgebra) = BracketAlgebraElem
Nemo.parent_type(::BracketAlgebraElem) = BracketAlgebra
Nemo.base_ring(b::BracketAlgebraElem) = Nemo.base_ring(Base.parent(b).R)
Nemo.base_ring(B::BracketAlgebra) = Nemo.base_ring(B.R)

Nemo.one(B::BracketAlgebra) = BracketAlgebraElem(B, one(B.R))
Nemo.zero(B::BracketAlgebra) = BracketAlgebraElem(B, zero(B.R))
(B::BracketAlgebra)(A::Vector{T}, m::Vector{Vector{Int}}) where {T<:Nemo.RingElem} = BracketAlgebraElem(B, B.R(A, m))
(B::BracketAlgebra)(p::Nemo.MPolyRingElem) = BracketAlgebraElem(B, p)

# Bracket expression from array: B([1,2,3,4]) = [1,2,3,4]
(B::BracketAlgebra)(bracket::Vector{<:Integer}) = length(unique(bracket)) == length(bracket) ? Nemo.sign(Nemo.Perm(Int.(indexin(bracket, sort(bracket)))))BracketAlgebraElem(B, B.variables[sort(bracket)]) : zero(B)

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
Nemo.exponent_vector(b::BracketAlgebra, n::Int) = Nemo.exponent_vector(b.polynomial, n)
Nemo.term(b::BracketAlgebraElem, n::Int) = parent(b)(Nemo.term(b.polynomial, n))
Nemo.leading_term(b::BracketAlgebraElem) = parent(b)(Nemo.leading_term(b.polynomial))
Nemo.leading_monomial(b::BracketAlgebraElem) = parent(b)(Nemo.leading_monomial(b.polynomial))
Nemo.leading_coefficient(b::BracketAlgebraElem) = Nemo.leading_coefficient(b.polynomial)

Nemo.factor(b::BracketAlgebraElem) = Nemo.factor(b.polynomial)

# return all brackets that appear in b as arrays
brackets(b::BracketAlgebraElem) = sort(collect(keys(parent(b).variables)), rev=true)[sum(Nemo.exponent_vectors(b)).>0]

Base.:*(n::Integer, b::BracketAlgebraElem) = BracketAlgebraElem(b.parent, n * b.polynomial)
Base.:*(a::BracketAlgebraElem, b::BracketAlgebraElem) = BracketAlgebraElem(a.parent, a.polynomial * b.polynomial)
Base.:+(a::BracketAlgebraElem, b::BracketAlgebraElem) = BracketAlgebraElem(a.parent, a.polynomial + b.polynomial)
Base.:-(a::BracketAlgebraElem, b::BracketAlgebraElem) = BracketAlgebraElem(a.parent, a.polynomial - b.polynomial)
Base.:-(b::BracketAlgebraElem) = BracketAlgebraElem(b.parent, -b.polynomial)
Base.:^(b::BracketAlgebraElem, n::Int) = BracketAlgebraElem(b.parent, b.polynomial^n)
Base.:>(a::BracketAlgebraElem, b::BracketAlgebraElem) = a.polynomial > b.polynomial
Base.:<(a::BracketAlgebraElem, b::BracketAlgebraElem) = a.polynomial < b.polynomial

function Nemo.evaluate(b::BracketAlgebraElem{T}, A::Vector{T}) where {T<:Union{Nemo.RingElem,Number}}
    Nemo.evaluate(b.polynomial, A)
end

function Nemo.evaluate(b::BracketAlgebraElem{T}, A::Vector{U}) where {T<:Union{Nemo.RingElem,Number},U<:Integer}
    Nemo.evaluate(b.polynomial, A)
end

function Nemo.evaluate(b::BracketAlgebraElem{T}, coordinization::AbstractMatrix{<:Union{Nemo.RingElem,Number}}) where {T<:Union{Nemo.RingElem,Number}}
    bracks = brackets(b)
    A = map(x -> x in bracks ? Nemo.det(T.(hcat(transpose(coordinization[:, x]), ones(parent(b).d + 1, 1)))) : 0, sort(collect(keys(parent(b).variables)), rev=true))
    Nemo.evaluate(b.polynomial, A)
end

mutable struct Tabloid
    matrix::Matrix{Int}
    ordering::Vector{Int}

    """
    Tabloid(matrix::AbstractMatrix{<:Integer}, ordering::AbstractVector{<:Integer}=collect(1:maximum(matrix)))

    Calculate the tabloid whose rows correspond to the rows of matrix with the ordering.
        Example:
        matrix = [4 2 1; 3 1 4], ordering = [3,1,2,4]
        => result: [3 1 4; 1 2 4] 
    """
    function Tabloid(matrix::AbstractMatrix{<:Integer}, ordering::AbstractVector{<:Integer}=collect(1:maximum(matrix)))
        rows = [matrix[i, :] for i in 1:size(matrix)[1]]
        sort!.(rows, by=(x -> indexin(x, ordering)[1]))
        sort!(rows, by=(row -> indexin(row, ordering)))
        return new(transpose(hcat(rows...)), ordering)
    end
end

Base.display(t::Tabloid) = display(t.matrix)

function Tabloid(rows::AbstractVector{<:AbstractVector{<:Integer}}, ordering::AbstractVector{<:Integer}=collect(1:maximum(vcat(rows...))))
    return Tabloid(transpose(hcat(rows...)), ordering)
end

Matrix(t::Tabloid) = t.matrix

Base.vcat(t1::Tabloid, t2::Tabloid) = t1.ordering == t2.ordering ? Tabloid(t1.ordering, vcat(t1.matrix, t2.matrix)) : error("Tabloids need to have same number of columns.")


"""
    tabloid(b::BracketAlgebraElem, ordering::Vector{Int}=collect(1:parent(b).n))

Return the tabloid corresponding to the bracket monomial b ordered by ordering. 
ordering determines the order in which 

Example:
b = [1,2,3][3,4,5], ordering = [4,3,1,2,5]
result: [4 3 5; 3 1 2]
"""
function Tabloid(b::BracketAlgebraElem, ordering::Vector{Int}=collect(1:parent(b).n))
    # see Sturmfels 2008, page 81 on how to build the tabloids
    if length(b) > 1
        error("Only tabloids of bracket monomials can be calculated.")
    end

    exps = collect(Nemo.exponent_vectors(b))[1]
    # all brackets that appear as rows 
    rows = [(repeat(sort(collect(keys(parent(b).variables)), rev=true)[i], exps[i]) for i in eachindex(exps))...]
    filter!(row -> length(row) > 0, rows)
    return Tabloid(rows, ordering)
end

"""
    standard_violation(t::Tabloid)

Return the index of the first violation to the standardness of t, i.e. the first index where t[i,j] > t[i+1, j] with regard to the ordering. Otherwise return nothing.
Example:
t = [1 2 3; 1 4 5; 1 5 6; 2 3 4] with ordering [1,2,3,4,5,6] => return (3,2)
t = [1 2 3; 1 2 4] with ordering [1,2,3,4] => return nothing
"""
function standard_violation(t::Tabloid)
    return findfirst([indexin(t.matrix[row, col], t.ordering)[1] > indexin(t.matrix[row+1, col], t.ordering)[1] for row in 1:size(t.matrix)[1]-1, col in 1:size(t.matrix)[2]])
end

is_standard(t::Tabloid) = isnothing(standard_violation(t))

function straightening_sizyge(α::Vector{<:Integer}, β::Vector{<:Integer}, γ::Vector{<:Integer}, B::BracketAlgebra)
    s = length(α) + 1
    d = B.d
    @assert length(β) == d + 2 "β needs to have length d + 2, but got $(length(β))."
    @assert length(γ) == d + 1 - s "γ needs to have length d + 1 - s, but got $(length(γ))."

    return sum(Nemo.sign(Nemo.Perm(vcat(setdiff(collect(1:d+2), τ), τ))) * B(vcat(α, β[setdiff(collect(1:d+2), τ)])) * B(vcat(β[τ], γ)) for τ in combinations(1:d+2, s))
end