mutable struct BracketAlgebra
    d::Int
    n::Int
    R::QQMPolyRing
    x::Vector{Nemo.QQMPolyRingElem}
    variable_dict::Dict{Vector{Int},Nemo.QQMPolyRingElem}
    ordering::DegRevLex{QQMPolyRingElem}
    groebner_basis::Union{Nothing,Vector{Nemo.QQMPolyRingElem}}

    function BracketAlgebra(d, n)
        variables = Nemo.AbstractAlgebra.variable_names(:x => combinations(1:n, d + 1))
        R, x = polynomial_ring(QQ, variables)
        variable_dict = Dict{Vector{Int},typeof(x[1])}()

        for (i, bracket) in enumerate(combinations(1:n, d + 1))
            variable_dict[bracket] = x[i]
        end

        ordering = Groebner.DegRevLex(x)

        return new(d, n, R, x, variable_dict, ordering, nothing)
    end
end

function sizyges(B::BracketAlgebra)
    # Sturmfels: Algorithms in invariant theory p.81 & p.84 exercise 3
    # in Sturmfels d is length of brackets, for us it is dimension. So every d in Sturmfels needs to be substituted by d+1

    n = B.n
    d = B.d
    R = B.R
    x = B.x
    variable_dict = B.variable_dict

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
    return prod(B.variable_dict[row] for row in t.rows)
end

function nonstandard_tabloids(B::BracketAlgebra)
    d = B.d
    n = B.n
    return (Tabloid([row1, row2]) for row1 in combinations(1:n, d + 1), row2 in combinations(1:n, d + 1) if (row1 < row2 && !is_standard(Tabloid([row1, row2]))))
end

function reduced_groebner_basis(B::BracketAlgebra)
    if !isnothing(B.groebner_basis)
        return B.groebner_basis
    end

    basis = sizyges_vector(B)
    tobereduced = [bracket_monomial(t, B) for t in collect(nonstandard_tabloids(B))]

    reduced = tobereduced .- Groebner.normalform(basis, tobereduced, ordering=B.ordering)
    filter!(b -> b != 0, reduced)

    B.groebner_basis = reduced
    return reduced
end