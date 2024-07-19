using Combinatorics
using Groebner
using Nemo

# number of vertices
n = 9
# 
d = 2

# big polynomial ring
variables = Nemo.AbstractAlgebra.variable_names(:x => combinations(1:n, d + 1))
R, x = polynomial_ring(QQ, variables)
variable_dict = Dict{Vector{Int},typeof(x[1])}()

for (i, bracket) in enumerate(combinations(1:n, d + 1))
    variable_dict[bracket] = x[i]
end

# monomial ordering
ordering = Groebner.DegRevLex(x)

function sizyges(d, n)
    # Sturmfels: Algorithms in invariant theory p.81 & p.84 exercise 3
    # in Sturmfels d is length of brackets, for us it is dimension. So every d in Sturmfels needs to be substituted by d+1

    variables = Nemo.AbstractAlgebra.variable_names(:x => combinations(1:n, d + 1))
    R, x = polynomial_ring(QQ, variables)
    variable_dict = Dict{Vector{Int},typeof(x[1])}()

    for (i, bracket) in enumerate(combinations(1:n, d + 1))
        variable_dict[bracket] = x[i]
    end

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

        return (sum(τ -> summand(α, β, γ, τ), combinations(1:d+2, s)) for α in combinations(1:n, s - 1), β in combinations(1:n, d + 2), γ in combinations(1:n, d + 1 - s) if α <= β)
    end

    return (sum(τ -> summand(α, β, γ, τ), combinations(1:d+2, s)) for s in 1:d for α in combinations(1:n, s - 1), β in combinations(1:n, d + 2), γ in combinations(1:n, d + 1 - s) if α <= β)
end

# collect(sizyges(d, n))