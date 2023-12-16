using Combinatorics

include("../Framework.jl")
include("../SimplicialSurface.jl")

function embeddingplan(f::AbstractSimpleGraph, d::Integer = 3; lookahead::Integer = 3)
    n = length(get_verts(f))

    # cost of embedding a vertex v if the vertices of a subset M are already embedded.
    function cost(M::Vector{<:Integer}, v::Integer)
        return min(length(M), d) - length(filter(e -> v in e && length(Base.intersect(M, e)) == 1, get_edges(f)))
    end

    # total cost of embedding the vertices in plan in the order they appear
    function cost(plan::Vector{<:Integer})
        return sum([max(cost(plan[1:i-1], plan[i]), 0) for i in 2:length(plan)])
    end
  
    ################################################################################################
    #       Dynamic Program to find the cost of a subset of vertices -> very expensive as complexity is O(n * 2^n)
    ################################################################################################

    # # convert a subset of the vertices of f into an integer id of that subset by translaing the characterictic function of the subset into a sum of powers of 2. [1,2,3] -> 2^0 + 2^1 + 2^2 = 7.
    # function id(M::Vector{<:Integer})
    #     return sum(2^(v-1) for v in M)
    # end

    # # lookup table for the cost of a subset.
    # lookup = Vector{Int}(undef, 2^n - 1)

    # # initialize lookop table with costs of 1-element subsets, which are zero.
    # for v in 1:n
    #     lookup[2^(v-1)] = 0
    # end

    # # fill the look up table with the costs of all subsets of size 2, 3, ..., n.
    # for l in 2:n
    #     for M in combinations(collect(1:n), l)
    #         lookup[id(M)] = min([lookup[id(setdiff(M, [v]))] + cost(setdiff(M, [v]), v) for v in M]...)
    #     end
    # end

    # return lookup[2^n-1]

    # ################################################################################################
    # #    Simulated Annealing approach to find a good embedding -> polynomial complexity O(n^2)
    # ################################################################################################

    # current = collect(1:n)  # start with a random solution
    # best = current
    # T = 20  # initial temperature
    # T_min = 0.001
    # alpha = 0.9

    # while T > T_min
    #     for i in 1:100
    #         new_solution = current
    #         # make a small change to the solution
    #         idx = rand(1:n-1)
    #         new_solution[[idx, idx+1]] = new_solution[[idx+1, idx]]

    #         cost_current = cost(current)
    #         cost_new = cost(new_solution)

    #         if cost_new < cost_current || exp(-(cost_new - cost_current) / T) > rand()
    #             current = new_solution
    #             if cost_new < cost(best)
    #                 best = new_solution
    #             end
    #         end
    #     end
    #     T = T * alpha
    # end

    # return best, cost(best)

    ################################################################################################
    #    Greedy approach to find a good embedding
    ################################################################################################
    # cost of embedding a vertex v if the vertices of a subset M are already embedded. lookahead is the number of steps the cost is evaluated greedily into the future.
    function cost(M::Vector{<:Integer}, v::Integer, lookahead::Integer)
        # if lookahead is too large, reduce it to the number of vertices left to embed
        lookahead = min(lookahead, n - length(M) - 1)

        # recursion base case
        if lookahead == 0
            return [cost(M, v)]
        end
        
        if length(setdiff(collect(1:n), vcat(M, [v]))) == 1
            w = setdiff(collect(1:n), vcat(M, [v]))[1]
            return vcat([cost(M,v)], cost(vcat(M, [v]), w, 0))
        end

        return vcat([cost(M,v)], min([cost(vcat(M, [v]), w, lookahead-1) for w in setdiff(collect(1:n), vcat(M, [v]))]...))
    end

    plan = Vector{Int}(undef, n)

    for i in 1:n-1
        plan[i] = argmin(v -> cost(plan[1:i-1], v, lookahead), setdiff(collect(1:n), plan[1:i-1]))
    end

    plan[n] = setdiff(collect(1:n), plan[1:n-1])[1]

    return plan, cost(plan)
end

function embeddingplan(f::AbstractEmbeddedGraph, d::Integer = 3; lookahead::Integer = 3)
    g = Graph(verts = collect(1:size(get_verts(f))[2]), edges = get_edges(f))

    return embeddingplan(g, d; lookahead = lookahead)
end