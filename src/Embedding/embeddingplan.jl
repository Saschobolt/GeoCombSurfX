# using Primes
# using Combinatorics

include("../Framework.jl")
include("../SimplicialSurface.jl")

function embeddingplan(f::AbstractSimpleGraph, d::Integer = 3)
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
    #   Branch and bound approach:
    #   Objective: Find plan with lowest number of new edges and cost.
    #   First find a solution to the problem with deep seach + greedy heuristic. Then do a deep search for 
    #   the global minimum  by going vertex by vertex and keeping track of the number of necessary edges to embed the vertex.
    #   Abandon the branch of the seach tree, if the req_edges_cost(branch) > req_edges_cost(heuristic) of the heuristic solution.
    ################################################################################################
    # n = length(get_verts(f))
    
    # # id of vertex subset M
    # id_set(M::Vector{<:Integer}) = sum([2^(m-1) for m in M])
    
    # ## lookup table for the cost of embedding vertex v to the subset M (decoded as an Integer id)
    # lookup_cost_set = Dict{Tuple{Int, Int}, Int}()
    
    # no_req_edges(M::Vector{<:Integer}, v::Integer) = length(filter(e -> v in e && length(Base.intersect(M, e)) == 1, get_edges(f)))
    
    # # cost of embedding a vertex v if the vertices of a subset M are already embedded.
    # cost(M::Vector{<:Integer}, v::Integer) = get!(lookup_cost_set, (id_set(M), v), min(length(M), d) - no_req_edges(M, v))
    
    # id_perm(M::Vector{<:Integer}) = prod(prime(i)^(M[i]) for i in eachindex(M))

    # # lookup table for the concatenation of edges and cost of a plan (decoded by an Integer id)
    # lookup_cost_perm = Dict{Int, Vector{Int}}()

    # # cost of a plan: concatenation of the cost of each step in the plan.
    # cost(plan::Vector{<:Integer}) = get!(lookup_cost_perm, (id_perm(plan)), vcat(0, Vector{Int}([cost(plan[1:i-1], plan[i]) for i in 2:length(plan)])))

    # # number of required edges for the plan to yield a rigid framework in every step
    # function no_req_edges(plan::Vector{<:Integer})
    #     a = cost(plan)
    #     return sum(a[a.>0])
    # end

    # # concatenation of no_req_edges(plan) and cost(plan)
    # req_edges_cost(plan::Vector{<:Integer}) = vcat(no_req_edges(plan), cost(plan))


    # # cost of embedding a vertex v if the vertices of a subset M are already embedded. lookahead is the number of steps the cost is evaluated greedily into the future.
    # function cost(M::Vector{<:Integer}, v::Integer, lookahead::Integer)
    #     # if lookahead is too large, reduce it to the number of vertices left to embed
    #     lookahead = min(lookahead, n - length(M) - 1)

    #     # recursion base case
    #     if lookahead == 0
    #         return get!(lookup_cost, (id(M), v, lookahead), [cost(M, v)])
    #     end
        
    #     if length(setdiff(collect(1:n), vcat(M, v))) == 1
    #         w = setdiff(collect(1:n), vcat(M, v))[1]
    #         return vcat([cost(M,v)], cost(vcat(M, v), w, 0))
    #     end

    #     return get!(lookup_cost, (id(M), v, lookahead), 
    #             vcat(cost(M,v), min([get!(lookup_cost, (id(vcat(M, v)), w, lookahead - 1), cost(vcat(M, v), w, lookahead-1)) for w in setdiff(collect(1:n), vcat(M, [v]))]...)))
    # end

    # # find a heuristic solution using greedy algorithm -> the next vertex to embed is the one with minimal cost in this step
    # heuristic = Vector{Int}(undef, n)

    # for i in 1:n-1
    #     heuristic[i] = argmin(v -> cost(heuristic[1:i-1], v), setdiff(collect(1:n), heuristic[1:i-1]))
    # end

    # heuristic[n] = setdiff(collect(1:n), heuristic[1:n-1])[1]


    # lookup_plans = Dict(Int, Vector{Int})()
    # optimal_costs = Dict(Int, Vector{Int})()

    # function last_vertex(M)
    #     argmin(v -> get!(optimal_costs, id(setdiff(M, [v])), get!(lookup_plans, id(setdiff(M, [v])), []))
    # end

    # function optimal_plan(M)
    #     if length(M) == 1
    #         return get!(lookup_plans, id(M), M)
    #     end

    #     best = argmin(N -> get!(optimal_costs, id(N), req_edges_cost(optimal_plan(N)), collect(combinations(M, length(M)-1))))

    #     get!(lookup_plans, id(M), argmin(N -> ))
    # end

    # queue = sort([[i] for i in 1:n], by = M -> vcat(n - length(M), req_edges_cost(get!(lookup_plans, id_set(M), M))))

    # # branch and bound method
    # branch(M::Vector{<:Integer}) = setdiff([sort(vcat(M, v)) for v in setdiff(collect(1:n), M)], queue)

    # bound(M::Vector{<:Integer}) = 

    # best = heuristic

    # queue = sort([[i] for i in 1:n], by = M -> vcat(n - length(M), bound(M)))

    # while length(queue) > 0
    #     # display(queue)
    #     M = popfirst!(queue)
    #     display(M)
    #     display(vcat(n - length(M), bound(M)))
    #     if bound(M) < req_edges_cost(best)
    #         for M_new in branch(M)
    #             if length(M_new) == n
    #                 if req_edges_cost(M_new) < req_edges_cost(best)
    #                     best = M_new
    #                 end
    #             else
    #                 pushfirst!(queue, M_new)
    #             end
    #         end
    #     end
    #     sort(queue, by = M -> vcat(n - length(M), bound(M)))
    # end

    # return best, cost(best), no_req_edges(best)


    # n = length(get_verts(f))

    # id(M::Vector{<:Integer}) = sum([2^(m-1) for m in M])

    # lookup_tuple = Dict(Tuple{Int, Int}, Vector{Int})()
    # lookup_plans = Dict(Int, Vector{Int})()

    # cost(M::Vector{<:Integer}, v::Integer) = get!(lookup_tuple, (id(M), v), min(length(M), d) - length(filter(e -> v in e && length(Base.intersect(M, e)) == 1, get_edges(f))))
    # cost(M::Vector{<:Integer}) = get!(lookup_plans, id(M), vcat(0, [cost(M[1:i-1], M[i]) for i in 2:length(M)]))

    # function plan(M::Vector{<:Integer})
    #     return          get!(lookup_plans, id(M), argmin(v -> cost(M, v), setdiff(collect(1:n), M)))
    # end




    no_req_edges(M::Vector{<:Integer}, v::Integer) = length(filter(e -> v in e && length(Base.intersect(M, e)) == 1, get_edges(f)))
    
    # cost of embedding a vertex v if the vertices of a subset M are already embedded.
    cost(M::Vector{<:Integer}, v::Integer) = min(length(M), d) - no_req_edges(M, v)

    # cost of a plan: concatenation of the cost of each step in the plan.
    cost(plan::Vector{<:Integer}) = vcat(0, Vector{Int}([cost(plan[1:i-1], plan[i]) for i in 2:length(plan)]))


    # number of required edges for the plan to yield a rigid framework in every step
    function no_req_edges(plan::Vector{<:Integer})
        a = cost(plan)
        return sum(a[a.>0])
    end

    # concatenation of no_req_edges(plan) and cost(plan)
    req_edges_cost(plan::Vector{<:Integer}) = vcat(no_req_edges(plan), cost(plan))

    S_comp = get_verts(f)
    E_comp = get_edges(f)
    S = [argmax(v -> count(x -> x == v, vcat(E_comp...)), S_comp)]
    E = filter(e -> S[1] in e, E_comp)
    E_comp = setdiff(E_comp, E)
    S_comp = setdiff(S_comp, S)

    while length(S_comp) > 0
        next = argmax(v -> (count(x -> x == v, vcat(E...)), count(x -> x == v, vcat(E_comp...))), S_comp)
        push!(S, next)
        push!(E, filter(e -> next in e, E_comp)...)
        E_comp = setdiff(E_comp, E)
        S_comp = setdiff(S_comp, S)
    end

    return S, cost(S), no_req_edges(S)
end

function embeddingplan(f::AbstractEmbeddedGraph)
    g = Graph(verts = collect(1:size(get_verts(f))[2]), edges = get_edges(f))

    return embeddingplan(g, dimension(f))
end