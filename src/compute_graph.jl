import LightGraphs
const LG = LightGraphs

struct ComputeGraph
    graph::LightGraphs.SimpleGraphs.SimpleDiGraph{Int64}
    name2id::Dict{Symbol, Int}
    id2name::Dict{Int, Symbol}
end

function ComputeGraph()
    ComputeGraph(LG.DiGraph(), Dict{Symbol, Int}(), Dict{Int, Symbol}())
end

struct Vertex
    id::Int
    name::Symbol
end

hasvertex(G::ComputeGraph, name::Symbol) = haskey(G.name2id, name)

function addvertex!(G::ComputeGraph, name::Symbol)
    if !hasvertex(G, name)
        LG.add_vertex!(G.graph)
        id = LG.nv(G.graph)
        G.name2id[name] = id
        G.id2name[id] = name

        return Vertex(id, name)
    else
        return Vertex(G.name2id[name], name)
    end
end
getvertex!(G::ComputeGraph, name::Symbol) = addvertex!(G, name)
Base.getindex(G::ComputeGraph, name::Symbol) = getvertex!(G, name)

function addedge!(G::ComputeGraph, from::Symbol, to::Symbol)
    v_from = getvertex!(G, from)
    v_to = getvertex!(G, to)

    LG.add_edge!(G.graph, LG.Edge(v_from.id, v_to.id))
end
addedge!(G::ComputeGraph, from::Vertex, to::Vertex) = LG.add_edge!(G.graph, LG.Edge(from.id, to.id))

in_neighbors(G::ComputeGraph, n::Symbol) = in_neighbors(G, G[n])
function in_neighbors(G::ComputeGraph, v::Vertex)
    nbs = LG.in_neighbors(G.graph, v.id)
    map(nb -> Vertex(nb, G.id2name[nb]), nbs)
end

function out_neighbors(G::ComputeGraph, v::Vertex)
    nbs = LG.out_neighbors(G.graph, v.id)
    map(nb -> Vertex(nb, G.id2name[nb]), nbs)
end
out_neighbors(G::ComputeGraph, n::Symbol) = out_neighbors(G, G[n])
