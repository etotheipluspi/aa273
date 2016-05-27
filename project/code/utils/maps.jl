include("intersect.jl")

abstract Building

type Map
    b::Vector{Building}
end
n_buildings(map::Map) = length(map.b)

function rectangle_map(;xs::Vector{Float64}=[3.], 
                        ys::Vector{Float64}=[1.],
                        offs::Matrix{Float64}=zeros(2,1))
    polys = Any[]
    for i = 1:length(xs)
        p = rectangle_poly(xs[i], ys[i], offs[:,i]) 
        push!(polys, p)
    end
    return Map(polys)
end
function rectangle_map(xs::Float64=3., 
                        ys::Float64=1.,
                        offs::Vector{Float64}=[2.,1])
    p = rectangle_poly(xs, ys, offs) 
    return Map([p])
end


type Polygon <: Building
    np::Int64            # n points
    pts::Matrix{Float64} # 2xN
    center::Vector{Float64}
end

function rectangle_poly(xs::Float64, ys::Float64, offset::Vector{Float64}=[0.,0])
    pts = [0. xs xs 0; 0 0 ys ys] .+ offset
    center = [mean(pts[1,:]), mean(pts[2,:])]
    return Polygon(4, pts, center)
end

function n_segments(p::Polygon)
    return p.np
end

function get_segment(p::Polygon, i::Int64)
    @assert i <= p.np "Segment input too large"
    if i == p.np
        return p.pts[:,i], p.pts[:,1]
    end
    return p.pts[:,i], p.pts[:,i+1]
end

function is_visible(map::Map, p1::Vector{Float64}, p2::Vector{Float64})
    for i = 1:n_buildings(map)
        is_visible(map.b[i], p1, p2) ? (nothing) : (return false)
    end
    return true
end

function is_visible(p::Polygon, p1::Vector{Float64}, p2::Vector{Float64})
    for i = 1:n_segments(p)
        q1, q2 = get_segment(p, i)  
        segments_intersect(q1, q2, p1, p2) ? (return false) : (nothing)
    end
    return true
end
