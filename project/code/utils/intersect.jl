# points is 2xN matrix
function is_visible(points::Matrix{Float64}, p::Vector{Float64}, p2::Vector{Float64})
    npts = size(points, 1) 
    for i = 1:npts

    end
end

# checks if 2D line segments p and q intersect
function segments_intersect(p::Vector{Float64}, p2::Vector{Float64}, q::Vector{Float64}, q2::Vector{Float64})
    r = p2 - p
    s = q2 - q
    uNum = cross2D(q - p, r)
    denom = cross2D(r, s)
    # check if lines are collinear
    if uNum == 0 && denom == 0
        # do they touch?
        if p == q || p == p2 || p2 == q || p2 == q2
            return true
        end
        # do they overlap?
        return ((q[1] - p[1] < 0) != (q[1] - p2[1] < 0) != (q2[1] - p[1] < 0) != (q2[1] - p2[1] < 0) ||
                (q[2] - p[2] < 0) != (q[2] - p2[2] < 0) != (q2[2] - p[2] < 0) != (q2[2] - p2[2] < 0))
    end
    # lines are parallel
    if denom == 0
        return false
    end
    u = uNum / denom
    # TODO: The line below is a bottleneck
    t = cross2D(q - p, s) / denom
    return (t >= 0) && (t <= 1) && (u >= 0) && (u <= 1)
end

# compute cross product in 2D
function cross2D(l1, l2)
    return l1[1] * l2[2] - l1[2] * l2[1]
end
