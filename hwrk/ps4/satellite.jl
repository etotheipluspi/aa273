using Distributions

"""
    Satellite data container that stores moments of inertia and covariance matrices
"""
type Satellite
    I::Vector{Float64}
    ix::Float64
    iy::Float64
    iz::Float64
    Q::Matrix{Float64}
    R::Matrix{Float64}
    c::Float64
    dt::Float64
    W::MvNormal
    V::MvNormal
end
Satellite(I::Vector{Float64}, Q::Matrix{Float64}, R::Matrix{Float64}, c::Float64, dt::Float64) = Satellite(I, (I[2]-I[3])/I[1], (I[3]-I[1])/I[2], (I[1]-I[2])/I[3], Q, R, c, dt, MvNormal(Q), MvNormal(R))

"""
    Satellite dynamics, return the new rotational velocities
"""
function step(s::Satellite, w::Vector{Float64})
    wp = zeros(3)
    dt = s.dt
    wp[1] = w[1] + dt * s.ix * w[2] * w[3]
    wp[2] = w[2] + dt * s.iy * w[3] * w[1]
    wp[3] = w[3] + dt * s.iz * w[1] * w[2]
    return wp 
end

"""
    Satellite observation, returns a saturated observation
"""
function observe(s::Satellite, w::Vector{Float64})
    wp = [om < s.c ? (om) : (s.c) for om in w]
    return wp 
end


proc_noise(s::Satellite) = rand(s.W)
obs_noise(s::Satellite) = rand(s.V)


"""
    Computes the A matrix (dynamics Jacobian)
    Inputs:
        - s is the satellite concrete type
        - w is the rotational velocity state estimates
"""
function compute_a(s::Satellite, w::Vector{Float64})
    A = eye(3) / s.dt
    A[1,2] = s.ix * w[3] 
    A[1,3] = s.ix * w[2]
    A[2,1] = s.iy * w[3]
    A[2,3] = s.iy * w[1]
    A[3,1] = s.iz * w[2] 
    A[3,2] = s.iz * w[1] 
    A = dt * A
    return A
end


"""
    Computes the C matrix (observation Jacobian)
    Inputs:
        - s is the satellite concrete type
        - w is the rotational velocity state estimates
"""
function compute_c(s::Satellite, w::Vector{Float64})
    C = eye(3)
    w[1] < s.c ? (nothing) : (C[1,1] = 0.)
    w[2] < s.c ? (nothing) : (C[2,2] = 0.)
    w[3] < s.c ? (nothing) : (C[3,3] = 0.)
    return C
end
