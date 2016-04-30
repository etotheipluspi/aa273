include("satellite.jl")

type EKF
    mu_up::Vector{Float64}
    sig_up::Matrix{Float64}
    mu_pred::Vector{Float64}
    sig_pred::Matrix{Float64}
end


function update!(f::EKF, s::Satellite, y::Vector{Float64})
    mu = f.mu_pred
    sig = f.sig_pred
    C = compute_c(s, mu)

    # update step
    K = sig * C * inv(C' * sig * C + s.R)
    f.mu_up = mu + K * (y - observe(s, mu))
    f.sig_up = sig - K * C * sig 
end


function predict!(f::EKF, s::Satellite)
    mu = f.mu_up
    sig = f.sig_up
    A = compute_a(s, mu)

    # predict step
    f.mu_pred = step(s, mu)
    f.sig_pred = A * sig * A' + s.Q
end
