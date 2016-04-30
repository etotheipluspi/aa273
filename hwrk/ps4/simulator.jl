using PyPlot

include("satellite.jl")
include("ekf.jl")


function simulate(s::Satellite, f::EKF, nsteps::Int64, x0::Vector{Float64})
    xh = zeros(3, nsteps)
    yh = zeros(3, nsteps)
    mh = zeros(3, nsteps)
    sh = zeros(3, 3, nsteps)

    x = x0
    cstep = 1
    while cstep <= nsteps
        # move system forward and make obs
        y = observe(s, x) + obs_noise(s)
        xp = step(s, x) + proc_noise(s)

        # record trajectories
        xh[:,cstep] = x
        yh[:,cstep] = y
        mh[:,cstep] = f.mu_pred
        sh[:,:,cstep] = f.sig_pred

        # apply Kalman filter
        update!(f, s, y)
        predict!(f, s)

        # to next time step
        x = xp
        cstep += 1
    end
    return xh, yh, mh, sh
end 


I = [1., 5, 5]
c = 10.
Q = 0.01 * eye(3)
R = 0.1 * eye(3)
dt = 0.01
w0 = [10., 0.1, 0.1]
mu0 = [10., 0, 0]
sig0 = eye(3)
nsteps = 1000

sat = Satellite(I, Q, R, c, dt)
filter = EKF(zeros(3), zeros(3,3), mu0, sig0)

xh, yh, mh, sh = simulate(sat, filter, nsteps, w0)

ax = 1
subplot(2,1,1)
plot(xh[ax,:]')
plot(yh[ax,:]')
plot(mh[ax,:]')
xlabel("time")
ylabel("omega")

sig = zeros(3, nsteps)

for i = 1:nsteps
    sig[1,i] = sh[1,1,i]
    sig[2,i] = sh[2,2,i]
    sig[3,i] = sh[3,3,i]
end

subplot(2,1,2)
plot(sig[ax,:]')
xlabel("time")

