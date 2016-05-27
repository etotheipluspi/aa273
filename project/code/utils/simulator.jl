
function simulate(m, f, nsteps::Int64, x0::Vector{Float64})
    xh = zeros(m.xdim, nsteps)
    yh = zeros(m.ydim, nsteps)
    mh = zeros(m.xdim, nsteps)
    sh = zeros(m.xdim, m.xdim, nsteps)

    x = x0
    cstep = 1
    t = 0.0
    while cstep <= nsteps
        # move system forward and make obs
        y = observe(m, x) 
        xp = step(m, x) 

        # record trajectories
        xh[:,cstep] = x
        yh[:,cstep] = y
        mh[:,cstep] = f.mu_pred
        sh[:,:,cstep] = f.sig_pred

        # apply Kalman filter
        predict!(f, m)
        update!(f, m, y)

        # to next time step
        x = xp
        cstep += 1
        t += m.dt
    end
    return xh, yh, mh, sh
end 


function simulate_no_noise(m, nsteps::Int64, x0::Vector{Float64})
    xh = zeros(m.xdim, nsteps)

    x = x0
    cstep = 1
    t = 0.0
    while cstep <= nsteps
        # move system forward and make obs
        y = observe(m, x)
        xp = step(m, x) 

        # record trajectories
        xh[:,cstep] = x

        # to next time step
        x = xp
        cstep += 1
        t += m.dt
    end
    return xh
end 


function simulate_map(m, f, nsteps::Int64, x0::Vector{Float64})
    xh = zeros(m.xdim, nsteps)
    yh = zeros(m.ydim, nsteps)
    mh = zeros(m.xdim, nsteps)
    ih = zeros(3, nsteps)
    sh = zeros(m.xdim, m.xdim, nsteps)

    x = x0
    cstep = 1
    t = 0.0
    while cstep <= nsteps
        # move system forward and make obs
        u = control(m, x, t)
        y = observe(m, x) + obs_noise(m)
        xp = step(m, x, u) + proc_noise(m)
        #y = observe(m, x) 
        #xp = step(m, x, u) 

        # record trajectories
        xh[:,cstep] = x
        yh[:,cstep] = y
        mh[:,cstep] = f.mu_pred
        sh[:,:,cstep] = f.sig_pred
        ih[:,cstep] = m.internal_state

        # apply Kalman filter
        predict!(f, m, u)
        update!(f, m, y)

        # to next time step
        x = xp
        cstep += 1
        t += m.dt
    end
    return xh, yh, mh, sh, ih
end 

function simulate_targets(m, f, nsteps::Int64, x0::Vector{Float64})
    xh = zeros(m.xdim, nsteps)
    yh = zeros(m.ydim, nsteps)
    mh1 = zeros(m.xdim, nsteps)
    mh2 = zeros(m.xdim, nsteps)
    sh1 = zeros(m.xdim, m.xdim, nsteps)
    sh2 = zeros(m.xdim, m.xdim, nsteps)
    wh = zeros(4, nsteps)

    x = x0
    cstep = 1
    t = 0.0
    while cstep <= nsteps
        # move system forward and make obs
        u = control(m, x, t)
        y = observe(m, x) + obs_noise(m)
        xp = step(m, x, u) + proc_noise(m)
        #y = observe(m, x) 
        #xp = step(m, x, u) 

        # record trajectories
        xh[:,cstep] = x
        yh[:,cstep] = y
        mh1[:,cstep] = f.mu1
        mh2[:,cstep] = f.mu2
        sh1[:,:,cstep] = f.sig1
        sh2[:,:,cstep] = f.sig2
        wh[:,cstep] = f.weights

        # apply Kalman filter
        predict!(f, m, u)
        update!(f, m, y)

        # to next time step
        x = xp
        cstep += 1
        t += m.dt
    end
    return xh, yh, mh1, mh2, sh1, sh2, wh
end 
