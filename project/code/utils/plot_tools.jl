using PyPlot

function plot_map(m::Map)
    for i = 1:n_buildings(m)
        plot_build(m.b[i])
    end
end

function plot_build(b::Building)
    for i = 1:n_segments(b)
        p1, p2 = get_segment(b, i) 
        plot([p1[1], p2[1]], [p1[2], p2[2]], "k", lw=1.5)
    end
end
