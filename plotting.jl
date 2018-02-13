# Plots using PyPlot

if !isdefined(:PyPlot)
    using PyPlot
end

function pcolor_coeffs( x, y, u, lim = 0, ax = 0, Cmap = "jet" )
    if (typeof(u) == Vector{Array{Float64}})||(typeof(u) == Vector{Array{Float64, 2}})
        sz = size(u, 1);
        pl = cheb_plan(size(u[1], 1));
        uu = Vector{Array{Float64, 2}}(sz);

        for i = 1:sz
        uu[i] = coeffs2vals(u[i], pl);
        end

        if lim == 0
            for i = 1:sz
                fm,_ = findmax(abs(uu[i]));
                if fm > lim
                    lim = fm;
                end
            end
        end

        plt[:ion]()
        for i = 1:sz
            pcolormesh(x[i], y[i], uu[i], vmin = -lim, vmax = lim, cmap = ColorMap(Cmap), shading = "gouraud");
        end
        plt[:ioff]()
    elseif (typeof(u) == Array{Array{Float64}, 2})||(typeof(u) == Array{Array{Float64, 2}, 2})
        sz = size(u);
        pl = cheb_plan(size(u[1], 1));
        uu = Array{Array{Float64, 2}, 2}(sz);

        for i = 1:sz[1]
            for j = 1:sz[2]
                uu[i, j] = coeffs2vals(u[i, j], pl);
            end
        end

        if lim == 0
            for i = 1:sz[1]
                for j = 1:sz[2]
                    fm,_ = findmax(abs(uu[i, j]));
                    if fm > lim
                        lim = fm;
                    end
                end
            end
        end

        plt[:ion]()
        for i = 1:sz[1]
            for j = 1:sz[2]
                pcolormesh(x[i, j], y[i, j], uu[i, j], vmin = -lim, vmax = lim, cmap = ColorMap(Cmap), shading = "gouraud");
            end
        end
        plt[:ioff]()
    else
        pl = cheb_plan(size(u, 1));
        uu = coeffs2vals(u, pl);
        if lim == 0
            lim,_ = findmax(abs(uu));
        end
        plt[:ion]()
        pcolormesh(x, y, uu, vmin = -lim, vmax = lim, cmap = ColorMap(Cmap), shading = "gouraud");
        plt[:ioff]()
    end
    axis("equal");
    if ax != 0;
        axis(ax);
    end
end

function surf_coeffs( x, y, u, lim, Cmap = "jet" )
    if (typeof(u) == Vector{Array{Float64}})||(typeof(u) == Vector{Array{Float64, 2}})
        sz = size(u, 1);
        pl = cheb_plan(size(u[1], 1));
        uu = Vector{Array{Float64, 2}}(sz);
        plt[:ion]()
        for i = 1:sz
            uu[i] = coeffs2vals(u[i], pl);
            plot_surface(x[i], y[i], uu[i], vmin = -lim, vmax = lim, cmap = ColorMap(Cmap), cstride = 1, rstride = 1);
        end
        plt[:ioff]()
    elseif (typeof(u) == Array{Array{Float64}, 2})||(typeof(u) == Array{Array{Float64, 2}, 2})
        sz = size(u);
        pl = cheb_plan(size(u[1]));
        uu = Array{Array{Float64, 2}, 2}(sz);
        plt[:ion]()
        for i = 1:sz[1]
            for j = 1:sz[2]
                uu[i, j] = coeffs2vals(u[i, j], pl);
                plot_surface(x[i, j], y[i, j], uu[i, j], vmin = -lim, vmax = lim, cmap = ColorMap(Cmap), cstride = 1, rstride = 1);
            end
        end
        plt[:ioff]()
    else
        pl = cheb_plan(size(u));
        uu = coeffs2vals(u, pl);
        plt[:ion]()
        plot_surface(x, y, uu, vmin = -lim, vmax = lim, cmap = ColorMap(Cmap), cstride = 1, rstride = 1);
        plt[:ioff]()
        
    end
    uu
end

function spy_big(A)
    m = A.m;
    n = A.n;
    colptr = A.colptr;
    vals = Array{Float64}(size(A.rowval, 1), 2);
    vals[:, 2] = A.rowval;
    mn = max(m, n);
    if mn > 2000
        mn /= 2000;
    else
        mn = 1;
    end

    for i = 1:n
        vals[colptr[i] : colptr[i+1]-1, 1] = floor(i/mn+0.5)*mn;
    end
    for i = 1:size(vals, 1)
        vals[i, 2] = floor(vals[i, 2]/mn+0.5)*mn;
    end
    vals = unique(vals, 1);
    plt[:ion]()
    scatter(vals[:, 1], vals[:, 2], color = "b", marker = ".", s = 1);
    axis(axis("equal")[[1, 2, 4, 3]])
    plt[:ioff]()
end
