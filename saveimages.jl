VER = v"0.7.0"
if VERSION >= VER
    error("Julia version < "*string(VER)*" required!")
end

# Saves hi-res images of Navier-Stokes solutions
using Images
using Colors
using ColorTypes
include("quadmesh.jl")

function mul2!(a::Array{Float64}, b::Vector{Float64}, bb::Vector{Float64})
    bb[1] = a[1, 1]*b[1] + a[1, 2]*b[2];
    bb[2] = a[2, 1]*b[1] + a[2, 2]*b[2];
end

function inv2s(a::Array{Float64}, b::Vector{Float64})
    bb = Vector{Float64}(2);
    invdet = (a[2, 2]*a[1, 1] - a[1, 2]*a[2, 1]);
    bb[1] = a[2, 2]*b[1] - a[1, 2]*b[2];
    bb[2] = a[1, 1]*b[2] - a[2, 1]*b[1];
    bb ./= invdet;
    bb
end

function jjet(x::Float64, interval::Vector{Float64} = [-1.0, 1.0])
    x -= interval[1];
    x /= interval[2] - interval[1];
    if x < 0
        x = 0;
    end
    if x > 1
        x = 1;
    end
    n=4*x
    r=min(max(min(n-1.5,-n+4.5),0),1);
    g=min(max(min(n-0.5,-n+3.5),0),1);
    b=min(max(min(n+0.5,-n+2.5),0),1);
    [r, g, b];
end

function chebyshevseries(x::Vector{Float64}, n::Int64)
    vals = zeros(n, length(x));
    vals[1, :] = 1;
    vals[2, :] = x;
    for i = 2:n-1
        vals[i+1, :] = 2.*x.*vals[i, :] - vals[i-1, :];
    end
    return vals
end

function chebyshevseries(x::Float64, n::Float64)
    vals = zeros(n);
    vals[1, :] = 1;
    vals[2, :] = x;
    for i = 2:n-1
        vals[i+1, :] = 2.*x.*vals[i, :] - vals[i-1, :];
    end
    return vals
end

function transform(params::Matrix{Float64}, x::Vector{Float64})
    return [params[1, 1] + params[2, 1]*x[1] + params[3, 1]*x[2] + params[4, 1]*x[1]*x[2],
            params[1, 2] + params[2, 2]*x[1] + params[3, 2]*x[2] + params[4, 2]*x[1]*x[2]]
end

function jacobian(params::Matrix{Float64}, x::Vector{Float64})
    return [params[2, 1] + params[4, 1]*x[1]  params[3, 1] + params[4, 1]*x[2];
            params[2, 2] + params[4, 2]*x[1]  params[3, 2] + params[4, 2]*x[2]];
end

function invtransform(params::Matrix{Float64}, xx::Vector{Float64})
    init = zeros(2);
    n = 0;
    for i = 1:5
        init -= inv2s(jacobian(params, init), transform(params, init)-xx)
    end
    init
end

function valueat(u::Matrix{Float64}, params::Matrix{Float64}, xx::Vector{Float64})
    xx = invtransform(params, xx);
    n = size(u, 1);
    chs = chebyshevseries(xx, n);
    val = 0.0;
    for i = 1:n
        for j = 1:n
            val += chs[i, 1]*chs[j, 2]*u[j, i];
        end
    end
    val
end

function pts_in_quad(quad::Matrix{Float64}, res::Float64, offset::Vector{Float64} = [0.0, 0.0])
    minval, minpt = findmin(quad[:, 2]);
    maxval, maxpt = findmax(quad[:, 2]);
    if quad[mod(minpt, 4)+1, 2] == minval
        minpt = mod(minpt, 4)+1;
    end
    if quad[mod(maxpt, 4)+1, 2] == maxval
        maxpt = mod(maxpt, 4)+1;
    end
    minindex = ceil(Int64, (minval-offset[2])/res);
    maxindex = floor(Int64, (maxval-offset[2])/res);
    yvalsrange = minindex:maxindex;
    xvalsrange = zeros(Int64, size(yvalsrange, 1), 2);
    edge = mod(maxpt-1:maxpt, 4)+1;
    slope = (quad[edge[2], 1] - quad[edge[1], 1])/(quad[edge[2], 2] - quad[edge[1], 2]);
    for i = size(yvalsrange, 1):-1:1
        if quad[edge[2], 2] > yvalsrange[i]*res + offset[2];
            edge = mod(edge, 4)+1;
            slope = (quad[edge[2], 1] - quad[edge[1], 1])/(quad[edge[2], 2] - quad[edge[1], 2]);
        end
        if !isfinite(slope);
            edge = mod(edge-2, 4)+1;
            slope = (quad[edge[2], 1] - quad[edge[1], 1])/(quad[edge[2], 2] - quad[edge[1], 2]);
        end
#         x = (y - y1)/m + x1;
#         here, m is 1/m
        xvalsrange[i, 1] = ceil(Int64, (quad[edge[1], 1] + slope*(yvalsrange[i]*res + offset[2] - quad[edge[1], 2])-offset[1])/res);
    end
    edge = mod(minpt-1:minpt, 4)+1;
    slope = (quad[edge[2], 1] - quad[edge[1], 1])/(quad[edge[2], 2] - quad[edge[1], 2]);
    for i = 1:size(yvalsrange, 1)
        if quad[edge[2], 2] < yvalsrange[i]*res + offset[2];
            edge = mod(edge, 4)+1;
            slope = (quad[edge[2], 1] - quad[edge[1], 1])/(quad[edge[2], 2] - quad[edge[1], 2]);
        end
        if !isfinite(slope);
            edge = mod(edge-2, 4)+1;
            slope = (quad[edge[2], 1] - quad[edge[1], 1])/(quad[edge[2], 2] - quad[edge[1], 2]);
        end
#         x = (y - y1)/m + x1;
#         here, m is 1/m
        xvalsrange[i, 2] = floor(Int64, (quad[edge[1], 1] + slope*(yvalsrange[i]*res + offset[2] - quad[edge[1], 2])-offset[1])/res);
    end
    return xvalsrange, yvalsrange;
end

bounding_box(xlist, ylist) = [findmin(xlist)[1],
                            findmax(xlist)[1],
                            findmin(ylist)[1],
                            findmax(ylist)[1]];

function bounding_box_overlap(bbox1::Vector{Float64}, bbox2::Vector{Float64})
    if any((bbox1[1:2] .>= bbox2[1]) & (bbox1[1:2] .<= bbox2[2])) && any((bbox1[3:4] .>= bbox2[3]) & (bbox1[3:4] .<= bbox2[4]))
        return true
    end
    if any((bbox2[1:2] .>= bbox1[1]) & (bbox2[1:2] .<= bbox1[2])) && any((bbox2[3:4] .>= bbox1[3]) & (bbox2[3:4] .<= bbox1[4]))
        return true
    end
    return false
end

function createimage(mesh::QuadMesh, u::Vector{Matrix{Float64}}, res::Float64, lim::Float64 = 1.0, bbox::Vector{Float64} = bounding_box(mesh.p[:, 1], mesh.p[:, 2]))
    dims = floor(Int64, [bbox[2] - bbox[1], bbox[4] - bbox[3]]./res);
    dims+=1;
    qnum = size(m.t, 1);
    vals = zeros(dims[2], dims[1]);
    filled = falses(dims[2], dims[1]);
    offset = [bbox[1], bbox[3]];
    for i = 1:qnum
        xlist, ylist = pts_in_quad(m.p[m.t[i, :], :], res, offset);
        bbox2 = bounding_box(res*xlist + offset[1], res*ylist + offset[2]);
        xlist+=1;
        ylist+=1;
        if bounding_box_overlap(bbox, bbox2)
            params = quad_to_params(m.p[m.t[i, :], :]);
            for j = 1:size(ylist, 1)
                for k = xlist[j, 1]:xlist[j, 2]
                    coords = [(k-1)*res + offset[1], (ylist[j]-1)*res + offset[2]];
                    if ylist[j] > 0 && ylist[j] <= dims[2] && k > 0 && k <= dims[1]
                        vals[ylist[j], k] = valueat(u[i], params, coords);
                        filled[ylist[j], k] = true;
                    end
                end
            end
        end
    end
    img = zeros(dims[2], dims[1], 3);
    for i = 1:dims[2]
        for j = 1:dims[1]
            if filled[i, j]
                img[dims[2]+1-i, j, :] = jjet(vals[i, j], [-lim, lim]);
            end
        end
    end
    img
end

function saveimage(path, img)
  imgm = colorim(img);
  save(path, imgm);
end
