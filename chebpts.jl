VER = v"0.7.0"
if VERSION >= VER
    error("Julia version < "*string(VER)*" required!")
end

# Handles Chebyshev points as well as the square to quadrilateral transformation

include("quadmesh.jl")
include("cheb_convert.jl")
include("basic_chebmat.jl")

# computes the area of a polygon via the Shoelace Theorem
# each row of the polygon matrix has an x coordinate and a y coordinate
function area( polygon::Matrix{Float64} )
    polygon = convert(Array{Float64}, polygon)
    n = size(polygon, 1);
    mod1 = vcat(2:n, 1);
    a = sum(polygon[:, 1].*polygon[mod1, 2]);
    a = a - sum(polygon[mod1, 1].*polygon[:, 2]);
    a = a*0.5;
    a
end

# Chebyshev points on an interval
function chebpts(n::Int64, interval::Vector{Float64} = [-1.0, 1.0])
    m = n-1;
    chpts = sin(pi*(-m:2:m)./(2*m));
    chpts.*=(0.5*diff(interval));
    chpts.+=(0.5*sum(interval));
    chpts
end

# generates a grid from a vector of x values and a vector of y values
function meshgrid(x::Vector{Float64}, y::Vector{Float64} = x)
    xx = ones(size(y))*x.';
    yy = y*ones(size(x)).';
    xx, yy
end

# A triangle can be broken into three quadrilaterals using the medians of the triangle.
# coords gives the coordinates of the quadrilateral with vertex m of the triangle
function coords(tri::Matrix{Float64}, m::Int64)
    v = mod(m-1:m+1, 3) + 1;
    c = hcat(tri[v[1], 1:2],
            (tri[v[1], 1:2] + tri[v[2], 1:2])/2,
            (tri[v[1], 1:2] + tri[v[2], 1:2] + tri[v[3], 1:2])/3,
            (tri[v[1], 1:2] + tri[v[3], 1:2])/2).';
    c
end

# Converts a quadrilateral to a transformation
# x = a1 + b1*r + c1*s + d1*r*s
# y = a2 + b2*r + c2*s + d2*r*s
# x and y are coordinates on the quadrilateral and r and s are horizontal and vertical coordinates on the square [-1, 1] x [-1, 1]
function quad_to_params(quad::Array{Float64})
    coords_to_params = [ 1  1  1  1;
                        1 -1 -1  1;
                        1  1 -1 -1;
                        1 -1  1 -1] / 4;
    return coords_to_params*quad;
end

# The determinant of the Jacobian of the transformation has the form rcoeff*r + scoeff*s + ccoeff
function det_coeffs(quad::Matrix{Float64})
    params = quad_to_params(quad);
    rcoeff = params[2, 1]*params[4, 2] - params[2, 2]*params[4, 1];
    scoeff = params[3, 2]*params[4, 1] - params[3, 1]*params[4, 2];
    ccoeff = params[2, 1]*params[3, 2] - params[2, 2]*params[3, 1];
    return [rcoeff, scoeff, ccoeff];
end

# Generates a grid of scaled chebyshev points for every quadrilateral in a quadrilateral mesh
function cheb_mesh_grid( mesh::QuadMesh, n::Int64 )
    coords_to_params = [ 1  1  1  1;
                        1 -1 -1  1;
                        1  1 -1 -1;
                        1 -1  1 -1] / 4;

    xx, yy = meshgrid(chebpts(n));

    qsize = size(mesh.t, 1);
    quads = Array{Array{Float64, 2}, 1}(qsize);
    prms = Array{Array{Float64, 2}, 1}(qsize);
    x = Array{Array{Float64, 2}, 1}(qsize);
    y = Array{Array{Float64, 2}, 1}(qsize);

    for i = 1:qsize
        quads[i] = mesh.p[mesh.t[i, :], :];
    end

    for i = 1:qsize
        params = coords_to_params * quads[i];
        prms[i] = params;

        # x, y on quadrilateral
        x[i] = params[1, 1] + xx.*params[2, 1] + yy.*params[3, 1] + xx.*yy.*params[4, 1];
        y[i] = params[1, 2] + xx.*params[2, 2] + yy.*params[3, 2] + xx.*yy.*params[4, 2];
    end
    x, y, quads, prms;
end
