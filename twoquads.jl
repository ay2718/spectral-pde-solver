# Provides matrices to stitch two quadrilaterals together

include("chebpts.jl")

function fastorder2boundary(m::QuadMesh, edge::Int64, n::Int64)
    Agi = [zeros(n, n*n), zeros(n, n*n)];
    Aig = [zeros(n*n, n), zeros(n*n, n)];

    lrownums = 1:n;
    rrownums = n+1:2*n;
    drownums = 2*n+1:n:n*(n-1)+1;
    urownums = 2*n+2:n:n*(n-1)+2;

    rownums = sort(vcat(lrownums, rrownums, drownums, urownums))

    a = m.p[m.e[edge, 2], 1] - m.p[m.e[edge, 1], 1];
    b = m.p[m.e[edge, 2], 2] - m.p[m.e[edge, 1], 2];
    r = sqrt(a^2 + b^2);
    a /= r;
    b /= r;
    nn = n-1;

#     Edges:  1 is top, 2 is left, 3 is bottom, 4 is right
#     Vertices:  1 is top right, then CCW ordering    
#     D is the normal derivative, and u_r is shorthand for du/dr (partial derivative)
#     det(r, s)*D = u_r(b y_s + a x_s) - u_s (b y_r + a x_r)
#     det(r, s)*D = u_r(q_c + q_d r) - u_s (q_b + q_d s)
#     q_b = B b_2 + A b_1, the boundary has slope B/A
#     q_c = B c_2 + A c_1
#     q_d = B d_2 + A d_1

#     if we have a boundary edge, we can't do anything
    if m.e[edge, 5] == 0
        return Aig, Agi, Agg;
    end

#     Quadrilateral coordinates and transformation parameters
    quads = [zeros(4, 2), zeros(4, 2)];
    params = [zeros(4, 2), zeros(4, 2)];

    quads[1] = m.p[m.t[m.e[edge, 3], :], :];
    quads[2] = m.p[m.t[m.e[edge, 5], :], :];

    for i = 1:2
        params[i] = quad_to_params(quads[i]);
    end

    q = zeros(3, 2);

    for i = 1:3
        for j = 1:2
            q[i, j] = a*params[j][i+1, 1] + b*params[j][i+1, 2];
        end
    end

#     u_r at (chpts(i), chpts(j))
    u_r(i::Int64, j::Int64) = kron(chviSD1[i, :], chv[j, :])
    u_s(i::Int64, j::Int64) = kron(chv[i, :], chviSD1[j, :])

    u_x(i::Int64, j::Int64, p::Matrix{Float64}, d::Matrix{Float64}) = d[i, j]*(u_r(i, j)*(p[3, 2]+chpts[i]*p[4, 2])-u_s(i, j)*(p[2, 2]+chpts[j]*p[4, 2]));

    u_y(i::Int64, j::Int64, p::Matrix{Float64}, d::Matrix{Float64}) = d[i, j]*(u_s(i, j)*(p[2, 1]+chpts[j]*p[4, 1])-u_r(i, j)*(p[3, 1]+chpts[i]*p[4, 1]));

    D1 = diffmat1(n);
    chv = cheb_vander(n);
    iS = invconvertmat0(n);
    
#     chv*iS*D1*u, where u is a Chebyshev series takes the derivative of u, converts the result back to Chebyshev coefficients, and then converts that result into values at Chebyshev points
    chviSD1 = chv*iS*D1;

    chpts = chebpts(n);
    r, s = meshgrid(chpts);
#     For each edge, we create an Agi block and an Aig block.  The Aig block handles Dirichlet conditions and Agi block handles Neumann conditions
    for i = 1:2
        dc = det_coeffs(quads[i]);
        detvals = dc[3] + dc[2]*r + dc[1]*s;
        detvals = 1./detvals;
        pts = zeros(Int64, n, 2)
#         Here, we figure out what orientation the edge is and then figure out what Chebyshev points to stitch to the other edge
        edgemod = m.e[edge, 2*i+2];
        if edgemod == 1
            pts[:, 1] = n:-1:1
            pts[:, 2] = n;
        elseif edgemod == 2
            pts[:, 1] = 1;
            pts[:, 2] = n:-1:1;
        elseif edgemod == 3
            pts[:, 1] = 1:n;
            pts[:, 2] = 1;
        else
            pts[:, 1] = n;
            pts[:, 2] = 1:n;
        end
#         In order to properly stitch together points on the boundary we must reverse the sign of the derivative on one edge.  Also, we have to reverse order the points on one edge
        if i == 1
            for j = 1:nn
                Aig[i][rownums[j+nn*(edgemod-1)], j] = -1
            end
            for j = 1:n
                Agi[i][j, :] = b*u_x(pts[j, 1], pts[j, 2], params[i], detvals) - a*u_y(pts[j, 1], pts[j, 2], params[i], detvals);
            end
        else
            for j = 1:nn
                Aig[i][rownums[j+nn*(edgemod-1)], n+1-j] = -1
            end
            for j = 1:n
                Agi[i][n+1-j, :] = -b*u_x(pts[j, 1], pts[j, 2], params[i], detvals) + a*u_y(pts[j, 1], pts[j, 2], params[i], detvals);
            end
        end
    end

    Aig, Agi
end
