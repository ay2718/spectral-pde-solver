VER = v"0.7.0"
if VERSION >= VER
    error("Julia version < "*string(VER)*" required!")
end

# Handles solving PDEs on ONE quadrilateral

include("quadmesh.jl")

# Produces a large matrix that multiplies a vector of 2d ultraspherical coefficients of parameter 2 by some function (given in Chebyshev coefficients)
function kronmultmatrix( coeffs::Matrix{Float64}, n::Int64, cutoff::Float64 = 1e-13 )
    mult_matrices = Vector{SparseMatrixCSC{Float64, Int64}}(n);
    mult_matrices[1] = speye(n);
    mult_matrices[2] = multmat2(n);
    for i = 3:n
        mult_matrices[i] = mult_matrices[2]*mult_matrices[i-1] - mult_matrices[i-2];
    end
    finalkronmatrix = spzeros(n*n, n*n);
    for i = 1:n
        ycoeffstotal = spzeros(n, n);
        for j = 1:n
            if abs(coeffs[i, j]) > cutoff
                ycoeffstotal+=coeffs[i, j]*mult_matrices[j];
            end
        end
        finalkronmatrix += kron(ycoeffstotal, mult_matrices[i]);
    end
    return finalkronmatrix
end

# Creates a matrix to solve a PDE with constant coefficients [u_xx, u_xy, u_yy, u_x, u_y, u] on a quadrilateral with vertices coords and polynomial degree n-1
# Rows corresponding to components of rhs with polynomial degreen n-1 or n-2 are eliminated to be replaced with boundary condition rows
function fastorder2matrix( coords::Matrix{Float64}, n::Int64, coeffs::Vector{Float64} )

    D1n = diffmat1(n);
    D2n = diffmat2(n);

#     Conversion matrices
    S0n = convertmat0(n);
    S1n = convertmat1(n);
    Sn = S1n*S0n;
#     Sn = S1n * S0n;

#     r, s on square
    (r, s) = meshgrid(chebpts(n));
    
    coords_to_params = [ 1  1  1  1;
                        1 -1 -1  1;
                        1  1 -1 -1;
                        1 -1  1 -1] / 4;

    params = coords_to_params * coords;
    
    a1 = params[1, 1];
    a2 = params[1, 2];
    
    b1 = params[2, 1];
    b2 = params[2, 2];
    
    c1 = params[3, 1];
    c2 = params[3, 2];
    
    d1 = params[4, 1];
    d2 = params[4, 2];
    
    # x, y on quadrilateral
    x = a1 + r.*b1 + s.*c1 + r.*s.*d1;
    y = a2 + r.*b2 + s.*c2 + r.*s.*d2;
    
    # Multiplication matrices for C2
    M2xn = multmat2(n);
    
    In = speye( n );
    Inn = speye( n*n );
#     dr/dx here is multiplied the determinant of the Jacobian (det)
    drdx = c2.*Inn + d2.*kron(M2xn, In);
    drdy = -(c1.*Inn + d1.*kron(M2xn, In));
    dsdx = -(b2.*Inn + d2.*kron(In, M2xn));
    dsdy = b1.*Inn + d1.*kron(In, M2xn);
    
    ddetdr = b1*d2 - b2*d1;
    ddetds = c2*d1 - c1*d2;
    
    det = (b1*c2 - b2*c1)*Inn + ddetdr*kron(M2xn, In) + ddetds*kron(In, M2xn);
    
#     d/dx det(r, s) is multiplied here by the det(r, s)
    ddetdx = ddetdr*drdx + ddetds*dsdx;
    ddetdy = ddetdr*drdy + ddetds*dsdy;
    
#     d^2 r/dx^2 is multiplied here by det(r, s)^3
    d2rdx =  d2*drdx*det - drdx*ddetdx;
    d2sdx = -d2*dsdx*det - dsdx*ddetdx;
    d2rdy = -d1*drdy*det - drdy*ddetdy;
    d2sdy =  d1*dsdy*det - dsdy*ddetdy;
    
    d2rdxy =  d1*drdx*det - drdy*ddetdx;
    d2sdxy = -d1*dsdx*det - dsdy*ddetdx;
    
#     operators for taking derivatives on the square
    urr = kron(D2n, Sn);
    urs = kron(S1n*D1n, S1n*D1n);
    uss = kron(Sn, D2n);
    ur = kron(S1n*D1n, Sn);
    us = kron(Sn, S1n*D1n);

#     operators for taking derivatives on the quadrilateral (and then multiplying everything by det(r, s)^3)
    uxx = det*(drdx^2)*urr +
            2*det*(drdx*dsdx)*urs +
            det*(dsdx^2)*uss +
            (d2rdx)*ur +
            (d2sdx)*us;

    uyy = det*(drdy^2)*urr +
            2*det*(drdy*dsdy)*urs +
            det*(dsdy^2)*uss +
            (d2rdy)*ur +
            (d2sdy)*us;

    uxy = det*(drdx*drdy)*urr +
            det*(drdx*dsdy + drdy*dsdx)*urs +
            det*(dsdx*dsdy)*uss +
            (d2rdxy)*ur +
            (d2sdxy)*us;

    ux = (det^2)*drdx*ur +
        (det^2)*dsdx*us;

    uy = (det^2)*drdy*ur +
        (det^2)*dsdy*us;

    uu = (det^3)*kron(Sn, Sn);

    A = coeffs[1]*uxx + coeffs[2]*uxy + coeffs[3]*uyy + coeffs[4]*ux + coeffs[5]*uy + coeffs[6]*uu;

    bulk_rows = zeros(n-2, n-2);

    for i = 1:n-2
        for j = 1:n-2
            bulk_rows[i, j] = n*(j-1) + i;
        end
    end

    bulk_rows = bulk_rows[:];
    bulk_rows = convert(Array{Int64}, bulk_rows);
    A = A[bulk_rows, :];
    precondition, _ = findmax(abs(A), 2);
    precondition = 1./precondition[:];
    precondition = Diagonal(precondition);
    A = precondition*A;

    A, precondition, x, y
end

# Same as above, but with variable coefficients (only one line is different)
function fastorder2matrix( coords::Matrix{Float64}, n::Int64, coeffs::Vector{Matrix{Float64}} )
    D1n = diffmat1(n);
    D2n = diffmat2(n);

    # Conversion matrices
    S0n = convertmat0(n);
    S1n = convertmat1(n);
    Sn = S1n*S0n;
    # Sn = S1n * S0n;

    # r, s on square
    (r, s) = meshgrid(chebpts(n));

    coords_to_params = [ 1  1  1  1;
                        1 -1 -1  1;
                        1  1 -1 -1;
                        1 -1  1 -1] / 4;

    params = coords_to_params * coords;

    a1 = params[1, 1];
    a2 = params[1, 2];

    b1 = params[2, 1];
    b2 = params[2, 2];

    c1 = params[3, 1];
    c2 = params[3, 2];

    d1 = params[4, 1];
    d2 = params[4, 2];

    # x, y on quadrilateral
    x = a1 + r.*b1 + s.*c1 + r.*s.*d1;
    y = a2 + r.*b2 + s.*c2 + r.*s.*d2;

    # Multiplication matrices for C2
    # C2_0 = 1, C2_1 = 4x
    M2xn = multmat2(n);

    In = speye( n );
    Inn = speye( n*n );
    drdx = c2.*Inn + d2.*kron(M2xn, In);
    drdy = -(c1.*Inn + d1.*kron(M2xn, In));
    dsdx = -(b2.*Inn + d2.*kron(In, M2xn));
    dsdy = b1.*Inn + d1.*kron(In, M2xn);

    ddetdr = b1*d2 - b2*d1;
    ddetds = c2*d1 - c1*d2;

    det = (b1*c2 - b2*c1)*Inn + ddetdr*kron(M2xn, In) + ddetds*kron(In, M2xn);

    ddetdx = ddetdr*drdx + ddetds*dsdx;
    ddetdy = ddetdr*drdy + ddetds*dsdy;

    d2rdx =  d2*drdx*det - drdx*ddetdx;
    d2sdx = -d2*dsdx*det - dsdx*ddetdx;
    d2rdy = -d1*drdy*det - drdy*ddetdy;
    d2sdy =  d1*dsdy*det - dsdy*ddetdy;

    d2rdxy =  d1*drdx*det - drdy*ddetdx;
    d2sdxy = -d1*dsdx*det - dsdy*ddetdx;

    urr = kron(D2n, Sn);
    urs = kron(S1n*D1n, S1n*D1n);
    uss = kron(Sn, D2n);
    ur = kron(S1n*D1n, Sn);
    us = kron(Sn, S1n*D1n);

    uxx = det*(drdx^2)*urr +
            2*det*(drdx*dsdx)*urs +
            det*(dsdx^2)*uss +
            (d2rdx)*ur +
            (d2sdx)*us;

    uyy = det*(drdy^2)*urr +
            2*det*(drdy*dsdy)*urs +
            det*(dsdy^2)*uss +
            (d2rdy)*ur +
            (d2sdy)*us;

    uxy = det*(drdx*drdy)*urr +
            det*(drdx*dsdy + drdy*dsdx)*urs +
            det*(dsdx*dsdy)*uss +
            (d2rdxy)*ur +
            (d2sdxy)*us;

    ux = (det^2)*drdx*ur +
        (det^2)*dsdx*us;

    uy = (det^2)*drdy*ur +
        (det^2)*dsdy*us;

    uu = (det^3)*kron(Sn, Sn);
    
#     The only line different from the above fastorder2matrix
    A = kronmultmatrix(coeffs[1], n)*uxx + kronmultmatrix(coeffs[2], n)*uxy + kronmultmatrix(coeffs[3], n)*uyy + kronmultmatrix(coeffs[4], n)*ux + kronmultmatrix(coeffs[5], n)*uy + kronmultmatrix(coeffs[6], n)*uu;

    bulk_rows = zeros(n-2, n-2);

    for i = 1:n-2
        for j = 1:n-2
            bulk_rows[i, j] = n*(j-1) + i;
        end
    end

    bulk_rows = bulk_rows[:];
    bulk_rows = convert(Array{Int64}, bulk_rows);
    A = A[bulk_rows, :];
    precondition, _ = findmax(abs(A), 2);
    precondition = 1./precondition[:];
    precondition = Diagonal(precondition);
    A = precondition*A;

    A, precondition, x, y
end

# Boundary condition rows for the matrix created by fastorder2matrix
# fastorder2ptsbcs can handle Dirichlet boundary conditions, Neumann bounday conditions, and Neumann boundary conditions that pretend that the normal vector changes along the boundary to avoid sharp changes at corners
function fastorder2ptsbcs(m::QuadMesh, t::Int64, n::Int64, neumann_list::BitVector, rounded::BitVector = falses(4))
#     List of pre-calculated stuff that we'll need later
#     Quadrilateral coordinates
    coords = m.p[m.t[t, :], :];
    chv = cheb_vander(n);
    iS0 = invconvertmat0(n);
    D1 = diffmat1(n);
    chviSD1 = chv*iS0*D1;
    chpts = chebpts(n);
    nn = n-1;
    leftrows(i) = i:n:n*n - n + i;
    downrows(i) = (n*i - n + 1):n*i;
#     BCS will eventually become the block of boundary conditions
    BCS = zeros(4*n-4, n*n);
#     Summing nbrow.*u (u is a Chebyshev series) gives u(-1)
    nbrow = (-1).^(0:n-1);
#     Summing u (u is a Chebyshev series) gives u(1)
    pbrow = ones(n);
#     Summing pdrow.*u gives u'(1)
    pdrow = (0:n-1).^2;
#     Summing ndrow.*u gives u'(-1)
    ndrow = -nbrow.*pdrow;
    r, s = meshgrid(chpts);
#     Determinant coefficients
    dc = det_coeffs(coords);
    detvals = dc[3] + dc[2]*r + dc[1]*s;
    detvals = 1./detvals;
    params = quad_to_params(coords);
#     u*v gives v(chpts(i), chpts(j))
    u(i::Int64, j::Int64) = kron(chv[i, :], chv[j, :])
#     u_r*v gives v_r(chpts(i), chpts(j)) (dv/dr)
    u_r(i::Int64, j::Int64) = kron(chviSD1[i, :], chv[j, :])
#     u_s*v gives v_s(chpts(i), chpts(j)) (dv/ds)
    u_s(i::Int64, j::Int64) = kron(chv[i, :], chviSD1[j, :])
#     u_x*v gives v_x(chpts(i), chpts(j)), d is the determinant, p is the list of parameters [a1 a2; b1 b2; c1 c2; d1 d2]
    u_x(i::Int64, j::Int64, p::Matrix{Float64}, d::Matrix{Float64}) = d[i, j]*(u_r(i, j)*(p[3, 2]+chpts[i]*p[4, 2])-u_s(i, j)*(p[2, 2]+chpts[j]*p[4, 2]));
#     u_y*v gives v_y(chpts(i), chpts(j))
    u_y(i::Int64, j::Int64, p::Matrix{Float64}, d::Matrix{Float64}) = d[i, j]*(u_s(i, j)*(p[2, 1]+chpts[j]*p[4, 1])-u_r(i, j)*(p[3, 1]+chpts[i]*p[4, 1]));

    for i = 1:4
#         The first n-1 rows are for the top edge, then the next n-1 rows are for the left edge, the next n-1 rows are for the bottom edge, and the final n-1 rows are for the right edge
        rownums = nn*(i-1)+(1:nn);
        pts = zeros(Int64, nn, 2)
#         Points go counterclockwise
        if i == 1
            pts[:, 1] = n:-1:2
            pts[:, 2] = n;
        elseif i == 2
            pts[:, 1] = 1;
            pts[:, 2] = n:-1:2;
        elseif i == 3
            pts[:, 1] = 1:nn;
            pts[:, 2] = 1;
        else
            pts[:, 1] = n;
            pts[:, 2] = 1:nn;
        end
#         If the edge given has Neumann conditions, figure out the normal vector and set the derivative across it
        if neumann_list[i]
            a = coords[mod(i, 4)+1, :] - coords[i, :];
            r = sqrt(a[1]^2 + a[2]^2);
            a /= r;
        
            en = m.e[next_edge(m, m.te[t, i]), 1:2];
            an = m.p[en[2], :] - m.p[en[1], :];
            r = sqrt(an[1]^2 + an[2]^2);
            an /= r;
    
            ep = m.e[prev_edge(m, m.te[t, i]), 1:2];
            ap = m.p[ep[2], :] - m.p[ep[1], :];
            r = sqrt(ap[1]^2 + ap[2]^2);
            ap /= r;
#             If the edge is rounded, then gradually change the normal vector over the edge
            if rounded[i]
                aa = zeros(nn, 2);
                for j = 1:div(nn, 2)
                    ratio = -0.5*chpts[i];
                    aa[j, :] = a*(1-ratio)+an*ratio;
                end
                for j = div(nn, 2)+1:nn
                    ratio = 0.5*chpts[i];
                    aa[j, :] = a*(1-ratio)+ap*ratio;
                end
                for j = 1:nn
                    aa[j, :] /= sqrt(aa[j, 1]^2 + aa[j, 2]^2);
                    BCS[rownums[j], :] = aa[j, 2]*u_x(pts[j, 1], pts[j, 2], params, detvals) - aa[j, 1]*u_y(pts[j, 1], pts[j, 2], params, detvals);
                end
#             If it is not rounded, then just keep the normal vector as is
            else
                for j = 1:nn
                    BCS[rownums[j], :] = a[2]*u_x(pts[j, 1], pts[j, 2], params, detvals) - a[1]*u_y(pts[j, 1], pts[j, 2], params, detvals);
                end
            end
#         Otherwise, the edge has Dirichlet conditions.
        else
            for j = 1:nn
                BCS[rownums[j], :] = u(pts[j, 1], pts[j, 2]);
            end
        end
    end
    BCS
end

# Multiplies u by the determinant of the jacobian of the transformation for a specific quadrilateral (defined by params)
function dets(u::Matrix{Float64}, params::Matrix{Float64}, M2xn::SparseMatrixCSC{Float64, Int64})
    ddetdxx = params[2, 1]*params[4, 2] - params[2, 2]*params[4, 1];
    ddetdyy = params[3, 2]*params[4, 1] - params[3, 1]*params[4, 2];
    ddetdzz = params[2, 1]*params[3, 2] - params[2, 2]*params[3, 1];
    ddetdzz *= 0.5;
    d = u*(ddetdxx*M2xn.' + ddetdzz*I) + (ddetdyy*M2xn + ddetdzz*I)*u;
    d
end

# Computes a right hand side (Lu = rhs) for a given function rhs in Chebyshev coefficients
function fastorder2rhs( coords::Matrix{Float64}, rhs::Matrix{Float64}, n::Int64 )

    Sn = convertmat1(n)*convertmat0(n);
    M2xn = multmat2(n);

    coords_to_params = [ 1  1  1  1;
                        1 -1 -1  1;
                        1  1 -1 -1;
                        1 -1  1 -1] / 4;

    params = coords_to_params * coords;

    # xx, yy = meshgrid(chebpts(n));
    #
    # # x, y on quadrilateral
    # x = params[1, 1] + xx.*params[2, 1] + yy.*params[3, 1] + xx.*yy.*params[4, 1];
    # y = params[1, 2] + xx.*params[2, 2] + yy.*params[3, 2] + xx.*yy.*params[4, 2];


    # Row numbers
    lrownums = n*(n - 2) + 1:n*n - n;
    rrownums = n*(n - 1) + 1:n*n;
    drownums = n*(1:n-2) - 1;
    urownums = n*(1:n-2);

    f = rhs;
#     rhs must be converted to Chebyshev coefficients
    f = Sn * f * Sn.';
#     rhs must also be multiplied by det(r, s)^3
    f = dets(dets(dets(f, params, M2xn), params, M2xn), params, M2xn);
#     High order coefficients must be zeroed out
    f[end-1:end, :] = 0;
    f[:, end-1:end] = 0;

    f
end
