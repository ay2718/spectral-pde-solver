VER = v"0.7.0"
if VERSION >= VER
    error("Julia version < "*string(VER)*" required!")
end

# Solves a PDE on an entire mesh

include("onematrix.jl")
include("twoquads.jl")
include("woodbury.jl")
include("sparse_solve.jl")

# PDE matrix container class, includes precomputed matrices for faster computation in fastorder2solve
type PartialMatrix
    m::QuadMesh
    pack::Vector{Woodbury}
    aig::Matrix{Matrix{Float64}}
    agi::Matrix{Matrix{Float64}}
    inv::SparseArrays.UMFPACK.UmfpackLU{Float64,Int64}
    norm::Diagonal{Float64}
    pre::Vector{Diagonal{Float64}}
    n::Int64
end

# There are two variants of this function:  this one works with constant coefficients.  fastorder2mesh_iter1 precomputes everything possible in order to solve a specific PDE on a specific domain
# Since most of the work is already done by fastorder2mesh_iter1, it is much faster to apply multiple right hand sides to the same equation
function fastorder2mesh_iter1( m::QuadMesh, n::Int64, coeffs::Vector{Float64}, neumann_list::BitMatrix = falses(size(m.bcs)), rounded_list::BitMatrix = falses(size(m.bcs)), robin::BitMatrix = falses(size(m.bcs)))
    qnum = size(m.t, 1);
    nbrow = (-1).^(0:n-1);
    pbrow = ones(n);
    pdrow = (0:n-1).^2;
    ndrow = -nbrow.*pdrow;

#     Boundary conditions in the corners
    cornerbcs = Vector{Vector{Float64}}(4);
    cornerbcs[1] = kron(pbrow, pbrow);
    cornerbcs[2] = kron(nbrow, pbrow);
    cornerbcs[3] = kron(nbrow, nbrow);
    cornerbcs[4] = kron(pbrow, nbrow);

    enum = size(m.e, 1);
    ienum = m.ienum;

    A = Vector{SparseMatrixCSC{Float64, Int64}}(qnum);

    Aig = Matrix{Matrix{Float64}}(qnum, 4)
    Agi = Matrix{Matrix{Float64}}(qnum, 4)

    Agg = spzeros(ienum*n, ienum*n);
    pre = Vector{Diagonal{Float64}}(qnum);
    x = Vector{Matrix{Float64}}(qnum);
    y = Vector{Matrix{Float64}}(qnum);

    Agit = [zeros(n, n*n), zeros(n, n*n)];
    Aigt = [zeros(n*n, n), zeros(n*n, n)];

#     Generate matrices for each mesh element
    for i = 1:qnum
        for j = 1:4
            Aig[i, j] = zeros(n*n, n);
            Agi[i, j] = zeros(n, n*n);
        end
        A[i], pre[i], x[i], y[i] = fastorder2matrix( m.p[m.t[i, :], :], n, coeffs);
    end

#     For each edge, set up boundary conditions
    for i = 1:enum
#         If the given edge is on the exterior, we don't need to do anything
        if m.e[i, 5] != 0
            ii = m.e[i, 7];
#             We set up the Dirichlet and Neumann rows
            Aigt, Agit = fastorder2boundary(m, i, n);

#             We take a vertex on the end of the edge that we are using
#             If the vertex is on the boundary, then we add a Dirichlet condition at that vertex
#             If the vertex is on the interior and the edge does not match m.vc of the vertex, then we add a Dirichlet condition at that vertex
#             Here's the rationale--if we have k quadrilaterals meeting at a point, k-1 Dirichlet constraints are needed at that point to ensure continuity.
#             If that point is on the exterior of the mesh, then there are only k-1 edges that meet at that point, so every edge must carry a Dirichlet condition.
#             If the point is on the interior of the mesh, then there are k edges meeting at the point.  We must choose one edge that does not carry a Dirichlet boundary condition for each interior vertex
            
            if m.vc[m.e[i, 1], 1] == 1 || m.vc[m.e[i, 1], 2] != i;
                Agit[1][1, :] = -cornerbcs[m.e[i, 4]];
                Agit[2][1, :] = cornerbcs[mod(m.e[i, 6], 4)+1];
            end

            if m.vc[m.e[i, 2], 1] == 1 || m.vc[m.e[i, 2], 2] != i;
                Agit[1][end, :] = cornerbcs[mod(m.e[i, 4], 4)+1];
                Agit[2][end, :] = -cornerbcs[m.e[i, 6]];
            end
            
#             Now we can put the smaller Aigt and Agit units into the larger Agi and Aig structure
            Aig[m.e[i, 3], m.e[i, 4]] = Aigt[1];
            Aig[m.e[i, 5], m.e[i, 6]] = Aigt[2];

            Agi[m.e[i, 3], m.e[i, 4]] = Agit[1];
            Agi[m.e[i, 5], m.e[i, 6]] = Agit[2];
        end
    end

    packets = Vector{Woodbury}(qnum);
    sigma = Agg;
    len = 0;
    for i = 1:qnum
#         Set up boundary conditions
        BCS = fastorder2ptsbcs(m, i, n, neumann_list[i, :], rounded_list[i, :]);
#         Bonus feature! Limited edition Robin conditions!
        if any(robin[i, :])
            BCS*=0.01
            BCS+=0.99*fastorder2ptsbcs(m, i, n, neumann_list[i, :]|robin[i, :], rounded_list[i, :]);
        end
#         Build Woodbury objects
        packets[i] = woodbury_precondition(A[i], BCS, n);
#         Build sigma matrix
        for j = 1:4
            if !m.bcs[i, j]
                intervalj = n*(m.e[m.te[i, j], 7]-1)+(1:n);
                rhsvec = packets[i]\Aig[i, j];
                for k = 1:4
                    if !m.bcs[i, k]
                        intervalk = n*(m.e[m.te[i, k], 7]-1)+(1:n);
                        sigma[intervalk, intervalj] -= Agi[i, k]*rhsvec;
                    end
                end
            end
        end
#         print output
        str = string("Step ", i, " / ", qnum);
        print("\b"^len, str);
        len = length(str);
    end

#     For added stability, normalize sigma before computing its lu factorization
    normalize_sigma, _ = findmax(abs(sigma), 2);
    normalize_sigma = normalize_sigma[:];
    normalize_sigma = 1./normalize_sigma;
    normalize_sigma = Diagonal(normalize_sigma);
    sigma = normalize_sigma*sigma;
    println("\nInverting Sigma...");
    @time sigmalu = lufact(sigma);

#     If sigma is sufficiently small, display the condition number
    if size(sigma, 1) < 2000
        @show cond(full(sigma))
    end
    data = PartialMatrix(m, packets, Aig, Agi, sigmalu, normalize_sigma, pre, n);
    data
end

# Same function as above, but handles variable coefficients
function fastorder2mesh_iter1( m::QuadMesh, n::Int64, coeffs::Matrix{Matrix{Float64}}, neumann_list::BitMatrix = falses(size(m.bcs)), rounded_list::BitMatrix = falses(size(m.bcs)))
    qnum = size(m.t, 1);
    nbrow = (-1).^(0:n-1);
    pbrow = ones(n);
    pdrow = (0:n-1).^2;
    ndrow = -nbrow.*pdrow;

    cornerbcs = Vector{Vector{Float64}}(4);
    cornerbcs[1] = kron(pbrow, pbrow);
    cornerbcs[2] = kron(nbrow, pbrow);
    cornerbcs[3] = kron(nbrow, nbrow);
    cornerbcs[4] = kron(pbrow, nbrow);

    enum = size(m.e, 1);
    ienum = m.ienum;

    A = Vector{SparseMatrixCSC{Float64, Int64}}(qnum);

    Aig = Matrix{Matrix{Float64}}(qnum, 4)
    Agi = Matrix{Matrix{Float64}}(qnum, 4)

    Agg = spzeros(ienum*n, ienum*n);
    pre = Vector{Diagonal{Float64}}(qnum);
    x = Vector{Matrix{Float64}}(qnum);
    y = Vector{Matrix{Float64}}(qnum);

    Agit = [zeros(n, n*n), zeros(n, n*n)];
    Aigt = [zeros(n*n, n), zeros(n*n, n)];


    for i = 1:qnum
        for j = 1:4
            Aig[i, j] = zeros(n*n, n);
            Agi[i, j] = zeros(n, n*n);
        end
        A[i], pre[i], x[i], y[i] = fastorder2matrix( m.p[m.t[i, :], :], n, coeffs[i, :]);
    end

    for i = 1:enum
        if m.e[i, 5] != 0
            ii = m.e[i, 7];
            Aigt, Agit = fastorder2boundary(m, i, n);

            if m.vc[m.e[i, 1], 1] == 1 || m.vc[m.e[i, 1], 2] != i;
                Agit[1][1, :] = -cornerbcs[m.e[i, 4]];
                Agit[2][1, :] = cornerbcs[mod(m.e[i, 6], 4)+1];
            end

        if m.vc[m.e[i, 2], 1] == 1 || m.vc[m.e[i, 2], 2] != i;
            Agit[1][end, :] = cornerbcs[mod(m.e[i, 4], 4)+1];
            Agit[2][end, :] = -cornerbcs[m.e[i, 6]];
        end
        Aig[m.e[i, 3], m.e[i, 4]] = Aigt[1];
        Aig[m.e[i, 5], m.e[i, 6]] = Aigt[2];

        Agi[m.e[i, 3], m.e[i, 4]] = Agit[1];
        Agi[m.e[i, 5], m.e[i, 6]] = Agit[2];
        end
    end

    packets = Vector{Woodbury}(qnum);
    sigma = Agg;
    len = 0;
    for i = 1:qnum
        BCS = fastorder2ptsbcs(m, i, n, neumann_list[i, :], rounded_list[i, :]);
        packets[i] = woodbury_precondition(A[i], BCS, n);
        for j = 1:4
            if !m.bcs[i, j]
                intervalj = n*(m.e[m.te[i, j], 7]-1)+(1:n);
                rhsvec = packets[i]\Aig[i, j];
                for k = 1:4
                    if !m.bcs[i, k]
                        intervalk = n*(m.e[m.te[i, k], 7]-1)+(1:n);
                        sigma[intervalk, intervalj] -= Agi[i, k]*rhsvec;
                    end
                end
            end
        end
        # sigma -= Agi[i]*(packets[i]\Aig[i]);
        str = string("Step ", i, " / ", qnum);
        print("\b"^len, str);
        len = length(str);
    end

    normalize_sigma, _ = findmax(abs(sigma), 2);
    normalize_sigma = normalize_sigma[:];
    normalize_sigma = 1./normalize_sigma;
    normalize_sigma = Diagonal(normalize_sigma);
    sigma = normalize_sigma*sigma;
    println("\nInverting Sigma...");
    @time sigmalu = lufact(sigma);

    if size(sigma, 1) < 2000
        @show cond(full(sigma))
    end
    data = PartialMatrix(m, packets, Aig, Agi, sigmalu, normalize_sigma, pre, n);
    data
end

# Takes some arbitrary function and returns a Chebyshev series for each element in the mesh (useful for testing)
function function_to_coeffs(m::QuadMesh, n::Int64, func::Function)
    qnum = size(m.t, 1);
    x, y = cheb_mesh_grid(m, n);
    coeffs = Vector{Matrix{Float64}}(qnum);
    chp = cheb_plan(n);
    for i = 1:qnum
        coeffs[i] = Matrix{Float64}(n, n);
        for j = 1:n*n
            coeffs[i][j] = func(x[i][j], y[i][j]);
        end
        coeffs[i] = vals2coeffs(coeffs[i], chp);
    end
    return coeffs
end

function function_to_rhs(m::QuadMesh, n::Int64, func::Function)
    coeffs = function_to_coeffs(m, n, func);
    qnum = size(m.t, 1);
    for i = 1:qnum
        coeffs[i] = fastorder2rhs(m.p[m.t[i, :], :], coeffs[i], n);
    end
    coeffs
end

function function_to_bcs(m::QuadMesh, n::Int64, func::Function)
    enum = size(m.e, 1);
    bcsnum = m.e[end, 7];
    bcs = Matrix{Float64}(bcsnum, n-1);
    bcscounter = 0;
    for i = 1:enum
        if m.e[i, 7] > bcscounter
            bcscounter = m.e[i, 7];
            chpts1 = chebpts(n, [m.p[m.e[i, 1], 1], m.p[m.e[i, 2], 1]])[1:end-1]
            chpts2 = chebpts(n, [m.p[m.e[i, 1], 2], m.p[m.e[i, 2], 2]])[1:end-1]
            bcs[bcscounter, :] = func(chpts1, chpts2)
        end
    end
    bcs
end

# Given the precomputed system, fastorder2solve can efficiently solve the given PDE on the given mesh for some arbitrary right hand side
function fastorder2solve( mesh::PartialMatrix, bcs::Matrix{Float64}, rhs::Vector{Matrix{Float64}})
    qnum = size(mesh.m.t, 1);
    n = mesh.n;
    m = mesh.m;
    nn = n-1;
    ff = Vector{Vector{Float64}}(qnum);

    lrownums = 1:n;
    rrownums = n+1:2*n;
    drownums = 2*n+1:2:4*n-5;
    urownums = 2*n+2:2:4*n-4;

    rownums = Vector{Vector{Int64}}(4);
    rownums[1] = urownums;
    rownums[2] = lrownums;
    rownums[3] = drownums;
    rownums[4] = rrownums;

    brownums = zeros(Int64, n-2, n-2);
    for i = 1:n-2
        for j = 1:n-2
            brownums[i, j] = n*(i-1) + j;
        end
    end
    brownums = brownums[:];

    f = Vector{Float64};

#     When performing multithreaded operations, we can't have all threads writing to the same variable
    intersum = Vector{Vector{Float64}}(Threads.nthreads());
    for i = 1:Threads.nthreads()
        intersum[i] = zeros(size(mesh.inv, 1));
    end

#     We fold boundary conditions into the right hand side
    for i = 1:qnum
        f = rhs[i][1:n-2, 1:n-2][:];
        f = mesh.pre[i]*f;
        bcs_row = zeros(4*n-4);

        for j = 1:4
            if mesh.m.e[mesh.m.te[i, j], 5] == 0;
                bcs_row[nn*(j-1)+(1:nn)] = bcs[mesh.m.e[mesh.m.te[i, j], 7], :];
            end
        end

        ff[i] = zeros(n*n);
        ff[i][mesh.pack[i].bulk] = f;
        ff[i][mesh.pack[i].wood] = bcs_row;
    end

#     Using Schur decomposition, we perform one step of computing the solution: forming the right hand side with which to invert sigma
    Threads.@threads for i = 1:qnum
        rhsvec = (mesh.pack[i]\ff[i])
        for j = 1:4
            if !m.bcs[i, j]
                intervalj = n*(m.e[m.te[i, j], 7]-1)+(1:n);
                intersum[Threads.threadid()][intervalj] -= mesh.agi[i, j]*rhsvec;
            end
        end
    end

    inter = sum(intersum);

#     inverting sigma...
    ug = mesh.inv\(mesh.norm*inter);
    u = Vector{Matrix{Float64}}(qnum);

#     Computing and reshaping the solution
    Threads.@threads for i = 1:qnum
        temp = ff[i]
        for j = 1:4
            if !m.bcs[i, j]
                intervalj = n*(m.e[m.te[i, j], 7]-1)+(1:n);
                temp -= mesh.aig[i, j]*ug[intervalj];
            end
        end
        temp = mesh.pack[i]\temp;
        u[i] = reshape(temp, n, n);
    end
    u
end
