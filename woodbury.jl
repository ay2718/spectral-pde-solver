VER = v"0.7.0"
if VERSION >= VER
    error("Julia version < "*string(VER)*" required!")
end

# We use the Woodbury Matrix Equation to use banded matrices (our matrices are almost banded with boundary condition rows)

import Base.\
include("sparse_solve.jl")

# Contains the LU factorization of the modified matrix and some other required precomputed matrices
type Woodbury
    a::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
    aa::SparseMatrixCSC{Float64, Int64}
    bulk::Vector{Int64}
    wood::Vector{Int64}
    n::Int64
    v::Matrix{Float64}
    corr::Matrix{Float64}
end

# Takes a banded matrix, boundary condition rows, and uses the Woodbury matrix equation to solve a banded matrix, not a nearly banded matrix
function woodbury_precondition(A::SparseMatrixCSC{Float64, Int64}, BCS::Matrix{Float64}, n::Int64)
#     The original matrix A must have its rows strategically placed to optimize the bandwidth
    bulk_rows = zeros(Int64, (n-2)^2);
#     The replacement boundary conditions correspond to constraining the lowest order coefficients of the solution
    woodbury_rows = zeros(Int64, 4*n-4);
    b_ind = 1;
    w_ind = 1;
    for i = 1:n
        for j = 1:n
            if (i > 2)&&(j > 2)
                bulk_rows[b_ind] = n*(i-1)+j;
                b_ind+=1;
            else
                woodbury_rows[w_ind] = n*(i-1)+j;
                w_ind+=1;
            end
        end
    end

#     Woodbury magic occurs here:
    UU = zeros(n^2, 4*n - 4);
    VV = UU.';

    UU[woodbury_rows, :] = eye(4*n - 4);
    VV[:, woodbury_rows] = eye(4*n - 4);
    VV = BCS - VV;

    AA = speye(n*n);
    AA[bulk_rows, :] = A;

    AA = lufact(AA);

    A_U = AA\UU;

    corr_factor = A_U/(eye(4*n-4) + VV*A_U);

    packet = Woodbury(AA, A, bulk_rows, woodbury_rows, n, VV, corr_factor);
    packet
end

# These functions allow the backslash operator to work on a Woodbury object
function (\)(A::Woodbury, B::Array{Float64})
    A_B = A.a\B;
    return (A_B - A.corr*(A.v*A_B))
end

function (\)(A::Woodbury, B::SparseMatrixCSC{Float64, Int64})
    A_B = A.a\B;
    return (A_B - convert(SparseMatrixCSC{Float64, Int64}, A.corr*(A.v*A_B)))
end
