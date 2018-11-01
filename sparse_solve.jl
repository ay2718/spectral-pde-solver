VER = v"0.7.0"
if VERSION >= VER
    error("Julia version < "*string(VER)*" required!")
end

# A \ B, where A is a sparse matrix UMFPACKLU object and B is a sparse matrix (this is useful if B has many zero columns)

import Base.\;

function (\){Tv,Ti}(A::SparseArrays.UMFPACK.UmfpackLU{Tv, Ti}, B::SparseMatrixCSC{Tv,Ti})
  mA, nA = size(A)
  mB, nB = size(B)
  if nA != mB; error("mismatched dimensions"); end
  # @show typeof(rowvalA)
  # @show size(rowvalA)
  # @show typeof(nzvalA)
  # @show size(nzvalA)
  colptrB = B.colptr; rowvalB = B.rowval; nzvalB = B.nzval
  colptrBB = unique(B.colptr);
  nBB = size(colptrBB, 1) - 1;
  BB = SparseMatrixCSC(mB, nBB, colptrBB, rowvalB, nzvalB);
  BB = full(BB);
  BB = A \ BB;
  BB = sparse(BB);
  colptrBB = BB.colptr;
  colptrBBB = zeros(Int64, nB + 1);
  j = 1;
  colptrBBB[1] = colptrBB[j];
  for i = 1:nB
    if colptrB[i] == colptrB[i+1];
      colptrBBB[i+1] = colptrBB[j];
    else
      j+=1;
      colptrBBB[i+1] = colptrBB[j];
    end
  end
  BBB = SparseMatrixCSC(mB, nB, colptrBBB, BB.rowval, BB.nzval);
  BBB
end
