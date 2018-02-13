# Conversion functions between chebyshev polynomial space and value space at chebyshev points

# Planned 2D Chebyshev conversion (similar to planned 2D DCT-1)
type chplan
    p::Base.DFT.FFTW.r2rFFTWPlan{Float64}
    mat::Matrix{Float64}
    imat::Matrix{Float64}
    dim::Int64
end

# Planned 2D Chebyshev conversion on a dim x dim grid (forward and backward)
function cheb_plan( dim::Int64 )
#     planned 2d DCT-1
    p = FFTW.plan_r2r(Matrix{Float64}(dim, dim), FFTW.REDFT00);
    mat = ones(dim, dim);
    mat[:, 2:end-1]*=2;
    mat[2:end-1, :]*=2;
    imat = 1./mat;
    d = copy(dim);
    dim -= 1;
    dim *= 2;
    dim*=dim;
    mat /= dim;
    pl = chplan(p, mat, imat, d);
    return pl
end

# Uses a planned 2D Chebyshev conversion to go from value space to coefficient space
function vals2coeffs( vals::Matrix{Float64}, plan::chplan )
    coeffs = plan.mat.*(plan.p*vals[end:-1:1, end:-1:1]);
    coeffs
end

# Uses a planned 2D Chebyshev conversion to go from coefficient space to value space
function coeffs2vals( coeffs::Matrix{Float64}, plan::chplan )
    vals = plan.p*(plan.imat.*coeffs);
    vals = vals[end:-1:1, end:-1:1];
    vals
end

# Multiplies matrices of 2D Chebyshev polynomial coefficients by converting to value space, multiplying the values, then converting back to coefficient space
function coeffs_multiply( a::Matrix{Float64}, b::Matrix{Float64}, plan::chplan )
    aa = zeros(plan.dim, plan.dim);
    bb = zeros(plan.dim, plan.dim);
    aa[1:size(a, 1), 1:size(a, 2)] = a;
    bb[1:size(b, 1), 1:size(b, 2)] = b;
    aa = coeffs2vals(aa, plan);
    bb = coeffs2vals(bb, plan);
    cc = aa.*bb;
    cc = vals2coeffs(cc, plan);
    cc
end

# Same functionality as above, but uses Vandermonde matrices
function coeffs_multiply( a::Matrix{Float64}, b::Matrix{Float64}, vand::Matrix{Float64}, ivand::Matrix{Float64})
    aa = zeros(size(vand));
    bb = zeros(size(vand));
    aa[1:size(a, 1), 1:size(a, 2)] = a;
    bb[1:size(b, 1), 1:size(b, 2)] = b;
    aa = vand*aa*vand.';
    bb = vand*bb*vand.';
    cc = aa.*bb;
    cc = ivand*cc*ivand.';
    cc
end

# Vandermonde matrix that converts a vector of Chebyshev coefficients to values at Chebyshev points
function cheb_vander( n::Int64 )
    vand = zeros(n, n);
    vand[1:n, 1] = 1;
    list = linspace(pi, 0, n);
    for i = 1:(n-1)
        vand[1:n, i + 1] = cos(list.*i);
    end
    vand
end

# Inverse Vandermonde matrix that converts a vectors of values at Chebyshev points to Chebyshev coefficients
function inv_cheb_vander( n::Int64 )
    vand = cheb_vander(n);
    maskvec = mod(1:n, 2)*2 - 1;
    maskvec[2:n-1]*=2;
    mask = Diagonal(maskvec);
    vand = mask*vand*mask/(2*n-2);
    vand
end
