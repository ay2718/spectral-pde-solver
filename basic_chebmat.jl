VER = v"0.7.0"
if VERSION >= VER
    error("Julia version < "*string(VER)*" required!")
end

# Basic conversion matrices for ultraspherical and Chebyshev coefficients

# S0, converts a vector of Chebyshev coefficients to ultraspherical coefficients of parameter 1
function convertmat0(n::Int64)
    S0n = spzeros(n, n);
    for i = 1:n
        S0n[i, i] = 0.5;
    end
    for i = 1:n-2
        S0n[i, i+2] = -0.5;
    end
    S0n[1, 1] = 1;
    S0n
end

#S0^-1, converts a vector of ultraspherical coefficients of parameter 1 to Chebyshev coefficients
function invconvertmat0(n::Int64)
    iS = zeros(n, n);
    for i = 1:n
        if mod(i, 2) == 1
            iS[1:2:i, i] = 2;
            iS[1, i] = 1;
        else
            iS[2:2:i, i] = 2;
        end
    end
    iS
end

# S1, converts a vector of ultraspherical coefficients of parameter 1 to ultraspherical coefficients of parameter 2
function convertmat1(n::Int64)
    S1n = spzeros(n, n);
    for i = 1:n
        S1n[i, i] = 1./i;
    end
    for i = 1:n-2
        S1n[i, i+2] = -1./(i+2);
    end
    S1n
end

# D1, converts a vector of Chebyshev coefficients to its derivative in ultraspherical coefficients of parameter 1
function diffmat1(n::Int64)
    D1n = spzeros(n, n);
    for i = 1:n-1
        D1n[i, i+1] = i;
    end
    D1n
end

# D2, converts a vector of Chebyshev coefficients to its second derivative in ultraspherical coefficients of parameter 2
function diffmat2(n::Int64)
    D2n = spzeros(n, n);
    for i = 1:n-2
        D2n[i, i+2] = 2*i+2;
    end
    D2n
end

# multiplies a vector of ultraspherical coefficients of parameter 1 by x
function multmat1(n::Int64)
    M1xn = spzeros(n, n);
    for i = 1:n-1
        M1xn[i, i+1] = 0.5;
        M1xn[i+1, i] = 0.5;
    end
    M1xn
end

# multiplies a vector of ultraspherical coefficients of parameter 2 by x
function multmat2(n::Int64)
    M2xn = spzeros(n, n);
    for i = 1:n-1
        M2xn[i, i+1] = 0.5*(3+i)./(2+i);
        M2xn[i+1, i] = 0.5*i./(1+i);
    end
    M2xn
end
