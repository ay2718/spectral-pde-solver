VER = v"0.7.0"
if VERSION >= VER
    error("Julia version < "*string(VER)*" required!")
end

println("Loading navier_stokes.jl...");
include("navier_stokes.jl");

println("Loading saveimages.jl...");
include("saveimages.jl");

res = 14;
println("Mesh resolution (default is ", res, "): ");
try
    res = parse(Int64, readline(STDIN));
end

len = 3;
println("Aspect ratio (default is ", len, "): ");
try
    len = parse(Int64, readline(STDIN));
end
len -= 1;

bound = 0.007
println("Boundary layer thickness (default is ", bound, "): ");
try
    bound = parse(Float64, readline(STDIN));
end

deg = 10;
println("Polynomial degree (default is ", deg, "): ");
try
    deg = parse(Int64, readline(STDIN));
end

step = 0.000005;
println("Time step (default is ", step, "): ");
try
    step = parse(Float64, readline(STDIN));
end

num_steps = 10000;
println("Number of steps (default is ", num_steps, "): ");
try
    num_steps = parse(Int64, readline(STDIN));
end

fluid_speed = 300.0;
println("Fluid speed (default is ", fluid_speed, "): ");
try
    fluid_speed = parse(Float64, readline(STDIN));
end

println("Generating mesh...");
m = skinny_ns_mesh(res, len);
m = boundary(m, bound);

println("Generating matrices...");
dat = navier_stokes_matrices(m, deg, step);

println("Running fluid flow simulation...");
U = navier_stokes_solve(dat, 10000, 300.0);

while true;
    println("Would you like to create an image? (Y/n)");
    yn = readline(STDIN);

    if yn == "n" || yn == "N"
        break;
    end

    println("Image file path: ");
    path = readline(STDIN);

    st = 1;
    println("Image step: ");
    try
        st = parse(Int64, readline(STDIN));
    end
    st = div(st, 10);
    if st < 1
        st = 1;
    end

    typ = 4;
    println("Enter 1 for x velocity, 2 for y velocity, 3 for pressure, 4 for vorticity (default is ", color_range, "):")
    try
        typ = parse(Int64, readline(STDIN));
    end
    if typ > 4
        typ = 4;
    end
    if typ < 1
        typ = 4;
    end

    color_range = 5000.0;
    println("Color range (Default is ", color_range ," for vorticity):");
    try
        color_range = parse(Float64, readline(STDIN));
    end
    
    pixel_size = 0.002;
    println("Pixel size (Default is ", pixel_size, "):");
    try
        pixel_size = parse(Float64, readline(STDIN));
    end

    println("Creating image...");
    img = createimage(m, U[st, typ], pixel_size, color_range);
    saveimage(path, img);
end
