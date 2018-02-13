# All sorts of fun and exciting mesh processing methods

# Contains everything you need to know about a quadrilateral mesh
type QuadMesh
#     t is a list of quadrilaterals (qnum x 4 matrix), where each row contains the vertex numbers contained by each element (in counterclockwise order)
    t::Matrix{Int64}
#     p is a list of coordinates (pnum x 2 matrix), where each row contains the x and y coordinates of a mesh node
    p::Matrix{Float64}
#     e is kind of complicated.  Each edge has its own row.  Columns 1 and 2 are the vertices on the edge's endpoints.
#     Column 3 is the quadrilateral number of a (or the) quadrilateral adjacent to the edge, and Column 4 is the local vertex number contained by that edge.
#     Columns 5 and 6 contain the second quadrilateral adjacent to the edge, and they are zeroed out of the edge is on the boundary
#     For instance, if edge 9 contained local vertices 2 and 3 of quadrilateral 6, then row 9 of e would look like [#, #, 6, 2, #, #, #] or [#, #, #, #, 6, 2, #]
#     Edges are either boundary edges (one quadrilateral neighbor) or interior edges (two quadrilateral neighbors)
#     Column 7 gives the interior edge number or boundary edge
    e::Matrix{Int64}
#     Each row of vc corresponds to a vertex.  Column 1 is zero iff the vertex is on th exterior of the mesh.  Column 2 gives the number of ONE edge that contains that vertex
    vc::Matrix{Int64}
#     The element in row i and column j is true iff local edge j of quadrilateral i is on the boundary (local edge j contains local vertices j and j+1 mod 4)
    bcs::BitMatrix
#     Row i and column j of te contains the global edge number of local edge j of quadrilateral i
    te::Matrix{Int64}
#     ienum is the number of interior edges in the mesh
    ienum::Int64
end

# generates a QuadMesh structure from a quadrilateral list and a vertex list
function QuadMesh(t::Matrix{Int64}, p::Matrix{Float64})
    n = size(t, 1);
    e = zeros(Int64, 4*n, 7);
#     local edge 1 of quadrilateral q has endpoints of local vertices 1 and 2 of quadrilateral q
    for i = 1:n
        for j = 1:4
            edge_coords = mod(j-1:j, 4) + 1;
            index = (i-1)*4+j;
#             for this step, the endpoints of each edge are sorted for convenience
            e[index, 1:2] = sort(t[i, edge_coords]);
            e[index, 3:4] = [i, j];
        end
    end
#     Here, we need to eliminate duplicate edges
#     First, we sort the edges...
    e_sorted = sortrows(e);
#     Make a copy...
    e_null = copy(e_sorted);
#     Zero out the quadrilateral information in the copy...
    e_null[:, 3:4] = 0;
#     ... and find the unique edges in the copy
    e = unique(e_null, 1);
    edge_size = size(e, 1);
#     index in the unique edge list
    un_index = 1;
    for i = 1:n
        for j = 1:4
#             index in the non-unique list
            ic_index = j + 4*(i-1);
#             If un_index is not out of bounds and element un_index matches element ic_index, then we enter the quadrilateral data into columns 3 and 4 of e and advance un_index (we are looking for the next unique edge)
#             Otherwise, we enter the quadrilateral data into columns 5 and 6 of e in row un_index-1
            if un_index <= edge_size && all(e_sorted[ic_index, 1:2] == e[un_index, 1:2])
                e[un_index, 3:4] = e_sorted[ic_index, 3:4];
                un_index+=1;
            else
                e[un_index-1, 5:6] = e_sorted[ic_index, 3:4];
            end
        end
    end
    
#     number of inside edges
    ins_edge_size = 0;
#     number of boundary edges
    bound_edge_size = 0;
    for i = 1:edge_size
#         If e[i, 5] or e[i, 6] is nonzero, then edge i has two neighboring quadrilaterals, so it is on the interior of the mesh
        if e[i, 5] != 0
            ins_edge_size+=1;
            e[i, 7] = ins_edge_size;
        else
            bound_edge_size+=1;
            e[i, 7] = bound_edge_size;
        end
    end
    
    # Each row of vc is a vertex.  The first column is nonzero iff the vertex is on the exterior.  The second column gives one edge that contains the vertex in that row (this is useful for boundary conditions)
    vc = zeros(Int64, size(p, 1), 2);
    
    for i = 1:edge_size
        if e[i, 5] == 0
            vc[e[i, 1], 1] = 1;
            vc[e[i, 2], 1] = 1;
        end
        vc[e[i, 1], 2] = i;
        vc[e[i, 2], 2] = i;
    end
    
#     If bcs_list[i, j] = 0, then local edge j of quadrilateral i is on the interior
    bcs_list = falses(n, 4);
#     quad_list[i, j] is the global edge number of local edge j of quadrilateral i
    quad_list = zeros(Int64, n, 4);
    for i = 1:edge_size
        e[i, 1] = t[e[i, 3], e[i, 4]]
        e[i, 2] = t[e[i, 3], mod(e[i, 4], 4) + 1]
        if e[i, 5] == 0
            bcs_list[e[i, 3], e[i, 4]] = true;
        else
            quad_list[e[i, 5], e[i, 6]] = i;
        end
        quad_list[e[i, 3], e[i, 4]] = i;
    end
    return QuadMesh(t, p, e, vc, bcs_list, quad_list, ins_edge_size);
end

# If two edges that share a quadrilateral have numbers differing by k, then this creates a nonzero region about k*p away from the diagonal in sigma, a matrix used in Schur decomposition (p is the polynomial degree used by the spectral element method)
# In a 2d mesh, we can get about O(N^0.5) differnece between edge numbers, where N is the total number of edges.
# Mesh_order essentially uses a breadth-first search to get an edge ordering (depth-first search would be less reliable at minimizing the maximum difference between adjacent edges)
function mesh_order( m::QuadMesh, seed::Int64 )
    curr_edges = zeros(Int64, size(m.e, 1));
    next_edges = zeros(Int64, size(m.e, 1));
    edge_rank = zeros(Int64, size(m.e, 1));
    edge_nums = zeros(Int64, size(m.e, 1));
    curr_edges[1:size(seed, 1)] = seed;
    
    csz = 1;
    nsz = 1;
    nums = 1;
    rank = 2;
    
    edge_nums[1] = seed;
    
    edge_rank[seed] = 1;
    
    while nsz > 0
        nsz = 0;
        for i = 1:csz
            ci = curr_edges[i];
            others = mod(m.e[ci, 4]:m.e[ci, 4]+2, 4)+1;
            for j = 1:3
                edgeval = m.te[m.e[ci, 3], others[j]];
                if edge_rank[edgeval] == 0
                    edge_rank[edgeval] = rank;
                    nsz += 1;
                    nums += 1;
                    edge_nums[nums] = edgeval;
                    next_edges[nsz] = edgeval;
                end
            end
            if m.e[ci, 5] != 0
                others = mod(m.e[ci, 6]:m.e[ci, 6]+2, 4)+1;
                for j = 1:3
                    edgeval = m.te[m.e[ci, 5], others[j]];
                    if edge_rank[edgeval] == 0
                        edge_rank[edgeval] = rank;
                        nsz += 1;
                        nums += 1;
            
                        edge_nums[nums] = edgeval;
                        next_edges[nsz] = edgeval;
                    end
                end
            end
        end

        csz = nsz;
        curr_edges[1:csz] = next_edges[1:csz];
        rank += 1;

    end

    m.e = m.e[edge_nums, :];
    m.vc = zeros(Int64, size(m.p, 1), 2);
    for i = 1:size(m.e, 1)
        if m.e[i, 5] == 0
            m.vc[m.e[i, 1], 1] = 1;
            m.vc[m.e[i, 2], 1] = 1;
        end
        m.vc[m.e[i, 1], 2] = i;
        m.vc[m.e[i, 2], 2] = i;
    end

    m.bcs = falses(size(m.t, 1), 4);
    m.te = zeros(Int64, size(m.t, 1), 4);
    for i = 1:size(m.e, 1)
        m.e[i, 1] = m.t[m.e[i, 3], m.e[i, 4]]
        m.e[i, 2] = m.t[m.e[i, 3], mod(m.e[i, 4], 4) + 1]
        if m.e[i, 5] == 0
            m.bcs[m.e[i, 3], m.e[i, 4]] = true;
        else
            m.te[m.e[i, 5], m.e[i, 6]] = i;
        end
        m.te[m.e[i, 3], m.e[i, 4]] = i;
    end

    ins_edge_size = 0;
    bound_edge_size = 0;
    for i = 1:size(m.e, 1);
        if m.e[i, 5] != 0
            ins_edge_size+=1;
            m.e[i, 7] = ins_edge_size;
        else
            bound_edge_size+=1;
            m.e[i, 7] = bound_edge_size;
        end
    end
    m
end

# Parametric equation of a tilted ellipse
function ellipse(major::Float64, minor::Float64, angle::Float64, t::Float64)
    coords = zeros(2);
    coords[1] = major*cos(t-angle)*cos(angle) - minor*sin(t-angle)*sin(angle);
    coords[2] = minor*sin(t-angle)*cos(angle) + major*cos(t-angle)*sin(angle);
    return coords
end

# Node spacing for a hole in a mesh
function exp_pts(outer::Vector{Float64}, inner::Vector{Float64}, res::Int64) # bar::Float64
    ratio = sqrt((outer[1]^2 + outer[2]^2)/(inner[1]^2 + inner[2]^2));
    pts = exp((-res:0)*log(ratio)/(res-1));
    pts -= pts[1];
    pts /= pts[end];
    pts2d = zeros(res+1, 2)
    pts2d[:, 1] = outer[1]*pts + inner[1]*(1-pts);
    pts2d[:, 2] = outer[2]*pts + inner[2]*(1-pts);
    return pts2d
end

# Constructs a long rectangular quadrilateral mesh with a tilted elliptical hole (Navier-Stokes equations test mesh)
function skinny_ns_mesh(res::Int64, length::Int64 = 2)
    ellipse_res = div(4*res, 5)+1;
    pptr = 4*(ellipse_res+1)*res;

    p = zeros(Float64, 4*(ellipse_res+1)*res + length*res*(res+1), 2);
    t = zeros(Int64, 4*ellipse_res*res + length*res*res, 4);

    major = 0.4;
    minor = 0.25;
    angle = 1.2;

    outerboxpts = zeros(Float64, 4*res, 2);
    innerboxpts = zeros(Float64, 4*res, 2);

    outerboxpts[1:res, 1] = 1.0;
    outerboxpts[1:res, 2] = -1+2*(0:res-1)/res;

    outerboxpts[res+1:2*res, 1] = 1-2*(0:res-1)/res;
    outerboxpts[res+1:2*res, 2] = 1.0;

    outerboxpts[2*res+1:3*res, 1] = -1.0;
    outerboxpts[2*res+1:3*res, 2] = 1-2*(0:res-1)/res;

    outerboxpts[3*res+1:4*res, 1] = -1+2*(0:res-1)/res;
    outerboxpts[3*res+1:4*res, 2] = -1.0;

    for i = 1:4*res
        innerboxpts[i, :] = ellipse(major, minor, angle, (i-1)*pi/(2*res)+7*pi/4);
    end

    for i = 1:4*res
        p[1+(i-1)*(ellipse_res+1):(ellipse_res+1)*i, :] = exp_pts(outerboxpts[i, :], innerboxpts[i, :], ellipse_res);
    end

    for i = 1:res+1
        p[pptr+(1:length*res)+length*res*(i-1), 1] = 1+2*(1:length*res)/res;
        p[pptr+(1:length*res)+length*res*(i-1), 2] = -1+2*(i-1)/res;
    end

    rings = zeros(Int64, 4*res+1, 2)
    k = 1;
    # bias = false;
    for i = 1:ellipse_res
        # if i == 1;
        #   bias = true;
        # else
        #   bias = false;
        # end
        rings[:, 1] = vcat((ellipse_res+1)*(0:4*res-1)+i, i);
        rings[:, 2] = vcat((ellipse_res+1)*(0:4*res-1)+i+1, i+1);
        for j = 1:4*res
            t[k, :] = vcat(rings[j+1:-1:j, 1], rings[j:j+1, 2]);
            k += 1;
        end
    end

    rings = zeros(Int64, length*res+1, 2)
    for i = 1:res
        rings[:, 1] = vcat((ellipse_res+1)*i, pptr+length*res*(i-1)+(1:length*res));
        rings[:, 2] = vcat((ellipse_res+1)*(i+1), pptr+length*res*i+(1:length*res));
        for j = 1:length*res
            t[k, :] = vcat(rings[j:j+1, 1], rings[j+1:-1:j, 2]);
            k += 1;
        end
    end
    m = QuadMesh(t, p);
    m = mesh_order(m, 1);
    return m
end

# next_edge takes a boundary edge and returns a boundary edge.  If we orient the mesh so that the original edge is horizontal and the interior of the mesh is down, then next_edge returns the boundary edge to the right of the original edge
function next_edge(m::QuadMesh, edge::Int64)
    if m.e[edge, 5] != 0
        return edge
    end
    ii = m.e[edge, 3:4];
    edge = m.te[ii[1], mod(ii[2]-2, 4)+1];
    while m.e[edge, 5] != 0
        if m.e[edge, 3] == ii[1]
            ii = m.e[edge, 5:6]
        else
            ii = m.e[edge, 3:4]
        end
        edge = m.te[ii[1], mod(ii[2]-2, 4)+1];
    end
    return edge
end

# Same as next_edge, but returns the edge to the left of the original edge
function prev_edge(m::QuadMesh, edge::Int64)
    if m.e[edge, 5] != 0
        return edge
    end
    ii = m.e[edge, 3:4];
    edge = m.te[ii[1], mod(ii[2], 4)+1];
    while m.e[edge, 5] != 0
        if m.e[edge, 3] == ii[1]
            ii = m.e[edge, 5:6]
        else
            ii = m.e[edge, 3:4]
        end
        edge = m.te[ii[1], mod(ii[2], 4)+1];
    end
    return edge
end

# Constructs a boundary of constant width around a mesh, keeping the original mesh boundary intact
function boundary(m::QuadMesh, width::Float64)
    count = 0;
    for i = 1:length(m.bcs)
        if m.bcs[i]
            count+=1;
        end
    end
    t = vcat(m.t, zeros(Int64, count, 4));
    p = vcat(m.p, zeros(count, 2));
    qnum = size(m.t, 1);
    pnum = size(m.p, 1);
    bcs = copy(m.bcs);

    count = 0;
    for i = 1:size(m.bcs, 1)
        for j = 1:4
            if bcs[i, j]
                index = copy(count);
                edge = m.te[i, j];
                ring = [0, 0, 0, 0];
                ii = [0, 0];
                jj = [0, 0];
                nextedge = next_edge(m, edge);
                nextedgestored = nextedge
                while nextedge != nextedgestored || count == index
                    count += 1
                    ii = m.e[edge, 3:4];
                    jj = m.e[nextedge, 3:4];
                    bcs[ii[1], ii[2]] = false;
                    ring = mod(jj[2]-1:jj[2]+2, 4)+1;
                    ec1 = m.p[m.t[ii[1], mod(ii[2], 4)+1], :] - m.p[m.t[ii[1], ii[2]], :];
                    ec2 = m.p[m.t[jj[1], mod(jj[2], 4)+1], :] - m.p[m.t[jj[1], jj[2]], :];
                    ec1 /= sqrt(ec1[1]^2 + ec1[2]^2);
                    ec2 /= sqrt(ec2[1]^2 + ec2[2]^2);
                    ec3 = (ec1 + ec2)/2;
                    ec3 /= (ec3[1]^2 + ec3[2]^2);
                    ec3 = [-ec3[2], ec3[1]];
                    t[qnum+count, :] = [pnum+count+1, pnum+count, m.t[jj[1], ring[2]], m.t[jj[1], ring[1]]];
                    p[m.t[ii[1], ii[2]], :] += ec3*width;
                    p[pnum+count, :] = m.p[m.t[ii[1], ii[2]], :];
                    edge = nextedge;
                    nextedge = next_edge(m, edge);
                end
                t[qnum + count, 1:2] = [pnum+index+1, pnum+count];
            end
        end
    end
    m = QuadMesh(t, p);
    m = mesh_order(m, 1);
    m
end

function check_mesh(m::QuadMesh)
    qnum = size(m.t, 1)
    
    for i = 1:qnum
        if area(m.p[m.t[i, :], :]) <= 0
            println("Reversed Element!  ", i);
        end
        for j = 1:4
            n = m.te[i, j];
            if m.e[n, 3] == i
                if m.e[n, 4] != j
                    println("Wrong numbering!  ", i, " ", j)
                end
            elseif m.e[n, 5] == i
                if m.e[n, 6] != j
                    println("Wrong numbering!  ", i, " ", j)
                end
            else
                println("Quad not connected to edge!")
            end
            if j == 1 || j == 3
                if !all(m.e[n, [2, 1]] == m.t[i, mod(j-1:j, 4)+1]);
                    println(m.e[n, [2, 1]], " ", m.t[i, mod(j-1:j, 4)+1]);
                end
            else
                if !all(m.e[n, 1:2] == m.t[i, mod(j-1:j, 4)+1]);
                    println(m.e[n, 1:2], " ", m.t[i, mod(j-1:j, 4)+1]);
                end
            end
        end
    end
end

if Pkg.installed("PyPlot") != nothing
    using PyPlot
    function disp_mesh( m::QuadMesh, vals = 0 )
        minax, _ = findmin(m.p, 1);
        maxax, _ = findmax(m.p, 1);
        axis([minax[1], maxax[1], minax[2], maxax[2]]);
        axis("equal");
        plt[:ion]()
        for i = 1:size(m.e, 1)
            if m.e[i, 3] == vals || m.e[i, 5] == vals;
                plot(m.p[m.e[i, 1:2], 1], m.p[m.e[i, 1:2], 2], color = "b", linewidth = 2, zorder = 10)
            else
            plot(m.p[m.e[i, 1:2], 1], m.p[m.e[i, 1:2], 2], color = "k")
        end
        plt[:ioff]()
    end

#     for i = 1:size(m.t, 1)
#       plot([mean(m.p[m.t[i, 1:2], 1]), mean(m.p[m.t[i, 1:4], 1])], [mean(m.p[m.t[i, 1:2], 2]), mean(m.p[m.t[i, 1:4], 2])], color = "r", linewidth = 2, zorder = 11);
#     end
    end
end
