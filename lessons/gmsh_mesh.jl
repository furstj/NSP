function gmsh_mesh(filename)
    
    if !isfile(filename)
        println("Error in gmsh.read():")
        println("  Could not find " * filename * "!")
        exit(1)
    end

    gmsh = open(filename)
    
    line = ""
    
    find_keyword(name) = begin
        while line != "\$"*name
            line = readline(gmsh)
        end
    end
    
    
    find_keyword("MeshFormat")
    mesh_format = readline(gmsh)
    
    
    physical_lines=Dict{Int,String}()
    find_keyword("PhysicalNames")
    nphys = parse(Int, readline(gmsh) )
    for i=1:nphys
        words = split(readline(gmsh))
        d   = parse(Int, words[1])
        num = parse(Int, words[2])
        name = words[3][2:end-1]
        if d==1
            physical_lines[num] = name
        end
    end
    
    
    find_keyword("Nodes")
    nnodes = parse(Int, readline(gmsh))
    nodes = Array{Vec2d,1}()
    for i=1:nnodes
        words = split(readline(gmsh))
        x = parse(Float64, words[2]) 
        y = parse(Float64, words[3]) 
        push!(nodes, Vec2d(x,y))
    end

    const Edge = Tuple{Int64,Int64}
    
    find_keyword("Elements")
    elements = Array{Array{Int,1},1}()
    boundary_lines = Dict{Int, Array{Edge,1}}()
    for (k,v) in physical_lines
        boundary_lines[k] = []
    end
    n = parse(Int, readline(gmsh))
    for i = 1:n
        numbers = [parse(Int,w) for w in split(readline(gmsh))]
        eid   = numbers[1]
        etype = numbers[2]
        ntags = numbers[3]

        if etype==1
            push!( boundary_lines[numbers[4]], (numbers[4+ntags],numbers[5+ntags]) )
        elseif etype==2 || etype==3
            push!(elements, numbers[4+ntags:end])
        else
            error("Unknown element type ", etype)
        end            
    end

    close(gmsh)

    #### Build data structures for Mesh
    face_buffer = Dict{Edge,Int64}()
    face_id = 0

    push_face!(flist, p1, p2, owner, neigh) = begin
        face_id += 1
        a, b = nodes[p1], nodes[p2]
        x = (a + b) / 2.0
        s = Vec2d(b[2]-a[2],a[1]-b[1])
        push!(flist, Face(face_id, x, s, owner, neigh))
    end
    
    # Internal faces
    faces = FaceList()
    for (owner,e) in enumerate(elements)
        cnodes = [e' e[1]]
        for i in 1:length(cnodes)-1
            p1 = cnodes[i]
            p2 = cnodes[i+1]

            if haskey(face_buffer, (p2,p1))
                neigh = face_buffer[(p2,p1)]
                delete!(face_buffer, (p2,p1))
                push_face!(faces, p1, p2, owner, neigh)
            else
                face_buffer[(p1,p2)] = owner
            end       
        end
    end
    
    # Boundary faces
    patches = PatchDict()
    fid = length(faces) + 1
    for (tag,edges) in boundary_lines
        name = physical_lines[tag]
        bfaces = FaceList()
        for e in edges
            if haskey(face_buffer,e)
                p1,p2 = e[1], e[2]
                owner = face_buffer[e]
                push_face!(bfaces, p1, p2, owner, 0)
                delete!(face_buffer, e)
            elseif haskey(face_buffer, (e[2],e[1]))
                p1, p2 = e[2], e[1]
                owner = face_buffer[(p1,p2)]
                push_face!(bfaces, p1, p2, owner, 0)
                delete!(face_buffer, (p1,p2))
            else
                error("Error in processing mesh faces!")
            end
        end
        patches[name] = bfaces
    end

    # Calculating cell volumes
    vol = zeros(length(elements))
    mass = zeros(Vec2d, length(elements))
    for f in faces
        o = f.owner
        n = f.neigh
        v = dot(f.x,f.s) / 2
        vol[o] += v
        vol[n] -= v
        mass[o] += 2.0/3.0*f.x*v
        mass[n] -= 2.0/3.0*f.x*v
    end

    for (n,flist) in patches
        for f in flist
            o = f.owner
            v = dot(f.x,f.s) / 2
            vol[o] += v
            mass[o] += 2.0/3.0*f.x*v
        end
    end

    cells = CellList()
    for cid in 1:length(vol)
        push!(cells, Cell(cid, mass[cid]/vol[cid], vol[cid]) )
    end

    return Mesh(nodes, faces, cells, patches)
end

#mesh = gmsh_mesh("ctverec.msh");
#println(mesh)

