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
    nodes = Array{Vector,1}()
    for i=1:nnodes
        words = split(readline(gmsh))
        x = parse(Float64, words[2]) 
        y = parse(Float64, words[3]) 
        push!(nodes, Vector(x,y))
    end

    Edge = Tuple{Int64,Int64}

    
    find_keyword("Elements")
    cell2point = LabelListList()
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

        if etype==1  # line
            push!( boundary_lines[numbers[4]], (numbers[4+ntags],numbers[5+ntags]) )
        elseif etype==2 || etype==3  # triangle or quad
            push!(cell2point, numbers[4+ntags:end])
        else
            error("Unknown element type ", etype)
        end            
    end

    close(gmsh)

    
    
    #### Build data structures for Mesh

    edge_owner = Dict{Edge,Label}()

    owner = LabelList()
    neighbor = LabelList()
    face2point = LabelListList()
    
    push_face(p1, p2, o, n=nothing) = begin
        push!(owner, o)
        if n != nothing; push!(neighbor, n); end
        push!(face2point, [p1, p2])
    end
    
    # Process all faces and push internal faces into owner/neighbor/face2point
    for o = 1:length(cell2point)
        pts = [cell2point[o]; cell2point[o][1]]
        for i = 1:length(pts)-1
            p1, p2 = pts[i], pts[i+1]
            
            if haskey(edge_owner, (p2,p1))
                n = edge_owner[(p2,p1)]
                delete!(edge_owner, (p2,p1))
                push_face(p1, p2, o, n)
            else
                edge_owner[(p1,p2)] = o
            end
        end
    end

    patch = PatchInfoList()
    for (tag, edges) in boundary_lines
        name = physical_lines[tag]
        start = length(owner)+1
        stop  = start + length(edges) - 1
        push!(patch, PatchInfo(name, start:stop))

        for e in edges
            p1, p2 = e[1], e[2]
            if haskey(edge_owner, (p1, p2))
                o = edge_owner[(p1, p2)]
                push_face(p1, p2, o)
            elseif haskey(edge_owner, (p2, p1))
                o = edge_owner[(p2, p1)]
                push_face(p2, p1, o)
            else
                error("Error in processing mesh faces!")
            end
        end
    end

    return Mesh(nodes, owner, neighbor, patch, cell2point, face2point)
end

#mesh = gmsh_mesh("ctverec.msh");
#println(mesh)

