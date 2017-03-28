using StaticArrays

typealias Vec2d SVector{2,Float64};

type CartesianMesh
    nx 
    ny 
    Δx
    Δy
end

CartesianMesh(nx, ny) = CartesianMesh(nx, ny, 1.0/nx, 1.0/ny);

immutable Cell
    id  :: Int
    x   :: Vec2d
    vol :: Float64
end

immutable Face
    x     :: Vec2d
    s     :: Vec2d
    owner :: Int
    neigh :: Int
end;

immutable Patch
    name
    faces
end


function cell(msh::CartesianMesh, id) 
    i,j = rem(id-1, msh.nx)+1, div(id-1, msh.nx)+1
    return Cell(id, Vec2d( (i-1/2)*msh.Δx, (j-1/2)*msh.Δy ), msh.Δx*msh.Δy )
end;


cells(msh::CartesianMesh) = @task begin
    id = 1
    vol = msh.Δx * msh.Δy
    for y in range(msh.Δy/2, msh.Δy, msh.ny)
        for x in range(msh.Δx/2, msh.Δx, msh.nx)
            produce( Cell(id, Vec2d(x,y), vol) )
            id += 1
        end
    end
end;


internalfaces(msh::CartesianMesh) = @task begin
    id(i,j) = i + (j-1)*msh.nx
    Δx, Δy = msh.Δx, msh.Δy
    
    s = Vec2d(Δy, 0.0)
    for j=1:msh.ny
        for i=1:msh.nx-1
            owner, neigh = id(i,j), id(i+1,j)
            x = Vec2d( i*Δx, (j-1/2)*Δy )
            produce( Face(x, s, owner, neigh) )
        end
    end
    
    s = Vec2d(0.0, Δx)
    for j=1:msh.ny-1
        for i=1:msh.nx
            owner, neigh = id(i,j), id(i,j+1)
            x = Vec2d( (i-1/2)*Δx, j*Δy )
            produce( Face(x, s, owner, neigh) )
        end
    end

end;


function patches(msh::CartesianMesh) 
    id(i,j) = i + (j-1)*msh.nx
    Δx, Δy = msh.Δx, msh.Δy
    
    leftfaces() = @task begin
        s = Vec2d(-Δy, 0.0)
        for j=1:msh.ny
            produce( Face(Vec2d(0.0, (j-1/2)*Δy), s, id(1,j), 0) )
        end
    end

    rightfaces() = @task begin
        s = Vec2d(Δy, 0.0)
        for j=1:msh.ny
            produce( Face(Vec2d(msh.nx*Δx, (j-1/2)*Δy), s, id(msh.nx,j), 0) )
        end
    end

    bottomfaces() = @task begin
        s = Vec2d(0.0,-Δx)
        for i=1:msh.nx
            produce( Face(Vec2d((i-1/2)*Δx, 0.0), s, id(i,1), 0) )
        end
    end

    topfaces() = @task begin
        s = Vec2d(0.0,Δx)
        for i=1:msh.nx
            produce( Face(Vec2d((i-1/2)*Δx, msh.ny*Δy), s, id(i,msh.ny), 0) )
        end
    end

    return [ Patch("left", leftfaces), Patch("right", rightfaces), 
        Patch("bottom", bottomfaces), Patch("top", topfaces) ]
end;

name(p::Patch) = p.name;
boundaryfaces(p::Patch) = p.faces();


