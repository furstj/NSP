using StaticArrays
import Base.-, Base.+, Base.*, Base./
import Base.getindex

const Vec2d = SVector{2,Float64};
const VecList = Array{Vec2d,1};


struct Cell
    id  :: Int
    x   :: Vec2d
    vol :: Float64
end

const CellList = Array{Cell,1};


struct Face
    id    :: Int
    x     :: Vec2d
    s     :: Vec2d
    owner :: Int
    neigh :: Int
end;

const FaceList = Array{Face,1};


const PatchDict = Dict{String,FaceList};


struct ConnectivityList
    start  :: Array{Int,1}
    target :: Array{Int,1}
end

getindex(c::ConnectivityList, i::Int) = view(c.target, c.start[i]:c.start[i+1]-1)

struct Mesh
    points  :: VecList
    faces   :: FaceList
    cells   :: CellList
    patches :: PatchDict
    cell2points :: ConnectivityList
end


function plot_mesh(m::Mesh, style="-k")
    for c in m.cells
        pids = m.cell2points[c.id]
        x = [ m.points[p][1] for p ∈ pids ]; push!(x, x[1])
        y = [ m.points[p][2] for p ∈ pids ]; push!(y, y[1])        
        plot(x,y, style)
    end
end


function triangles(m::Mesh)
    tri = Array{Array{Int,1},1}()
    for c in m.cells
        pids = m.cell2points[c.id]
        for i=3:length(pids)
            push!(tri, [pids[1], pids[i-1], pids[i]])
        end
    end
    return tri
end


function cartesian_mesh(nx, ny)
    Δx, Δy = 1.0/nx, 1.0/ny
    
    pid(i,j) = i + j*(nx+1) + 1
    cid(i,j) = i + (j-1)*nx

    fid_  = 0 
    function fid() 
        fid_+=1
        return fid_
    end
    
    # mesh points
    points = VecList()
    for j=0:ny,i=0:nx
        push!(points, Vec2d(i*Δx, j*Δy))
    end
    
    faces = FaceList()
    for j=1:ny, i=1:nx
        if i<nx
            x = Vec2d(i*Δx,(j-0.5)*Δy)
            s = Vec2d(Δy,0)
            push!(faces, Face(fid(), x, s, cid(i,j), cid(i+1,j)))
        end
        
        if j<ny
            x = Vec2d((i-0.5)*Δx,j*Δy)
            s = Vec2d(0,Δx)
            push!(faces, Face(fid(), x, s, cid(i,j), cid(i,j+1)))
        end
            
    end
    
    # mesh cells
    cells = CellList()
    for j=1:ny,i=1:nx
        push!(cells, Cell(cid(i,j), Vec2d((i-0.5)*Δx,(j-0.5)*Δy), Δx*Δy))
    end
    
    # boundary patches
    patches = PatchDict(
        "bottom" => [ Face( fid(), Vec2d((i-0.5)*Δx,0), Vec2d(0,-Δy), cid(i,1), 0) for i=1:nx],
        "top"    => [ Face( fid(), Vec2d((i-0.5)*Δx,1), Vec2d(0,Δy), cid(i,ny), 0) for i=1:nx],
        "left"   => [ Face( fid(), Vec2d(0,(j-0.5)*Δy), Vec2d(-Δx,0), cid(1,j), 0) for j=1:ny],
        "right"  => [ Face( fid(), Vec2d(1,(j-0.5)*Δy), Vec2d(Δx,0),  cid(nx,j), 0) for j=1:ny],
    )

    c2p_start = [1+4*(i-1) for i=1:nx*ny+1]
    c2p_target = []
    for j=1:ny, i=1:nx
        append!(c2p_target, [pid(i-1,j-1), pid(i,j-1), pid(i,j), pid(i-1,j)])
    end
    
    return Mesh(points, faces, cells, patches, ConnectivityList(c2p_start, c2p_target) )
end