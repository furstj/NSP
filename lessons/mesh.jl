using Plots
using StaticArrays
using LinearAlgebra
import Base.-, Base.+, Base.*, Base./
import Base.getindex

const Scalar = Float64;
const Vector = SVector{2,Float64};
const Label  = Int;


const ScalarList = Array{Scalar,1};
const VectorList = Array{Vector,1};
const LabelList  = Array{Label,1};
const LabelListList = Array{LabelList,1};


×(u::Vector, v::Vector) = u[1]*v[2] - u[2]*v[1];


struct PatchInfo
    name  :: String
    range :: UnitRange{Label}
end;


const PatchInfoList = Array{PatchInfo,1};


struct Mesh
    point      :: VectorList     # Souradnice vrcholu site 
    centre     :: VectorList     # Souradnice stredu bunek 
    volume     :: ScalarList     # Velikosti (objemy) bunek
    surface    :: VectorList     # Normalovy vektor na stenu 
    facecentre :: VectorList     # Souradnice stredu strany
    owner      :: LabelList      # Vlastnik steny
    neighbor   :: LabelList      # Soused steny
    patch      :: PatchInfoList  # Seznam okrajovych podminek
    
    cell2point :: LabelListList  # Indexy vrcholu bunek
    face2point :: LabelListList  # Indexy vrcholu sten
end


"""
    Mesh(point, owner, neighbor, patch, cell2point, face2point)

    Vytvori sit ze zadanych souradnic vrcholu a informaci o konektivite
"""
function Mesh(point::VectorList, owner::LabelList, neighbor::LabelList, patch::PatchInfoList, cell2point::LabelListList, face2point::LabelListList)
    m = Mesh(point, VectorList(), ScalarList(), VectorList(), VectorList(), owner, neighbor, patch, cell2point, face2point)
    update!(m)
    return m
end


"""
    update!(mesh::Mesh)

    Vypocte objemy bunek, souradnice stredu bunek, normalove vektory na steny a stredy sten.
"""
function update!(mesh::Mesh)
    
    # Recalculate cell centres & volumes
    ncells = size(mesh.cell2point,1)
    empty!(mesh.centre)
    empty!(mesh.volume)
    for c=1:ncells
        pts = mesh.point[mesh.cell2point[c]]
        
        v = Scalar(0)
        x = Vector(0,0)
        p1 = pts[end]
        for i=1:size(pts,1)-2
            p2 = pts[i]
            p3 = pts[i+1]
            vv = (p2 - p1) × (p3 - p1) / 2
            v += vv
            x += vv*(p1 + p2 + p3)/3
        end
        push!(mesh.volume, v)
        push!(mesh.centre, x/v)
    end
    
    # Recalculate face sizes & positions
    empty!(mesh.surface)
    empty!(mesh.facecentre)
    nfaces = size(mesh.face2point,1)
    for f=1:nfaces
        pts = mesh.point[mesh.face2point[f]]
        push!(mesh.surface, Vector(pts[2][2]-pts[1][2], pts[1][1]-pts[2][1]))
        push!(mesh.facecentre, (pts[1]+pts[2])/2)
    end 
end

"""
    cartesian_mesh(nx, ny)

    Vytvori kartezskou sit o nx*ny bunkach v jednotkovem ctverci. Okrajove podminky pojmenuje "left", "right", "top" a "bottom"
"""
function cartesian_mesh(nx, ny)
    Δx, Δy = 1.0/nx, 1.0/ny
    
    pid(i,j) = i + j*(nx+1) + 1
    cid(i,j) = i + (j-1)*nx


    # Mesh points
    point = VectorList() 
    for j=0:ny, i=0:nx
        push!(point, Vector(i*Δx, j*Δy))
    end
    
    # Mesh cells
    cell2point = LabelListList()
    for j=1:ny, i=1:nx
        pts = [pid(i-1,j-1), pid(i,j-1), pid(i,j), pid(i-1,j)] |> LabelList
        push!(cell2point, pts)
    end
    
    # Mesh faces
    owner = LabelList()
    neighbor = LabelList()
    face2point = LabelListList()
    
    # Internal faces
    for j=1:ny, i=1:nx-1
        push!(owner, cid(i,j))
        push!(neighbor, cid(i+1,j))
        push!(face2point, [pid(i,j-1), pid(i,j)])
    end

    for j=1:ny-1, i=1:nx
        push!(owner, cid(i,j))
        push!(neighbor, cid(i,j+1))
        push!(face2point, [pid(i,j), pid(i-1,j)])
    end
    
    # Boundary patches
    # - bottom
    patch = PatchInfoList()
    j=0
    start = size(owner,1)+1
    for i=1:nx
        push!(owner, cid(i,j+1))
        push!(face2point, [pid(i-1,j), pid(i,j)])
    end
    push!(patch, PatchInfo("bottom", UnitRange(start, size(owner,1))))

    # - right
    i=nx
    start = size(owner,1)+1
    for j=1:ny
        push!(owner, cid(i,j))
        push!(face2point, [pid(i,j-1), pid(i,j)])
    end
    push!(patch, PatchInfo("right", UnitRange(start, size(owner,1))))
    
    # - top
    j=ny
    start = size(owner,1)+1
    for i=1:nx
        push!(owner, cid(i,j))
        push!(face2point, [pid(i,j), pid(i-1,j)])
    end
    push!(patch, PatchInfo("top", UnitRange(start, size(owner,1))))

    # - left
    i=0
    start = size(owner,1)+1
    for j=1:ny
        push!(owner, cid(i+1,j))
        push!(face2point, [pid(i,j), pid(i,j-1)])
    end
    push!(patch, PatchInfo("left", UnitRange(start, size(owner,1))))
    
    return Mesh(point, owner, neighbor, patch, cell2point, face2point)
end


cells(m::Mesh) = UnitRange(1,length(m.volume))

internal_faces(m::Mesh) = UnitRange(1,length(m.neighbor))

all_faces(m::Mesh) = UnitRange(1,length(m.owner))

boundary_patches(m::Mesh) = UnitRange(1,length(m.patch))

patch_names(m::Mesh) = [ m.patch[p].name for p in boundary_patches(m)]

patch_by_name(m::Mesh, name::String) = findfirst(x->x==name, patch_names(m))

patch_faces(m::Mesh, patch::Label) = m.patch[patch].range


function plot_mesh(m::Mesh)
    plt = plot(aspect_ratio=:equal, legend=:none)

    for f in internal_faces(m)
        pts = m.point[m.face2point[f]]
        x = [pts[1][1], pts[2][1]]
        y = [pts[1][2], pts[2][2]]
        plot!(plt, x, y, c="black")
    end
    
    for p in boundary_patches(m)
        for f in patch_faces(m, p)
            pts = m.point[m.face2point[f]]
            x = [pts[1][1], pts[2][1]]
            y = [pts[1][2], pts[2][2]]
            plot!(plt, x, y, c="black", width=2)
        end
    end
    
    return plt
end
