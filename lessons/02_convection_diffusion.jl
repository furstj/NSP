
using Plots
using StaticArrays
import Base.-, Base.+, Base.*, Base./
using ForwardDiff
using LinearAlgebra

pyplot();

const Scalar = Float64;
const Vector = SVector{2,Float64};
const Label  = Int;

×(u::Vector, v::Vector) = u[1]*v[2] - u[2]*v[1];

const ScalarList = Array{Scalar,1};
const VectorList = Array{Vector,1};
const LabelList  = Array{Label,1};
const LabelListList = Array{LabelList,1};

struct PatchInfo
    name  :: String
    range :: UnitRange{Label}
end;

const PatchInfoList = Array{PatchInfo,1};

struct Mesh
    point      :: VectorList
    centre     :: VectorList
    volume     :: ScalarList
    surface    :: VectorList
    facecentre :: VectorList
    owner      :: LabelList
    neighbor   :: LabelList
    patch      :: PatchInfoList
    
    cell2point :: LabelListList
    face2point :: LabelListList
end


Mesh() = Mesh(VectorList(), VectorList(), ScalarList(), VectorList(), VectorList(), LabelList(), LabelList(), PatchInfoList(), LabelListList(), LabelListList())

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

function cartesian_mesh(nx, ny)
    Δx, Δy = 1.0/nx, 1.0/ny
    
    pid(i,j) = i + j*(nx+1) + 1
    cid(i,j) = i + (j-1)*nx

    mesh = Mesh() 

    # Mesh points
    for j=0:ny, i=0:nx
        push!(mesh.point, Vector(i*Δx, j*Δy))
    end
    
    # Mesh cells
    for j=1:ny, i=1:nx
        pts = [pid(i-1,j-1), pid(i,j-1), pid(i,j), pid(i-1,j)] |> LabelList
        push!(mesh.cell2point, pts)
    end
    
    # Internal faces
    for j=1:ny, i=1:nx-1
        push!(mesh.owner, cid(i,j))
        push!(mesh.neighbor, cid(i+1,j))
        push!(mesh.face2point, [pid(i,j-1), pid(i,j)])
    end

    for j=1:ny-1, i=1:nx
        push!(mesh.owner, cid(i,j))
        push!(mesh.neighbor, cid(i,j+1))
        push!(mesh.face2point, [pid(i,j), pid(i-1,j)])
    end
    
    # Boundary patches
    # - bottom
    j=0
    start = size(mesh.owner,1)+1
    for i=1:nx
        push!(mesh.owner, cid(i,j+1))
        push!(mesh.face2point, [pid(i-1,j), pid(i,j)])
    end
    push!(mesh.patch, PatchInfo("bottom", UnitRange(start, size(mesh.owner,1))))

    # - right
    i=nx
    start = size(mesh.owner,1)+1
    for j=1:ny
        push!(mesh.owner, cid(i,j))
        push!(mesh.face2point, [pid(i,j-1), pid(i,j)])
    end
    push!(mesh.patch, PatchInfo("right", UnitRange(start, size(mesh.owner,1))))
    
    # - top
    j=ny
    start = size(mesh.owner,1)+1
    for i=1:nx
        push!(mesh.owner, cid(i,j))
        push!(mesh.face2point, [pid(i,j), pid(i-1,j)])
    end
    push!(mesh.patch, PatchInfo("top", UnitRange(start, size(mesh.owner,1))))

    # - left
    i=0
    start = size(mesh.owner,1)+1
    for j=1:ny
        push!(mesh.owner, cid(i+1,j))
        push!(mesh.face2point, [pid(i,j), pid(i,j-1)])
    end
    push!(mesh.patch, PatchInfo("left", UnitRange(start, size(mesh.owner,1))))

    update!(mesh)
    
    return mesh
end

internal_faces(m::Mesh) = UnitRange(1,size(m.neighbor,1))

boundary_patches(m::Mesh) = UnitRange(1,size(m.patch,1))

patch_names(m::Mesh) = [ m.patch[p].name for p in boundary_patches(m)]

patch_by_name(m::Mesh, name::String) = findfirst(x->x==name, patch_names(m))

patch_faces(m::Mesh, patch::Label) = m.patch[patch].range

m2 = cartesian_mesh(2,2);

for p in m2.point
    println("p = ", p)
end

for c in m2.centre
    println("c = ", c)
end

print("volume[:] = ", m2.volume)

for f in internal_faces(m2)
    println("S=", m2.surface[f], "\t o=", m2.owner[f], "\t n=",m2.neighbor[f])
end

for p in boundary_patches(m2)
    println("PATCH ", m2.patch[p].name)
    for f in patch_faces(m2, p)
        println("S=", m2.surface[f], "\t o=", m2.owner[f])
    end
    println()
end

mutable struct Field
    values
    mesh
    boundaries
end

Field(m::Mesh) = Field( zeros(length(m.centre)), m, Dict{String,Any}());
Field(values, m::Mesh) = Field(values, m, Dict{String,Any}());

for op ∈ [:+, :-, :*, :/]
    @eval begin
        function ($op)(a::Field, b::Field)
            @assert a.mesh == b.mesh
            Field(broadcast($op,a.values,b.values), a.mesh)
        end
    end
end

a2 = Field(m2); a2.values = ones(4);
b2 = Field(m2); b2.values = 3*ones(4);

c2 = a2 * b2;
c2.values

for op ∈ [:+, :-, :*, :/]
    @eval begin
        function ($op)(a::Field, b::Number)
            Field(broadcast($op,a.values,b), a.mesh)
        end
    end
end      

for op ∈ [:+, :-, :*, :/]
    @eval begin
        function ($op)(a::Field, b::Array{Float64,1})
            @assert size(a.values) == size(b) 
            Field(broadcast($op,a.values,b), a.mesh)
        end
    end
end      

-(a::Field) = Field(-a.values, a.mesh);

c2 = a2 / 3;
c2.values

function ←(a::Field, b::Field)
    @assert a.mesh == b.mesh
    a.values = copy(b.values)
end

b2.values = 3*ones(4)
c2 ← b2;
d2 = b2;
b2.values[1] = 5;
println(c2.values)
println(d2.values)

struct BoundaryCondition
    c1    :: ScalarList
    c2    :: ScalarList
    field :: Field
    patch :: Label
end

function boundary_coeffs(bc::BoundaryCondition, f::Label)
    field = bc.field
    mesh = field.mesh
    start = first(mesh.patch[bc.patch].range)
    i = f + 1 - start
    
    c1 = bc.c1[i]
    c2 = bc.c2[i]

    return c1, c2
end

function boundary_value(bc::BoundaryCondition, f::Label)
    field = bc.field
    mesh = field.mesh
    u1 = field.values[mesh.owner[f]]
    c1, c2 = boundary_coeffs(bc, f)
    return c1 + c2*u1
end

function set_dirichlet_patch!(field::Field, name::String, val::Scalar)
    mesh = field.mesh
    p = patch_by_name(mesh, name)
    c1 = [val for f in patch_faces(mesh, p)]
    c2 = zero(c1)
    field.boundaries[name] = BoundaryCondition(c1, c2, field, p)
end;

function set_neumann_patch!(field::Field, name::String)
    mesh = field.mesh
    p = patch_by_name(mesh, name)
    c2 = [1.0 for f in patch_faces(mesh, p)]
    c1 = zero(c2)
    field.boundaries[name] = BoundaryCondition(c1, c2, field, p)
end;

m3 = cartesian_mesh(3,3);

u3 = Field(m3);
set_dirichlet_patch!(u3, "left", 1.0);
set_dirichlet_patch!(u3, "bottom", 0.0);
set_neumann_patch!(u3, "right");
set_neumann_patch!(u3, "top");

v = Vector(1, 0.5);

# Tok typu upwind
H(ul, ur, s) = dot(v,s)>0 ? dot(v,s)*ul : dot(v,s)*ur

function R(u::Field)
    mesh = u.mesh
    r = zero(u.values)    
    
    for f in internal_faces(mesh)
        o, n = mesh.owner[f], mesh.neighbor[f]
        S = mesh.surface[f]
        flux = H(u.values[o], u.values[n], S)
        r[o] += flux
        r[n] -= flux
    end
    
    for p in boundary_patches(mesh)
        name = mesh.patch[p].name
        bc = u.boundaries[name]
        
        for f in patch_faces(mesh, p)
            o = mesh.owner[f]
            S = mesh.surface[f]
            ub = boundary_value(bc, f)
            flux = H(u.values[o], ub, S)
            r[o] += flux
        end
    end

    r ./= mesh.volume
    
    return r
end

R(u3)

m = cartesian_mesh(50,50);

u = Field(m)
set_dirichlet_patch!(u, "left", 1.0);
set_dirichlet_patch!(u, "bottom", 0.0);
set_neumann_patch!(u, "right");
set_neumann_patch!(u, "top");

Δt = 0.01
for t = 0 : Δt : 2
    u ← u - Δt * R(u)
end

plot(reshape(u.values, 50,50)', st=:contour, color=:lightrainbow, fill=true, aspect_ratio=:equal)

using SparseArrays;

struct Equation
    A
    x
    b
end

# Scitani a odcitani rovnic 
for op ∈ [:+, :-]
    @eval begin
        function ($op)(e1::Equation, e2::Equation)
            @assert e1.x == e2.x
            Equation(($op)(e1.A, e2.A), e1.x, $(op)(e1.b, e2.b))
        end
    end
end

# Pricteni/odecteni pole o spravne velikosti 
for op ∈ [:+, :-]
    @eval begin
        function ($op)(e1::Equation, f::Array{Float64,1})
            @assert size(e1.x) == size(f)
            Equation(copy(e1.A), e1.x, $(op)(e1.b, f))
        end
    end
end

# Vynasobeni rovnice konstantou
*(a::Float64, eq::Equation) = Equation(a*eq.A, eq.x, a*eq.b);

function solve!(eqn::Equation)
    x = eqn.A \ eqn.b
    for i in eachindex(eqn.x)
        eqn.x[i] = -x[i]
    end
end;

function ddt(f::Field,Δt)
    n = length(f.values)
    Equation(I/Δt, f.values, -f.values/Δt)
end;    

ddt(u3,0.1)

u3.values = zeros(9)
for t=1:10
    solve!(ddt(u3,0.1)-ones(9))
    print(u3.values[1], "\t")
end

function ∇(v::Vector, u::Field)
    A = spzeros(length(u.values),length(u.values))
    b = zero(u.values)
    
    msh = u.mesh
    
    for f in internal_faces(msh)
        o = msh.owner[f]
        n = msh.neighbor[f]
        S = msh.surface[f]
        uo = u.values[o]
        un = u.values[n]
    
        ϕ = dot(v,S)
        α = max(ϕ, 0.0)
        β = min(ϕ, 0.0)
        
        A[o,o] += α / msh.volume[o]
        A[o,n] += β / msh.volume[o]
            
        A[n,o] -= α / msh.volume[n]
        A[n,n] -= β / msh.volume[n]
    end

    for p in boundary_patches(msh)
        name = msh.patch[p].name
        bc = u.boundaries[name]
        
        for f in patch_faces(msh, p)
            o = msh.owner[f]
            S = msh.surface[f]
            
            ϕ = dot(v,S)
            α = max(ϕ, 0.0)
            β = min(ϕ, 0.0)
            
            # ub = c1 + c2*u[o], flux=alpha*u[o] + beta*ub
            c1, c2 = boundary_coeffs(bc, f)

            A[o,o] += (α + β*c2) / msh.volume[o]
            b[o]   += β*c1 / msh.volume[o]
        end
    end
    
    return Equation(A, u.values, b)
end     

Matrix(∇(v,u3).A)

u.values = zero(u.values)

Δt = 0.05
for t=0:Δt:0.5
    solve!( ddt(u,Δt) + ∇(v,u) ); 
end

plot(reshape(u.values, 50,50)', st=:contour, color=:lightrainbow, fill=true, aspect_ratio=:equal)

function Δ(μ::Float64, u::Field)
    mesh = u.mesh
    
    A = spzeros(length(u.values),length(u.values))
    b = zero(u.values)

    
    for f in internal_faces(mesh)
        o, n = mesh.owner[f], mesh.neighbor[f]
        xo, xn = mesh.centre[o], mesh.centre[n]
        S = mesh.surface[f]
        
        g = μ * norm(S) / norm(xn - xo)

        A[o,o] -= g / mesh.volume[o]
        A[o,n] += g / mesh.volume[o]
            
        A[n,o] += g / mesh.volume[n]
        A[n,n] -= g / mesh.volume[n]
    end
    
    
    for p in boundary_patches(mesh)
        name = mesh.patch[p].name
        bc = u.boundaries[name]
        
        for f in patch_faces(mesh, p)
            o = mesh.owner[f]
            S = mesh.surface[f]
            xo, xf = mesh.centre[o], mesh.facecentre[f]
            c1, c2 = boundary_coeffs(bc, f)
            
            # normal derivative at boundary dudn = (ub - uin) / delta
            # with ub = c[1] + c[2]*uin, i.e.
            # dudn = c[1]/delta + (c[2]-1)/delta*uin
            
            δ = dot(xf - xo, S)/norm(S)
            A[o,o] += μ * (c2-1)/δ * norm(S) / mesh.volume[o]
            b[o] += μ * c1/δ * norm(S) / mesh.volume[o]
        end
    end

    Equation(A, u.values, b)
end


D=Matrix(Δ(1.0,u3).A)

A = Matrix((ddt(u3,Δt) + ∇(v,u3) - Δ(1.e-1,u3)).A)

u.values = zero(u.values)

Δt = 0.1
for t=0:Δt:0.5
    solve!( ddt(u,Δt) + ∇(v,u) - Δ(1.e-1,u) ); 
end

plot(reshape(u.values, 50,50)', st=:contour, color=:lightrainbow, fill=true, aspect_ratio=:equal)
