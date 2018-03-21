
using PyPlot
using StaticArrays
import Base.-, Base.+, Base.*, Base./
using ForwardDiff

const Vec2d = SVector{2,Float64};
const VecList = Array{Vec2d,1};

immutable Cell
    id  :: Int
    x   :: Vec2d
    vol :: Float64
end

const CellList = Array{Cell,1};

immutable Face
    id    :: Int
    x     :: Vec2d
    s     :: Vec2d
    owner :: Int
    neigh :: Int
end;

const FaceList = Array{Face,1};

const PatchDict = Dict{String,FaceList};

type Mesh
    points  :: VecList
    faces   :: FaceList
    cells   :: CellList
    patches :: PatchDict
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
    
    return Mesh(points, faces, cells, patches)
end

m2 = cartesian_mesh(2,2);

for p in m2.points
    println(p)
end

for f in m2.faces
    println(f)
end

for c in m2.cells
    println(c)
end

for (name,faces) in m2.patches
    println("PATCH ", name)
    for f in faces
        println(f)
    end
    println()
end

type Field
    values
    mesh
    boundaries
end

Field(m::Mesh) = Field( zeros(length(m.cells)), m, Dict{String,Any}());
Field(values, m::Mesh) = Field(values, m, Dict{String,Any}());

for op ∈ [:+, :-, :*, :/]
    @eval begin
        function ($op)(a::Field, b::Field)
            assert(a.mesh == b.mesh)
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
            assert( size(a.values) == size(b) )
            Field(broadcast($op,a.values,b), a.mesh)
        end
    end
end      

-(a::Field) = Field(-a.values, a.mesh);

c2 = a2 / 3;
c2.values

function ←(a::Field, b::Field)
    assert(a.mesh == b.mesh)
    a.values = copy(b.values)
end

b2.values = 3*ones(4)
c2 ← b2;
d2 = b2;
b2.values[1] = 5;
println(c2.values)
println(d2.values)

function boundary_value(patch, i::Int)
    c = boundary_coeffs(patch,i)
    u1 = patch.field.values[patch.faces[i].owner]
    return c[1] + c[2]*u1
end;

type DirichletPatch
    value :: Array{Float64,1}
    field :: Field
    faces :: FaceList
end

boundary_coeffs(p::DirichletPatch, i::Int) = (p.value[i], 0.0)

# Pomocna funkce pro nastaveni Dirichletovy okrajove podminky s jednou zadanou hodnotou
function set_dirichlet_patch!(u::Field, name::String, val::Float64)
    value = [val for i=1:length(u.mesh.patches[name])]
    u.boundaries[name] = DirichletPatch(value, u, u.mesh.patches[name])
end;

type NeumannPatch
    deriv :: Array{Float64,1}
    field :: Field
    faces :: FaceList
end


function boundary_coeffs(p::NeumannPatch, i::Int)
    u = p.field
    m = u.mesh
    f = p.faces[i]
    δ = dot(f.x - m.cells[f.owner].x, f.s) / norm(f.s)
    return (p.deriv[i] * δ, 1.0)
end

function set_neumann_patch!(u::Field, name::String, der::Float64)
    u.boundaries[name] = NeumannPatch([der for i=1:length(u.mesh.patches[name])], u, u.mesh.patches[name])
end;


m3 = cartesian_mesh(3,3);

u3 = Field(m3);
set_dirichlet_patch!(u3, "left", 1.0);
set_dirichlet_patch!(u3, "bottom", 0.0);
set_neumann_patch!(u3, "right", 0.0);
set_neumann_patch!(u3, "top", 0.0);

v = Vec2d(1, 0.5);

# Tok typu upwind
H(ul, ur, s) = dot(v,s)>0 ? dot(v,s)*ul : dot(v,s)*ur

function R(u::Field)
    m = u.mesh
    r = zeros(u.values)    
    
    for f ∈ m.faces
        flux = H( u.values[f.owner], u.values[f.neigh], f.s)
        r[f.owner] += flux
        r[f.neigh] -= flux
    end
    
    for (name,faces) ∈ m.patches
        bc = u.boundaries[name]
        for (i,f) ∈ enumerate(bc.faces)
            ub = boundary_value(bc, i)
            flux = H( u.values[f.owner], ub, f.s)
            r[f.owner] += flux
        end
    end
    
    for c in m.cells
        r[c.id] /= c.vol
    end
    
    return r
end

R(u3)

m = cartesian_mesh(50,50);

u = Field(m)
set_dirichlet_patch!(u, "left", 1.0);
set_dirichlet_patch!(u, "bottom", 0.0);
set_neumann_patch!(u, "right", 0.0);
set_neumann_patch!(u, "top", 0.0);

Δt = 0.01
for t = 0 : Δt : 2
    u ← u - Δt * R(u)
end

contourf(reshape(u.values,50,50)'); 
colorbar(); axis("equal");

type Equation
    A
    x
    b
end

# Scitani a odcitani rovnic 
for op ∈ [:+, :-]
    @eval begin
        function ($op)(e1::Equation, e2::Equation)
            assert(e1.x == e2.x)
            Equation(($op)(e1.A, e2.A), e1.x, $(op)(e1.b, e2.b))
        end
    end
end

# Pricteni/odecteni pole o spravne velikosti 
for op ∈ [:+, :-]
    @eval begin
        function ($op)(e1::Equation, f::Array{Float64,1})
            assert(size(e1.x) == size(f))
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
    Equation(speye(n)/Δt, f.values, -f.values/Δt)
end;    

ddt(u3,0.1)

u3.values = zeros(9)
for t=1:10
    solve!(ddt(u3,0.1)-ones(9))
    print(u3.values[1], "\t")
end

function ∇(v::Vec2d, u::Field)
    A = spzeros(length(u.values),length(u.values))
    b = zeros(length(u.values))
    
    msh = u.mesh
    
    for f ∈ msh.faces
        o = f.owner
        n = f.neigh
        uo = u.values[o]
        un = u.values[n]
    
        ϕ = dot(v,f.s)
        α = max(ϕ, 0.0)
        β = min(ϕ, 0.0)
        
        A[o,o] += α / msh.cells[o].vol
        A[o,n] += β / msh.cells[o].vol
            
        A[n,o] -= α / msh.cells[n].vol
        A[n,n] -= β / msh.cells[n].vol
    end

    for (name,faces) ∈ msh.patches
        bc = u.boundaries[name]
        for (i,f) ∈ enumerate(bc.faces)
            o = f.owner
            coeff = boundary_coeffs(bc, i)  # tj. ub = coef[1] + coef[2]*uin
            ϕ = dot(v,f.s)
            α = max(ϕ, 0.0)
            β = min(ϕ, 0.0)

            c = boundary_coeffs(bc, i)  # tj. ub = c[1] + c[2]*uin
            A[o,o] += (α + β*c[2]) / msh.cells[o].vol
            b[o]   += β*c[1] / msh.cells[o].vol
        end
    end

    
    return Equation(A, u.values, b)
end     

full(∇(v,u3).A)

u.values = zeros(u.values)

Δt = 0.05
for t=0:Δt:0.5
    solve!( ddt(u,Δt) + ∇(v,u) ); 
end

contourf(reshape(u.values,50,50)'); 
colorbar(); axis("equal");

function Δ(μ::Float64, u::Field)
    mesh = u.mesh
    
    A = spzeros(length(u.values),length(u.values))
    b = zeros(length(u.values))

    
    for f in mesh.faces
        o, n = f.owner, f.neigh
        co, cn = mesh.cells[o], mesh.cells[n]
        
        g = μ * norm(f.s) / norm(cn.x-co.x)

        A[o,o] -= g / co.vol
        A[o,n] += g / co.vol
            
        A[n,o] += g / cn.vol
        A[n,n] -= g / cn.vol
    end
    
    for (name,faces) ∈ mesh.patches
        bc = u.boundaries[name]
        for (i,f) ∈ enumerate(bc.faces)
            o = f.owner
            co = mesh.cells[o]
    
            # normal derivative at boundary dudn = (ub - uin) / delta
            # with ub = c[1] + c[2]*uin, i.e.
            # dudn = c[1]/delta + (c[2]-1)/delta*uin
            
            c = boundary_coeffs(bc, i)  # tj. ub = c[1] + c[2]*uin            
            δ = dot(f.x - co.x, f.s) / norm(f.s)

            A[o,o] += μ * (c[2]-1)/δ * norm(f.s) / co.vol
            b[o] += μ * c[1]/δ * norm(f.s) / co.vol
        end
    end
    
    Equation(A, u.values, b)
end


D=full(Δ(1.0,u3).A)

A = full((ddt(u3,Δt) + ∇(v,u3) - Δ(1.e-1,u3)).A)

sum(A,2)'

u.values = zeros(u.values)

Δt = 0.1
for t=0:Δt:0.5
    solve!( ddt(u,Δt) + ∇(v,u) - Δ(1.e-1,u) ); 
end

contourf(reshape(u.values,50,50)'); 
colorbar(); axis("equal");
