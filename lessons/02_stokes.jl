
using PyPlot
using StaticArrays
import Base.-, Base.+, Base.*

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

m2 = CartesianMesh(2,2)

for c in cells(m2)
    println(c)
end
println()
println(cell(m2,3))

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

for f in internalfaces(m2)
    println(f)
end

immutable Patch
    name
    faces
end

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

name(p::Patch) = p.name
boundaryfaces(p::Patch) = p.faces()

for p in patches(m2)
    println("PATCH ", p.name)
    for f in boundaryfaces(p)
        println(f)
    end
    println()
end

type Equation
    A
    x
    b
end

-(eq::Equation) = Equation(-eq.A, eq.x, -eq.b);
-(eq::Equation, b::Array{Float64,1}) = Equation(eq.A, eq.x, eq.b + b);
+(eq::Equation, b::Array{Float64,1}) = Equation(eq.A, eq.x, eq.b - b);
*(a::Float64, eq::Equation) = Equation(a*eq.A, eq.x, a*eq.b);

function relax!(eqn::Equation, α)
    D = diag(eqn.A)
    for i=1:length(D) 
        eqn.A[i,i] /= α
    end
    eqn.b += (1-α)/α * D .* eqn.x
end;

function solve!(eqn::Equation)
    x = eqn.A \ eqn.b
    for i in eachindex(eqn.x)
        eqn.x[i] = x[i]
    end
end;

type Field
    value
    mesh
    boundarycondition
end

Field(m :: CartesianMesh, bc) = Field( zeros(m.nx*m.ny), m, bc);

asmatrix(f::Field) = reshape(f.value, (f.mesh.nx,f.mesh.ny) );

type Dirichlet
    value
end

type Neumann
    value
end

bndvalue(uin, Δ, bc::Dirichlet) = bc.value;
bndvalue(uin, Δ, bc::Neumann) = uin + Δ * bc.value;

ddncoeffs(Δ, bc::Dirichlet) = (-1/Δ, bc.value/Δ);
ddncoeffs(Δ, bc::Neumann) = (0, bc.value);


msh3 = CartesianMesh(3,3)

function createfields(msh)
    u = Field(msh, Dict( "top"=>Dirichlet(1), "left"=>Dirichlet(0), "right"=>Dirichlet(0), "bottom"=>Dirichlet(0)) );
    v = Field(msh, Dict( "top"=>Dirichlet(0), "left"=>Dirichlet(0), "right"=>Dirichlet(0), "bottom"=>Dirichlet(0)) );
    p = Field(msh, Dict( "top"=>Neumann(0), "left"=>Neumann(0), "right"=>Neumann(0), "bottom"=>Neumann(0)) );
    return (u,v,p)
end

function laplace(ν, u)
    mesh = u.mesh
    dims = (mesh.nx, mesh.ny)
    n  = prod(dims)
    A = spzeros(n,n)
    b = zeros(n)
    
    for f in internalfaces(mesh)
        owner, neigh = f.owner, f.neigh
        
        νf = (ν[owner]+ν[neigh]) / 2.0
        
        co = cell(mesh, owner)
        cn = cell(mesh, neigh)
        
        g = νf * norm(f.s) / norm(cn.x-co.x)

        A[owner, owner] -= g / co.vol
        A[owner, neigh] += g / co.vol
            
        A[neigh, owner] += g / cn.vol
        A[neigh, neigh] -= g / cn.vol
    end
    
    for patch in patches(mesh)

        bc = u.boundarycondition[ name(patch) ]
        for f in boundaryfaces(patch)
            owner = f.owner
            co = cell(mesh, owner)
            
            Δ = norm(f.x - co.x)
            νf = ν[owner]
            
            a,val = ddncoeffs(Δ, bc)
            A[owner,owner] += νf * a * norm(f.s) / co.vol
            b[owner] -= νf * val * norm(f.s) / co.vol
        end
    end
    
    Equation(A, u.value, b)
end

laplace(u) = laplace(ones(u.value), u)

function ddxi(p, dir)
    mesh = p.mesh
    dp = zeros(p.value)
    
    for f in internalfaces(mesh)
        owner = f.owner
        neigh = f.neigh
        
        pf = (p.value[owner] + p.value[neigh]) / 2.0
        
        dp[owner] += pf * f.s[dir]
        dp[neigh] -= pf * f.s[dir]
    end

    for patch in patches(mesh)
        bc = p.boundarycondition[ name(patch) ]
        for f in boundaryfaces(patch)
            owner = f.owner
            co = cell(mesh, owner)
            Δ = norm(f.x - co.x)
            pf = bndvalue(p.value[owner], Δ, bc)
            dp[owner] += pf * f.s[dir]
        end
    end
    
    for c in cells(mesh)
        dp[c.id] /= c.vol
    end
    
    return dp
end;

ddx(p) = ddxi(p,1);
ddy(p) = ddxi(p,2);

uu,_,_ = createfields(msh3);

uuEqn = laplace(uu)
full(uuEqn.A)

uuEqn.x

uuEqn.b

solve!(uuEqn);
uu.value

msh100 = CartesianMesh(100,100);

u100,_,_ = createfields(msh100);

solve!( laplace(u100) - ones(100*100));

contourf(asmatrix(u100)'); colorbar();

minimum(u100.value), maximum(u100.value)

ν = 0.01

msh10 = CartesianMesh(10,10);
u,v,p = createfields(msh10);

uEqn = (-ν) * laplace(u);
vEqn = (-ν) * laplace(v);

solve!(uEqn + ddx(p));
solve!(vEqn + ddy(p));

Ac(eq::Equation) = diag(eq.A);
H(eq::Equation) =  eq.b + Ac(eq) .* eq.x -  eq.A * eq.x;

ra = 1 ./ Ac(uEqn);

uBar = Field(ra .* H(uEqn), u.mesh, u.boundarycondition);
vBar = Field(ra .* H(vEqn), v.mesh, v.boundarycondition);

pEqn = laplace(ra, p) - (ddx(uBar) + ddy(vBar));
pEqn.A[1,1] += 1;

solve!(pEqn)

contourf(asmatrix(p)'); colorbar();

u.value = uBar.value - ra .* ddx(p);
v.value = vBar.value - ra .* ddy(p);

quiver(asmatrix(u)', asmatrix(v)');

reshape(ddx(u) + ddy(v), (u.mesh.nx,u.mesh.ny))

msh25 = CartesianMesh(25,25)

u,v,p = createfields(msh25)

α = 0.7
β = 0.3
ν = 1.0

for iter = 0:50
    
    uOld, vOld, pOld = copy(u.value), copy(v.value), copy(p.value)
    
    uEqn = (-ν)*laplace(u)
    vEqn = (-ν)*laplace(v)

    relax!(uEqn, α)
    relax!(vEqn, α)
    
    solve!(uEqn + ddx(p))
    solve!(vEqn + ddy(p))
    
    ra = 1 ./ Ac(uEqn);
    
    uBar = Field(ra .* H(uEqn), u.mesh, u.boundarycondition);
    vBar = Field(ra .* H(vEqn), v.mesh, v.boundarycondition);
    
    pEqn = laplace(ra, p) - (ddx(uBar) + ddy(vBar));
    pEqn.A[1,1] -= 1/(p.mesh.Δx*p.mesh.Δy)
    #pEqn.A[1,1] *= 2
    solve!(pEqn)
    
    p.value = (1-β) * pOld + β * p.value
    u.value = uBar.value - ra .* ddx(p)
    v.value = vBar.value - ra .* ddy(p)
    
    if rem(iter,5)==0
        nxny = u.mesh.nx*u.mesh.ny
        pRez = norm(pOld - p.value) / nxny
        uRez = norm(uOld - u.value) / nxny
        vRez = norm(vOld - v.value) / nxny
        println(iter, "\t", pRez, "\t", uRez, "\t", vRez)
    end
end

quiver(asmatrix(u)', asmatrix(v)');

contourf(asmatrix(p)'); colorbar();

contourf(√(asmatrix(u).^2 + asmatrix(v).^2)'); colorbar();




