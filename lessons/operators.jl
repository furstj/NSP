type Equation{T}
    A :: SparseMatrixCSC{Float64,Int64}
    x :: Array{T,1}
    b :: Array{T,1}
end

function *(A::SparseMatrixCSC{Float64,Int64}, x::Array{Vec2d,1})
    n,m = size(A)
    assert(length(x) == m)

    y = zeros(Vec2d, n)
    rows = rowvals(A)
    vals = nonzeros(A)
    for c=1:m
        for k ∈ nzrange(A,c)
            r = rows[k]
            a = vals[k]
            y[r] += a*x[c]
        end
    end
    return y
end

# Scitani a odcitani rovnic 
for op ∈ [:+, :-]
    @eval begin
        function ($op){T}(e1::Equation{T}, e2::Equation{T})
            assert(e1.x == e2.x)
            Equation{T}(($op)(e1.A, e2.A), e1.x, $(op)(e1.b, e2.b))
        end
    end
end

# Pricteni/odecteni pole o spravne velikosti 
for op ∈ [:+, :-]
    @eval begin
        function ($op){T}(e1::Equation{T}, f::Array{T,1})
            assert(size(e1.x) == size(f))
            Equation{T}(copy(e1.A), e1.x, $(op)(e1.b, f))
        end
    end
end

# Vynasobeni rovnice konstantou
*{T}(a::Float64, eq::Equation{T}) = Equation{T}(a*eq.A, eq.x, a*eq.b);

# Unarni minus
-{T}(eq::Equation{T}) = Equation{T}(-eq.A, eq.x, -eq.b);


function solve!(eqn::Equation{Float64})
    x = eqn.A \ eqn.b
    for i in eachindex(eqn.x)
        eqn.x[i] = -x[i]
    end
end;


function solve!(eqn::Equation{Vec2d})
    rhs = zeros(length(eqn.b),2)
    for i in eachindex(eqn.b)
        rhs[i,1], rhs[i,2] = eqn.b[i][1], eqn.b[i][2]
    end
    x = eqn.A \ rhs
    for i in eachindex(eqn.x)
        eqn.x[i] = -Vec2d(x[i,1], x[i,2])
    end
end;


function relax!{T}(eqn::Equation{T}, α)
    D = diag(eqn.A)
    for i=1:length(D) 
        eqn.A[i,i] /= α
    end
    eqn.b -= (1-α)/α * D .* eqn.x
end;


# Aproximace casove derivace
function ddt{T}(f::Field{T},Δt)
    n = length(f.values)
    Equation{T}(speye(n)/Δt, f.values, -f.values/Δt)
end;    

# Laplaceuv operator
function Δ{T}(μ::Array{Float64,1}, u::Field{T})
    mesh = u.mesh
    
    A = spzeros(length(u.values),length(u.values))
    b = zeros(T,length(u.values))

    
    for f in mesh.faces
        o, n = f.owner, f.neigh
        co, cn = mesh.cells[o], mesh.cells[n]

        μf = (μ[o] + μ[n]) / 2
        g = μf * norm(f.s) / norm(cn.x-co.x)

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
    
            c = boundary_coeffs(bc, i)  
            δ = dot(f.x - co.x, f.s) / norm(f.s)

            μf = μ[o]
            
            A[o,o] += μf * (c[2]-1)/δ * norm(f.s) / co.vol
            b[o] += μf * c[1]/δ * norm(f.s) / co.vol
        end
    end
    
    Equation{T}(A, u.values, b)
end


Δ{T}(μ::Float64, u::Field{T}) = Δ(μ*ones(length(u.values)), u::Field{T})


function ∂{T}(u::Field{T},dir)
    A = spzeros(length(u.values),length(u.values))
    b = zeros(u.values)
    
    mesh = u.mesh
    
    for f ∈ mesh.faces
        o = f.owner
        n = f.neigh
        uo = u.values[o]
        un = u.values[n]
    
        ϕ = f.s[dir]
        α = ϕ / 2
        β = ϕ / 2
        
        A[o,o] += α / mesh.cells[o].vol
        A[o,n] += β / mesh.cells[o].vol
            
        A[n,o] -= α / mesh.cells[n].vol
        A[n,n] -= β / mesh.cells[n].vol
    end

    for (name,faces) ∈ mesh.patches
        bc = u.boundaries[name]
        for (i,f) ∈ enumerate(bc.faces)
            o = f.owner
            coeff = boundary_coeffs(bc, i)  
            ϕ = f.s[dir]
            α = ϕ / 2
            β = ϕ / 2

            c = boundary_coeffs(bc, i)  # tj. ub = c[1] + c[2]*uin
            A[o,o] += (α + β*c[2]) / mesh.cells[o].vol
            b[o]   += β*c[1] / mesh.cells[o].vol
        end
    end

    
    return Equation{T}(A, u.values, b)
end     

∂x{T}(u::Field{T}) = ∂(u::Field{T},1)
∂y{T}(u::Field{T}) = ∂(u::Field{T},2)


# Explicitni vypocet divergence
function ∇(u::VectorField)
    m = u.mesh
    r = zeros(length(u.values))    
    
    for f ∈ m.faces
        o, n = f.owner, f.neigh
        uf = (u.values[o] + u.values[n]) / 2.0
        flux = dot(uf, f.s)
        r[f.owner] += flux
        r[f.neigh] -= flux
    end
    
    for (name,faces) ∈ m.patches
        bc = u.boundaries[name]
        for (i,f) ∈ enumerate(bc.faces)
            ub = boundary_value(bc, i)
            flux = dot(ub, f.s)
            r[f.owner] += flux
        end
    end
    
    for c in m.cells
        r[c.id] /= c.vol
    end
    
    return r
end

# Explicitni vypocet gradientu
function ∇(p::ScalarField)
    m = p.mesh
    r = zeros(Vec2d,length(u.values))    
    
    for f ∈ m.faces
        o, n = f.owner, f.neigh
        pf = (p.values[o] + p.values[n]) / 2.0
        flux = pf * f.s
        r[f.owner] += flux
        r[f.neigh] -= flux
    end
    
    for (name,faces) ∈ m.patches
        bc = p.boundaries[name]
        for (i,f) ∈ enumerate(bc.faces)
            pb = boundary_value(bc, i)
            flux = pb * f.s
            r[f.owner] += flux
        end
    end
    
    for c in m.cells
        r[c.id] /= c.vol
    end
    
    return r
end
