using SparseArrays
using LinearAlgebra

struct Equation{T}
    A :: SparseMatrixCSC{Scalar,Label}
    x :: Array{T,1}
    b :: Array{T,1}
end


Ac(eq::Equation{T}) where {T} = Array(diag(eq.A));
H(eq::Equation{T}) where {T} =  -eq.b .+ Ac(eq) .* eq.x .-  eq.A * eq.x;


function *(A::SparseMatrixCSC{Scalar, Label}, x::VectorList)
    n,m = size(A)
    @assert length(x) == m

    y = zeros(Vector, n)
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
        function ($op)(e1::Equation{T}, e2::Equation{T}) where {T}
            @assert e1.x == e2.x
            Equation{T}(($op)(e1.A, e2.A), e1.x, $(op)(e1.b, e2.b))
        end
    end
end


# Pricteni/odecteni pole o spravne velikosti 
for op ∈ [:+, :-]
    @eval begin
        function ($op)(e1::Equation{T}, f::Array{T,1}) where {T}
            @assert size(e1.x) == size(f)
            Equation{T}(copy(e1.A), e1.x, $(op)(e1.b, f))
        end
    end
end

# Vynasobeni rovnice konstantou
*(a::Scalar, eq::Equation{T}) where {T} = Equation{T}(a*eq.A, eq.x, a*eq.b);

# Unarni minus
-(eq::Equation{T}) where {T} = Equation{T}(-eq.A, eq.x, -eq.b);


function solve!(eqn::Equation{Scalar})
    x = eqn.A \ eqn.b
    for i in eachindex(eqn.x)
        eqn.x[i] = -x[i]
    end
end;


function solve!(eqn::Equation{Vector})
    rhs = zeros(length(eqn.b),2)
    for i in eachindex(eqn.b)
        rhs[i,1], rhs[i,2] = eqn.b[i][1], eqn.b[i][2]
    end
    x = eqn.A \ rhs
    for i in eachindex(eqn.x)
        eqn.x[i] = -Vector(x[i,1], x[i,2])
    end
end;


function relax!(eqn::Equation{T}, α) where {T}
    D = diag(eqn.A)
    for i=1:length(D) 
        eqn.A[i,i] /= α
    end
    eqn.b[:] -= (1-α)/α * D .* eqn.x
end;


# Aproximace casove derivace
function ddt(f::Field{T}, Δt) where {T} 
    n = length(f.values)
    Equation{T}(sparse(I/Δt,n,n), f.values, -f.values/Δt)
end;    


# Laplaceuv operator
function Δ(μ::Array{Float64,1}, u::Field{T}) where {T}
    mesh = u.mesh
    
    A = spzeros(length(u.values),length(u.values))
    b = zeros(T,length(u.values))

    
    for f in internal_faces(mesh)
        o, n = mesh.owner[f], mesh.neighbor[f]
        xo, xn = mesh.centre[o], mesh.centre[n]
        S = mesh.surface[f]
        
        μf = (μ[o] + μ[n]) / 2
        g = μf * norm(S) / norm(xn - xo)
        
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
            μf = μ[o]

            A[o,o] += μf * (c2-1)/δ * norm(S) / mesh.volume[o]
            b[o] += μf * c1/δ * norm(S) / mesh.volume[o]
        end
    end

    return Equation{T}(A, u.values, b)

end


Δ(μ::Float64, u::Field{T}) where {T} = Δ(μ*ones(length(u.values)), u::Field{T})


function ∂(u::Field{T},dir) where {T}
    A = spzeros(length(u.values),length(u.values))
    b = zeros(T, length(u.values))
    
    mesh = u.mesh
    
    for f in internal_faces(mesh)
        o = mesh.owner[f]
        n = mesh.neighbor[f]
        uo = u[o]
        un = u[n]
    
        ϕ = mesh.surface[f][dir]
        α = ϕ / 2
        β = ϕ / 2
        
        A[o,o] += α / mesh.volume[o]
        A[o,n] += β / mesh.volume[o]
            
        A[n,o] -= α / mesh.volume[o]
        A[n,n] -= β / mesh.volume[o]
    end

    for p in boundary_patches(mesh)
        name = mesh.patch[p].name
        bc = u.boundaries[name]
        
        for f in patch_faces(mesh, p)
            o = mesh.owner[f]
            ϕ = mesh.surface[f][dir]
            α = ϕ / 2
            β = ϕ / 2

            c1, c2 = boundary_coeffs(bc, f)  # tj. ub = c[1] + c[2]*uin
            A[o,o] += (α + β*c2) / mesh.volume[o]
            b[o]   += β*c1 / mesh.volume[o]
        end
    end

    return Equation{T}(A, u.values, b)
end     


∂x(u::Field{T}) where {T} = ∂(u::Field{T},1) 
∂y(u::Field{T}) where {T} = ∂(u::Field{T},2) 

# Explicitni vypocet divergence
function div(u::VectorField)
    mesh = u.mesh
    r = zeros(length(cells(mesh)))
    
    for f in internal_faces(mesh)
        o, n = mesh.owner[f], mesh.neighbor[f]
        uf = (u[o] + u[n]) / 2.0
        flux = dot(uf, mesh.surface[f])
        r[o] += flux
        r[n] -= flux
    end

    for p in boundary_patches(mesh)
        name = mesh.patch[p].name
        bc = u.boundaries[name]
        
        for f in patch_faces(mesh, p)
            ub = boundary_value(bc, f)
            flux = dot(ub, mesh.surface[f])
            r[mesh.owner[f]] += flux
        end
    end
    
    for c in cells(mesh)
        r[c] /= mesh.volume[c]
    end
    
    return r
end


# Explicitni vypocet gradientu
function grad(u::ScalarField)
    mesh = u.mesh
    r = zeros(Vector,length(cells(mesh)))
    
    for f in internal_faces(mesh)
        o, n = mesh.owner[f], mesh.neighbor[f]
        uf = (u[o] + u[n]) / 2.0
        flux = uf * mesh.surface[f]
        r[o] += flux
        r[n] -= flux
    end
    
    for p in boundary_patches(mesh)
        name = mesh.patch[p].name
        bc = u.boundaries[name]
        
        for f in patch_faces(mesh, p)
            o = mesh.owner[f]
            ub = boundary_value(bc, f)
            flux = ub * mesh.surface[f]
            r[o] += flux
        end
    end

    for c in cells(mesh)
        r[c] /= mesh.volume[c]
    end
    
    return r
end

