type Field{T}
    mesh
    value :: Array{T,1}
    boundaryfield :: Dict{String, Array{T,1}}
    boundarycondition
end

function Field(T::Type, m::CartesianMesh, bc)
    internalfield = zeros(T, m.nx*m.ny)
    boundaryfield = Dict(
                         "top" => zeros(T, m.nx),
                         "bottom" => zeros(T, m.nx),
                         "left" => zeros(T, m.ny),
                         "right" => zeros(T, m.ny)
    )
    Field{T}(m, internalfield, boundaryfield, bc)
end
    
typealias ScalarField Field{Float64};
typealias VectorField Field{Vec2d};

ScalarField(m :: CartesianMesh, bc) = Field(Float64, m, bc);
VectorField(m :: CartesianMesh, bc) = Field(Vec2d, m, bc);

asmatrix(f::ScalarField) = reshape(f.value, (f.mesh.nx,f.mesh.ny) );

function asmatrix(f::VectorField)
    m = zeros(f.mesh.nx, f.mesh.ny, 2)
    for j=1:f.mesh.ny, i=1:f.mesh.nx
        m[i,j,:] = f.value[i+(j-1)*f.mesh.nx]
    end
    return m
end;


type FaceField{T}
    mesh
    internalfield :: Array{T,1}
    boundaryfield :: Dict{String, Array{T,1}}
end

typealias ScalarFaceField FaceField{Float64}


#
# α uf + β (uf-u1) / Δ = g   => (αΔ + β) uf = gΔ + β u1
#
type Robin
    α
    β
    g
end

# uf = c1 + c2 * u1 = ...
bndcoeffs(uin, Δ, bc::Robin) = ( bc.g*Δ / (bc.α*Δ + bc.β), bc.β / (bc.α*Δ + bc.β) );

function bndvalue(uin, Δ, bc::Robin)
    c1, c2 = bndcoeffs(uin, Δ, bc)
    c1 + c2*uin
end

# dudn = (uf - u1) / Δ = g/(αΔ + β) - α/(β + αΔ) * u1
ddncoeffs(uin, Δ, bc::Robin) = ( bc.g / (bc.α*Δ + bc.β), -bc.α / (bc.α*Δ + bc.β) );


Dirichlet(value) = Robin(1.0, 0.0, value);
Neumann(value) = Robin(0.0, 1.0, value);


function correctboundary!{T}(u::Field{T})
    mesh = u.mesh

    for patch in patches(mesh)
        bc = u.boundarycondition[ name(patch) ]
        up = u.boundaryfield[ name(patch) ]

        for fi in zip(boundaryfaces(patch), eachindex(up))
            f,i = fi
            owner = f.owner
            co = cell(mesh, owner)
            Δ = norm(f.x - co.x)
            up[i] = bndvalue(u.value[owner], Δ, bc)
        end
    end
end    
