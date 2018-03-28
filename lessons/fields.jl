import Base.-, Base.+, Base.*, Base./
import Base.getindex, Base.setindex!

type Field{T}
    values :: Array{T,1}
    mesh
    boundaries
end

getindex{T}(f::Field{T}, i::Int) = f.values[i]
setindex!{T}(f::Field{T}, v::T, i::Int) = begin f.values[i] = v end

const ScalarField = Field{Float64}
ScalarField(m::Mesh) = ScalarField( zeros(Float64,length(m.cells)), m, Dict{String,Any}());


const VectorField = Field{Vec2d}
VectorField(m::Mesh) = VectorField( zeros(Vec2d,length(m.cells)), m, Dict{String,Any}());


for op ∈ [:+, :-, :*, :/]
    @eval begin
        function ($op){T}(a::Field{T}, b::Field{T})
            assert(a.mesh == b.mesh)
            Field{T}(broadcast($op,a.values,b.values), a.mesh, Dict{String,Any}());
        end
    end
end


for op ∈ [:+, :-, :*, :/]
    @eval begin
        function ($op){T}(a::Field{T}, b::T)
            Field{T}(broadcast($op,a.values,b), a.mesh, Dict{String,Any}());
        end
    end
end      


for op ∈ [:+, :-, :*, :/]
    @eval begin
        function ($op){T}(a::Field{T}, b::Array{T,1})
            assert( size(a.values) == size(b) )
            Field(broadcast($op,a.values,b), a.mesh, Dict{String,Any}());
        end
    end
end      

*{T}(a::Number, f::Field{T}) = Field(a*f.values, f.mesh, Dict{String,Any}());
-{T}(a::Field{T}) = Field{T}(-a.values, a.mesh);


function ←{T}(a::Field{T}, b::Field{T})
    assert(a.mesh == b.mesh)
    a.values = copy(b.values)
end


function boundary_value(patch, i::Int)
    c = boundary_coeffs(patch,i)
    u1 = patch.field.values[patch.faces[i].owner]
    return c[1] + c[2]*u1
end;


type DirichletPatch
    value 
    field 
    faces :: FaceList
end


boundary_coeffs(p::DirichletPatch, i::Int) = (p.value[i], 0.0)


# Pomocna funkce pro nastaveni Dirichletovy okrajove podminky s jednou zadanou hodnotou
function set_dirichlet_patch!{T}(u::Field{T}, name::String, val::T)
    value = [val for i=1:length(u.mesh.patches[name])]
    u.boundaries[name] = DirichletPatch(value, u, u.mesh.patches[name])
end;

type NeumannPatch
    deriv 
    field 
    faces :: FaceList
end


function boundary_coeffs(p::NeumannPatch, i::Int)
    u = p.field
    m = u.mesh
    f = p.faces[i]
    δ = dot(f.x - m.cells[f.owner].x, f.s) / norm(f.s)
    return (p.deriv[i] * δ, 1.0)
end

function set_neumann_patch!{T}(u::Field{T}, name::String, der)
    deriv = [der for i=1:length(u.mesh.patches[name])]
    u.boundaries[name] = NeumannPatch(deriv, u, u.mesh.patches[name])
end;


