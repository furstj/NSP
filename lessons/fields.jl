import Base.-, Base.+, Base.*, Base./
import Base.getindex, Base.setindex!

struct Field{T}
    values :: Array{T,1}
    mesh
    boundaries
end

struct BoundaryCondition{T}
    c1    :: Array{T,1}
    c2    :: Array{Scalar,1}
    field :: Field{T}
    patch :: Label
end


getindex(f::Field{T}, i::Int) where {T} = f.values[i]
setindex!(f::Field{T}, v::T, i::Int) where {T} = begin f.values[i] = v end

const ScalarField = Field{Scalar}
ScalarField(m::Mesh) = ScalarField( zeros(Scalar,length(m.volume)), m, Dict{String,Any}());


const VectorField = Field{Vector}
VectorField(m::Mesh) = VectorField( zeros(Vector,length(m.volume)), m, Dict{String,Any}());


for op ∈ [:+, :-, :*, :/]
    @eval begin
        function ($op)(a::Field{T}, b::Field{T}) where {T}
            assert(a.mesh == b.mesh)
            Field{T}(broadcast($op,a.values,b.values), a.mesh, Dict{String,Any}());
        end
    end
end


for op ∈ [:+, :-, :*, :/]
    @eval begin
        function ($op)(a::Field{T}, b::T) where {T}
            Field{T}(broadcast($op,a.values,b), a.mesh, Dict{String,Any}());
        end
    end
end      


for op ∈ [:+, :-, :*, :/]
    @eval begin
        function ($op)(a::Field{T}, b::Array{T,1}) where {T}
            @assert size(a.values) == size(b) 
            Field(broadcast($op,a.values,b), a.mesh, Dict{String,Any}());
        end
    end
end      

*(a::Number, f::Field{T}) where {T} = Field(a*f.values, f.mesh, Dict{String,Any}());
-(a::Field{T}) where {T} = Field{T}(-a.values, a.mesh);


function ←(a::Field{T}, b::Field{T}) where {T}
    @assert a.mesh == b.mesh
    for i in eachindex(a.values)
        @inbounds a.values[i] = b.values[i]
    end
end

function ←(a::Field{T}, b::Float64) where {T}
    a.values[:] .= b
end

function ←(a::VectorField, b::Vector) where {T}
    for i in eachindex(a.values)
        @inbounds a.values[i] = b
    end
end

function ←(a::Field{T}, b::Array{T,1}) where {T}
    @assert length(a.values) == length(b)
    for i in eachindex(a.values)
        @inbounds a.values[i] = b[i]
    end
end


function boundary_coeffs(bc::BoundaryCondition{T}, f::Label) where {T}
    field = bc.field
    mesh = field.mesh
    start = first(mesh.patch[bc.patch].range)
    i = f + 1 - start
    
    c1 = bc.c1[i]
    c2 = bc.c2[i]

    return c1, c2
end


function boundary_value(bc::BoundaryCondition{T}, f::Label) where {T}
    field = bc.field
    mesh = field.mesh
    u1 = field.values[mesh.owner[f]]
    c1, c2 = boundary_coeffs(bc, f)
    return c1 + c2*u1
end


function set_dirichlet_patch!(field::Field{T}, name::String, val::T) where {T}
    mesh = field.mesh
    p = patch_by_name(mesh, name)
    c1 = [val for f in patch_faces(mesh, p)]
    c2 = zeros(Scalar, length(c1))
    field.boundaries[name] = BoundaryCondition(c1, c2, field, p)
end


function set_neumann_patch!(field::Field{T}, name::String) where {T}
    mesh = field.mesh
    p = patch_by_name(mesh, name)
    c2 = ones(length(patch_faces(mesh, p)))
    c1 = zeros(T, length(c2))
    field.boundaries[name] = BoundaryCondition(c1, c2, field, p)
end



function interpolate_to_points(f::Field{T}) where{T}
    mesh = f.mesh
    fn = zeros(T, length(mesh.point))
    vn = zeros(Scalar, length(mesh.point))
    for c in cells(mesh)
        for pid in mesh.cell2point[c]
            fn[pid] += f[c] * mesh.volume[c]
            vn[pid] += mesh.volume[c]
        end
    end
    return fn ./ vn
end


# Vykresleni pomoci matplotlib
function plot_contour(f::ScalarField, kwargs...)
    mesh = f.mesh
    x = [p[1] for p in mesh.point]
    y = [p[2] for p in mesh.point]
    tri  = triangles(mesh) .- 1
    fn = interpolate_to_points(f)

    PyPlot.tricontour(x, y, tri, fn; kwargs...)
end


function plot_contour(f::VectorField; kwargs...)
    mag = ScalarField(f.mesh)
    mag ← [norm(v) for v in f.values]
    plot_contour(mag; kwargs...)
end


function plot_contourf(f::ScalarField; kwargs...)
    mesh = f.mesh
    x = [p[1] for p in mesh.point]
    y = [p[2] for p in mesh.point]
    tri  = triangles(mesh)
    for i=1:length(tri)
        tri[i] .-= 1
    end
    fn = interpolate_to_points(f)

    PyPlot.tricontourf(x, y, tri, fn; kwargs...)
end


function plot_contourf(f::VectorField; kwargs...)
    mag = ScalarField(f.mesh)
    mag ← [norm(v) for v in f.values]
    plot_contourf(mag; kwargs...)
end


function plot_arrows(f::VectorField; kwargs...)
    mesh = f.mesh
    x = [c[1] for c in mesh.centre]
    y = [c[2] for c in mesh.centre]
    u = [v[1] for v in f.values]
    v = [v[2] for v in f.values]

    PyPlot.quiver(x, y, u, v; kwargs...)
end
