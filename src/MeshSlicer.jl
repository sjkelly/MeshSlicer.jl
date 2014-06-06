module MeshSlicer

using ImmutableArrays
using DataStructures
import Base.push!
import Base.isapprox

export Bounds, Face, PolygonMesh, LineSegment, PolygonSlice,
       update!, rotate!, rotate

type Bounds{T<:Number}
    xmax::T
    ymax::T
    zmax::T
    xmin::T
    ymin::T
    zmin::T
end

type Face
    vertidx::Vector3{Int64}
    normal::Vector3{Float64}
end

type PolygonMesh
    bounds::Bounds{Float64}
    vertices::Array{Vector3{Float64}}
    faces::Array{Face}
end

#type LineSegment
#    start::Vector2{Float64}
#    finish::Vector2{Float64}
#    normal::Vector3{Float64}
#end

#type Polygon
#    segments::LinkedList{LineSegments}
#end

#type MeshSlice
#    polygons::Array{Polygon}
#    layer::Float64
#end

################################################################################
#
# MeshSlice
#
################################################################################

#function MeshSlice(mesh::PolygonMesh, height::Float64)

#    segmentlist = LineSegment[]

#    for face in mesh.faces
#        zmin, zmax = extrema([face.vertices[j].e3 for j=1:3])
#        if height > zmax
#            break
#        elseif zmin <= height
#            seg = LineSegment(face, height)
#            if !is(seg, nothing)
#                push!(segmentlist, seg)
#            end
#        end
#    end

#    return MeshSlice(segmentlist, height)
#end

#function MeshSlice(mesh::PolygonMesh, heights::Array{Float64})
#    # slice a mesh at heights given in a
#    # monotonically increasing array of heights

#    slices = PolygonSlice[]

#    # Preinitialize the array
#    for height in heights
#        push!(slices, PolygonSlice(LineSegment[],height))
#    end

#    for face in mesh.faces
#        zmin, zmax = extrema([face.vertices[j].e3 for j=1:3])
#        i = 1
#        for height in heights
#            if height > zmax
#                break
#            elseif zmin <= height
#                seg = LineSegment(face, height)
#                if !is(seg, nothing)
#                    push!(slices[i].segments, seg)
#                end
#            end
#            i = i + 1
#        end
#    end

#    return slices
#end


################################################################################
#
# Polygon
#
################################################################################

#Polygon() = Polygon(nil(LineSegment))

#function push!(poly::Polygon, f::LineSegment)
#    poly.segments = cons(f, poly.segments)
#end

#function Polygon(lines::Array{LineSegment})
#end


################################################################################
#
# PolygonMesh
#
################################################################################

PolygonMesh() = PolygonMesh(Bounds(Float64), Vector3{Float64}[], Face[])

function isapprox(a::Vector3, b::Vector3)
    return (isapprox(a[1],b[1])
           && isapprox(a[2],b[2])
           && isapprox(a[3],b[3]))
end

function push!(poly::PolygonMesh, verts::Array, normal::Vector3)
    vertidx = [0,0,0]
    for i = 1:length(poly.vertices), j = 1:3
        if isapprox(poly.vertices[i], verts[j])
            vertidx[j] = i
        end
    end

    for i = 1:3
        if vertidx[i] == 0
            push!(poly.vertices, verts[i])
            vertidx[i] = length(poly.vertices)
        end
    end

    update!(poly.bounds, [poly.vertices[i] for i in vertidx])
    push!(poly.faces, Face(Vector3(vertidx), normal))

end


function PolygonMesh(path::String)
    # create a mesh representation from an STL file location 
    file = open(path, "r")

    mesh = PolygonMesh()

    try
        # Discover file type
        if endswith(path, ".stl")
            header = ascii(readbytes(file, 5))
            if lowercase(header) == "solid"
                filetype = :ascii_stl
            else
                readbytes(file, 75) # throw out header
                filetype = :binary_stl
                read(file, Uint32) # throwout triangle count
            end
        end

        # Construct mesh
        mesh = PolygonMesh(file, filetype)

    finally
        close(file)
    end

    return mesh
end

function PolygonMesh(m::IOStream, filetype::Symbol)
    # Parses a file based on the given symbol and return a PolygonMesh
    # Symbol can be :ascii_stl, :binary_stl
    #  facet normal -1 0 0
    #    outer loop
    #      vertex 0 0 10
    #      vertex 0 10 10
    #      vertex 0 0 0
    #    endloop
    #  endfacet

    mesh = PolygonMesh()

    if filetype == :ascii_stl
        while !eof(m)
            line = split(lowercase(readline(m)))
            if line[1] == "facet"
                normal = Vector3(float64(line[3:5]))
                readline(m) # Throw away outerloop
                indices = [0,0,0]
                push!(mesh, [Vector3(float64(split(readline(m))[2:4])) for i=1:3], normal)
                readline(m) # throwout endloop
                readline(m) # throwout endfacet
            end
        end

    elseif filetype == :binary_stl
        while !eof(m)
            normal = Vector3([float64(read(m, Float32)) for i = 1:3])
            push!(mesh, [Vector3([float64(read(m, Float32)) for i = 1:3]) for j = 1:3], normal)
            read(m, Uint16) # throwout attribute
        end
    end

    # If we can't find anything, return nothing
    return mesh
end

function rotate!(mesh::PolygonMesh, angle::Float64, axis::Array{Float64}, through::Array{Float64})
    axis = axis/norm(axis) # normalize axis
    a, b, c = through
    u, v, w = axis
    for face in mesh.faces
        for i = 1:3
            x, y, z = face.vertices[i]
            face.vertices[i] = rotate(x, y, z, a, b, c, u, v, w, angle)
        end
        update!(mesh.bounds, face)
    end
end

function rotate(x, y, z, a, b, c, u, v, w, angle)
    # See: http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/#x1-10011
    return Vector3((a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(angle))+x*cos(angle)+(-c*v+b*w-w*y+v*z)*sin(angle),
            (b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(angle))+y*cos(angle)+(c*u-a*w+w*x-u*z)*sin(angle),
            (c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(angle))+z*cos(angle)+(-b*u+a*v-v*x+u*y)*sin(angle))

end

################################################################################
#
# Face
#
################################################################################

function (==)(a::Face, b::Face)
    return (a.vertidx == b.vertidx &&
            a.normal == b.normal)
end

################################################################################
#
# LineSegment
#
################################################################################

#function LineSegment(f::Face, z::Float64)

#    p0 = f.vertices[1]
#    p1 = f.vertices[2]
#    p2 = f.vertices[3]

#    if p0.e3 < z && p1.e3 >= z && p2.e3 >= z
#        return LineSegment(p0, p2, p1, z, f.normal)
#    elseif p0.e3 > z && p1.e3 < z && p2.e3 < z
#        return LineSegment(p0, p1, p2, z, f.normal)
#    elseif p1.e3 < z && p0.e3 >= z && p2.e3 >= z
#        return LineSegment(p1, p0, p2, z, f.normal)
#    elseif p1.e3 > z && p0.e3 < z && p2.e3 < z
#        return LineSegment(p1, p2, p0, z, f.normal)
#    elseif p2.e3 < z && p1.e3 >= z && p0.e3 >= z
#        return LineSegment(p2, p1, p0, z, f.normal)
#    elseif p2.e3 > z && p1.e3 < z && p0.e3 < z
#        return LineSegment(p2, p0, p1, z, f.normal)
#    else
#        return nothing
#    end

#end

#function LineSegment(p0::Vector3, p1::Vector3, p2::Vector3, z::Float64, normal::Vector3)
#    start = Vector2(p0.e1 + (p1.e1 - p0.e1) * (z - p0.e3) / (p1.e3 - p0.e3),
#                    p0.e2 + (p1.e2 - p0.e2) * (z - p0.e3) / (p1.e3 - p0.e3))
#    finish = Vector2(p0.e1 + (p2.e1 - p0.e1) * (z - p0.e3) / (p2.e3 - p0.e3),
#                     p0.e2 + (p2.e2 - p0.e2) * (z - p0.e3) / (p2.e3 - p0.e3))
#    return LineSegment(start, finish, normal);
#end


#function (==)(a::LineSegment, b::LineSegment)
#    return (a.start == b.start &&
#            a.finish == b.finish)
#end

################################################################################
#
# Bounds
#
################################################################################

function Bounds(T)
    return Bounds{T}(-Inf,-Inf,-Inf,Inf,Inf,Inf)
end

function update!(box::Bounds, verts::Array)
    # update the bounds against a face
    xmin, xmax = extrema([verts[i].e1 for i = 1:3])
    ymin, ymax = extrema([verts[i].e2 for i = 1:3])
    zmin, zmax = extrema([verts[i].e3 for i = 1:3])

    box.xmin = min(xmin, box.xmin)
    box.ymin = min(ymin, box.ymin)
    box.zmin = min(zmin, box.zmin)
    box.xmax = max(xmax, box.xmax)
    box.ymax = max(ymax, box.ymax)
    box.zmax = max(zmax, box.zmax)
end

function (==)(a::Bounds, b::Bounds)
    return (a.xmax == b.xmax &&
            a.ymax == b.ymax &&
            a.zmax == b.zmax &&
            a.xmin == b.xmin &&
            a.ymin == b.ymin &&
            a.zmin == b.zmin)
end


end # module
