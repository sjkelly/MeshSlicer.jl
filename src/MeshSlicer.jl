# [https://github.com/sjkelly/MeshSlicer.jl](https://github.com/sjkelly/MeshSlicer.jl)
# ![](../img/sliced_cylinder.png)

module MeshSlicer

using ImmutableArrays
import Base.push!
import Base.show

export Bounds3, Bounds2, Face, PolygonMesh, LineSegment,
       update!, rotate!, rotate, MeshSlice, Polygon

type Bounds3{T<:Number}
    xmax::T
    ymax::T
    zmax::T
    xmin::T
    ymin::T
    zmin::T
end

type Bounds2{T<:Number}
    xmax::T
    ymax::T
    xmin::T
    ymin::T
end

type Face
    vertices::Array{Vector3{Float64}}
    normal::Vector3{Float64}
    zmin::Float64
    zmax::Float64
end

function Face(vertices::Array, normal::Vector3{Float64})
    zmax = max(max(vertices[1].e3, vertices[2].e3), vertices[3].e3)
    zmin = min(min(vertices[1].e3, vertices[2].e3), vertices[3].e3)
    return Face(vertices, normal, zmin, zmax)
end


type PolygonMesh
    bounds::Bounds3
    faces::Array{Face}
    patios::Array{Float64}
end

type LineSegment
    start::Vector2{Float64}
    finish::Vector2{Float64}
    normal::Vector3{Float64}
end

type Polygon
    bounds::Bounds2
    segments::Array{LineSegment}
end

type MeshSlice
    bounds::Bounds2
    polygons::Array{Polygon}
    layer::Float64
end


function show(io::IO, poly::Polygon)
    println("Polygon:")
    println("\tBounds:", poly.bounds)
    println("\tSegments:")
    for seg in poly.segments
        println("\t\tStart:", seg.start, " Finish:", seg.finish)
    end
end

# ##MeshSlice(*mesh::PolygonMesh, heights::Array{Float64}*)
#
# Create an array of MeshSlice at heights given by a
# monotonically increasing array of `Float64`.
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------------- --------------------------------------------------
# `mesh`               The `PolygonMesh` to slice.
#
# `heights`            An array of heights at which to slice.
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Returns:
# -------------------- --------------------------------------------------
# `Array{MeshSlice}`   An `Array` of `MeshSlice` at the requested height.
# ----------------------------------------------------------------------------
function MeshSlice(mesh::PolygonMesh, heights::Array{Float64}; eps=0.00001, autoeps=true)
    slices = [LineSegment[] for i = 1:length(heights)]
    bounds = [Bounds2() for i = 1:length(heights)]

    for face in mesh.faces
        i = 1
        for height in heights
            if height > face.zmax
                break
            elseif face.zmin <= height
                seg = LineSegment(face, height)
                if seg != nothing
                    push!(slices[i], seg)
                    update!(bounds[i], seg)
                end
            end
            i = i + 1
        end
    end

    polys = MeshSlice[]

    for i = 1:length(heights)
        push!(polys, MeshSlice(bounds[i], Polygon(slices[i], eps=eps, autoeps=autoeps), heights[i]))
    end

    return polys
end

# ##Polygon()
#
# ----------------------------------------------------------------------------
# Returns:
# -------------------- --------------------------------------------------
# `Polygon`            An empty `Polygon`.
# ----------------------------------------------------------------------------
Polygon() = Polygon(Bounds2(), LineSegment[])

# ##Polygon(*lines::Array{LineSegment}, eps::Real*)
#
# Construct an `Array` of `Polygon` from an `Array` of `LineSegments`.
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------------- --------------------------------------------------
# `lines`              An unorder `Array` of `LineSegment`.
#
# `eps`                The tolerance in which to count starting and finishing points as identical.
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Returns:
# -------------------- --------------------------------------------------
# `Array{Polygon}`     An `Array` of `Polygon`.
# ----------------------------------------------------------------------------
function Polygon(lines::Array{LineSegment}; eps=0.00001, autoeps=true)
    n = length(lines)
    if n == 0
        return [Polygon()]
    end
    polys = Polygon[]
    paired = [false for i = 1:n]
    start = 1
    seg = 1
    paired[seg] = true

    if autoeps
        for segment in lines
            eps = min(eps, norm(segment.start-segment.finish)/2)
        end
    end

    while true
        #Start new polygon with seg
        poly = Polygon()
        push!(poly, lines[seg])

        #Pair lines until we get to start point
        lastseg = seg
        while norm(lines[start].start - lines[seg].finish) >= eps
            lastseg = seg

            for i = 1:n
                if !paired[i]
                    if norm(lines[seg].finish - lines[i].start) <= eps
                        push!(poly, lines[i])
                        paired[i] = true
                        seg = i
                    end
                end
            end

            if (seg == start #We couldn't pair the segment
                || seg == lastseg) #The polygon can't be closed
                break
            end
        end

        if length(poly.segments) > 2
            closed = true
            if poly.segments[1].start != poly.segments[end].finish
                closed = false
            end
            for i = 1:length(poly.segments)-2
                if closed
                    break
                end
                for j = i+2:length(poly.segments)
                    if poly.segments[i].start == poly.segments[j].finish
                        poly.segments = poly.segments[i:j]
                        closed = true
                        break
                    end
                end
            end
            push!(polys,poly)
        end
        #start new polygon
        for i = 1:length(lines)
            if !paired[i] #Find next unpaired seg
                start = i
                paired[i] = true
                seg = start
                break
            elseif i == length(lines) #we have paired each segment
                return polys
            end
        end
    end
end

# ##PolygonMesh()
#
# ----------------------------------------------------------------------------
# Returns:
# -------------------- --------------------------------------------------
# `PolygonMesh`        An empty `PolygonMesh`.
# ----------------------------------------------------------------------------
PolygonMesh() = PolygonMesh(Bounds3(), Face[], Float64[])

function binarySTLvertex(file)
    return Vector3(float64(read(file, Float32)),
                   float64(read(file, Float32)),
                   float64(read(file, Float32)))
end

function binarySTLvertices(file)
    return Vector3[binarySTLvertex(file),
            binarySTLvertex(file),
            binarySTLvertex(file)]
end


# ##PolygonMesh(*path::String*)
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------------- --------------------------------------------------
# `mesh`               The `PolygonMesh` to slice.
#
# `height`             The height at which to slice.
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Returns:
# -------------------- --------------------------------------------------
# `MeshSlice`          A MeshSlice at the requested height.
# ----------------------------------------------------------------------------
function PolygonMesh(path::String)

    file = open(path, "r")

    mesh = PolygonMesh()

    try
        #Discover file type and construct mesh
        #STL
        if endswith(path, ".stl")
            header = ascii(readbytes(file, 5))

            #ASCII STL
            #https://en.wikipedia.org/wiki/STL_%28file_format%29#ASCII_STL
            if lowercase(header) == "solid"
                while !eof(file)
                    line = split(lowercase(readline(file)))
                    if line[1] == "facet"
                        normal = Vector3(float64(line[3:5]))
                        readline(file) # Throw away outerloop
                        vertices = [Vector3(float64(split(readline(file))[2:4])) for i = 1:3]
                        readline(file) # throwout endloop
                        readline(file) # throwout endfacet
                        push!(mesh, Face(vertices,normal))
                    end
                end

            #Binary STL
            #https://en.wikipedia.org/wiki/STL_%28file_format%29#Binary_STL
            else
                readbytes(file, 75) # throw out header
                read(file, Uint32) # throwout triangle count
                while !eof(file)
                    normal = binarySTLvertex(file)
                    vertices = binarySTLvertices(file)
                    skip(file, 2) # throwout 16bit attribute
                    push!(mesh, Face(vertices,normal))
                end
            end
        end

    finally
        close(file)
    end

    return mesh
end

# ##rotate!(*mesh::PolygonMesh, angle::Float64, axis::Array{Float64}*, through::Array{Float64})
#
# Rotates a `PolygonMesh` around an arbitrary axis.
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------------- --------------------------------------------------
# `mesh`               The `PolygonMesh` to rotate.
#
# `angle`              The rotation angle in radians.
#
# `axis`               The axis  to rotate around.
#
# `through`            The `[x, y, z]` point the `axis` passes through.
# ----------------------------------------------------------------------------
function rotate!(mesh::PolygonMesh, angle::Float64, axis::Array{Float64}, through::Array{Float64})
    axis = axis/norm(axis) # normalize axis
    a, b, c = through
    u, v, w = axis
    for face in mesh.faces
        face.vertices = [begin
                            x, y, z = face.vertices[i];
                            rotate(x, y, z, a, b, c, u, v, w, angle)
                         end for i = 1:3]
        face.zmax = max(max(face.vertices[1].e3, face.vertices[2].e3), face.vertices[3].e3)
        face.zmin = min(min(face.vertices[1].e3, face.vertices[2].e3), face.vertices[3].e3)
        update!(mesh.bounds, face)
    end
end

# ##rotate(*x, y, z, a, b, c, u, v, w, angle*)
#
# Create a MeshSlice from a PolygonMesh at a given height.
# See: [http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/#x1-10011](http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/#x1-10011)
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------------- --------------------------------------------------
# `x, y, z`            The position to rotate.
#
# `u, v, w`            The axis to rotate around.
#
# `a, b, c`            The point to rotate through.
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Returns:
# -------------------- --------------------------------------------------
# `Vector3`            The rotation of `x, y, z`.
# ----------------------------------------------------------------------------
function rotate(x, y, z, a, b, c, u, v, w, angle)
    return Vector3((a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(angle))+x*cos(angle)+(-c*v+b*w-w*y+v*z)*sin(angle),
            (b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(angle))+y*cos(angle)+(c*u-a*w+w*x-u*z)*sin(angle),
            (c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(angle))+z*cos(angle)+(-b*u+a*v-v*x+u*y)*sin(angle))

end

# ##push!(*mesh::PolygonMesh, f::Face*)
#
# Add a `Face` to a `PolygonMesh` object.
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------------- --------------------------------------------------
# `mesh`               The `PolygonMesh` to mutate.
#
# `f`                  The `Face` to add to mesh.
# ----------------------------------------------------------------------------
function push!(mesh::PolygonMesh, f::Face)
    if abs(f.normal.e3) == 1 && !in(f.vertices[1].e3, mesh.patios)
        push!(mesh.patios,f.vertices[1].e3)
        sort!(mesh.patios) #Keep monotonic
    end
    push!(mesh.faces, f)
    update!(mesh.bounds, f)
end

function push!(poly::Polygon, seg::LineSegment)
    push!(poly.segments, seg)
    update!(poly.bounds, seg)
end

# ##LineSegment(*f::Face, z::Float64*)
#
# Get a `LineSegment` in the X-Y plane from a `Face` at a requested height.
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------------- --------------------------------------------------
# `f`                  The `Face` to slice.
#
# `z`                  The height to slice the face.
# ----------------------------------------------------------------------------
# Returns:
# -------------------- --------------------------------------------------
# `LineSegment`        A `LineSegment` at the requested height in the slice plane.
# ----------------------------------------------------------------------------
function LineSegment(f::Face, z::Float64)
    p0 = f.vertices[1]
    p1 = f.vertices[2]
    p2 = f.vertices[3]

    if p0.e3 < z && p1.e3 >= z && p2.e3 >= z
        return LineSegment(p0, p2, p1, z, f.normal)
    elseif p0.e3 > z && p1.e3 < z && p2.e3 < z
        return LineSegment(p0, p1, p2, z, f.normal)
    elseif p1.e3 < z && p0.e3 >= z && p2.e3 >= z
        return LineSegment(p1, p0, p2, z, f.normal)
    elseif p1.e3 > z && p0.e3 < z && p2.e3 < z
        return LineSegment(p1, p2, p0, z, f.normal)
    elseif p2.e3 < z && p1.e3 >= z && p0.e3 >= z
        return LineSegment(p2, p1, p0, z, f.normal)
    elseif p2.e3 > z && p1.e3 < z && p0.e3 < z
        return LineSegment(p2, p0, p1, z, f.normal)
    end

end

# ##LineSegment(*p0::Vector3, p1::Vector3, p2::Vector3, z::Float64, normal::Vector3*)
#
# Get a `LineSegment` at `z` generated by the rays $\overline{p0p1}$ and $\overline{p0p2}$.
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------------- --------------------------------------------------
# `p0`                 The base point.
#
# `p1`                 A point such that `z` is between `p0[3]` and `p1[3]`.
#
# `p2`                 A point such that `z` is between `p0[3]` and `p2[3]`.
#
# `z`                  The height to slice the rays $\overline{p0p1}$ and $\overline{p0p2}$.
#
# `normal`             The normal of the face, retained for insetting and offsetting.
# ----------------------------------------------------------------------------
# Returns:
# -------------------- --------------------------------------------------
# `LineSegment`        A `LineSegment` at the requested height in the slice plane.
# ----------------------------------------------------------------------------
function LineSegment(p0::Vector3, p1::Vector3, p2::Vector3, z::Float64, normal::Vector3)
    start = Vector2(p0.e1 + (p1.e1 - p0.e1) * (z - p0.e3) / (p1.e3 - p0.e3),
                    p0.e2 + (p1.e2 - p0.e2) * (z - p0.e3) / (p1.e3 - p0.e3))
    finish = Vector2(p0.e1 + (p2.e1 - p0.e1) * (z - p0.e3) / (p2.e3 - p0.e3),
                     p0.e2 + (p2.e2 - p0.e2) * (z - p0.e3) / (p2.e3 - p0.e3))
    return LineSegment(start, finish, normal);
end

# ##==(*a::LineSegment, b::LineSegment*)
#
# ----------------------------------------------------------------------------
# Returns:
# -------------------- --------------------------------------------------
# `Bool`               Whether two `LineSegment`s are equal.
# ----------------------------------------------------------------------------
function ==(a::LineSegment, b::LineSegment)
    return (a.start == b.start &&
            a.finish == b.finish)
end

# ##Bounds3()
#
# ----------------------------------------------------------------------------
# Returns:
# -------------------- --------------------------------------------------
# `Bound3`              An empty `Bounds3` of type `Float64`.
# ----------------------------------------------------------------------------
Bounds3() = Bounds3{Float64}(-Inf,-Inf,-Inf,Inf,Inf,Inf)
Bounds2() = Bounds2{Float64}(-Inf,-Inf,Inf,Inf)
# ##update!(*box::Bounds, face::Face*)
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------------- --------------------------------------------------
# `box`                The `Bounds` to update.
#
# `face`               The `Face` to update the bounds against.
# ----------------------------------------------------------------------------
function update!(box::Bounds3, face::Face)
    #update the bounds against a face
    x1, y1, z1 = face.vertices[1]
    x2, y2, z2 = face.vertices[2]
    x3, y3, z3 = face.vertices[3]

    box.xmin = min(min(min(box.xmin, x1), x2), x3)
    box.ymin = min(min(min(box.ymin, y1), y2), y3)
    box.zmin = min(box.zmin, face.zmin)
    box.xmax = max(max(max(box.xmax, x1), x2), x3)
    box.ymax = max(max(max(box.ymax, y1), y2), y3)
    box.zmax = max(box.zmax, face.zmax)
end

function update!(box::Bounds2, line::LineSegment)
    #update the bounds against a line segement
    box.xmin = min(box.xmin, min(line.start[1], line.finish[1]))
    box.ymin = min(box.ymin, min(line.start[2], line.finish[2]))
    box.xmax = max(box.xmax, max(line.start[1], line.finish[1]))
    box.ymax = max(box.xmax, max(line.start[2], line.finish[2]))
end

function update!(b1::Bounds2, b2::Bounds2)
    b1.xmin = min(b1.xmin, b2.xmin)
    b1.ymin = min(b1.ymin, b2.ymin)
    b1.xmax = max(b1.xmax, b2.xmax)
    b1.ymax = max(b1.ymax, b2.ymax)
end

# ##==(*a::Bounds3, b::Bounds3*)
#
# ----------------------------------------------------------------------------
# Returns:
# -------------------- --------------------------------------------------
# `Bool`               Whether two `Bounds3` are equal.
# ----------------------------------------------------------------------------
function ==(a::Bounds3, b::Bounds3)
    return (a.xmax == b.xmax &&
            a.ymax == b.ymax &&
            a.zmax == b.zmax &&
            a.xmin == b.xmin &&
            a.ymin == b.ymin &&
            a.zmin == b.zmin)
end

function ==(a::Bounds2, b::Bounds2)
    return (a.xmax == b.xmax &&
            a.ymax == b.ymax &&
            a.xmin == b.xmin &&
            a.ymin == b.ymin)
end

end # module
