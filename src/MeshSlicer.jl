module MeshSlicer

type Bounds
    xmax::Float64
    ymax::Float64
    zmax::Float64
    xmin::Float64
    ymin::Float64
    zmin::Float64
end

type Face
    vertices::Array{Array}
    normal::Array{Float64}

    Face(v, n) = new(sort!(v, by=x->x[3]), n/norm(n))
end

type PolygonMesh
    bounds::Bounds
    faces::Array{Face}
end

type LineSegment
    start::Array{Float64}
    finish::Array{Float64}
    normal::Array{Float64}
end

type PolygonSlice
    segments::Array{LineSegment}
    layer::Float64
end


function slice(path::String, thickness::Float64)

    mesh = PolygonMesh(path)

    startZ = mesh.bounds.zmin

    #We can only print an integer number of layers
    layercount = round((mesh.bounds.zmax - mesh.bounds.zmin)/thickness)

    #Adjust sliceheight
    sliceheight = (mesh.bounds.zmax - mesh.bounds.zmin)/layercount

    layers = [mesh.bounds.zmin:sliceheight:mesh.bounds.zmax]

    segmentlist = Array(PolygonSlice,convert(Int64,layercount))

    for i = 1:layercount
        segmentlist[i] = PolygonSlice(LineSegment[],layers[i])
    end

    #println(segmentlist)
    for face in mesh.faces

        initialSlice = convert(Int64, floor((face.vertices[1][3] - mesh.bounds.zmin)/sliceheight))
        finalSlice = convert(Int64, floor((face.vertices[3][3] - mesh.bounds.zmin)/sliceheight))

        locallayer = layers[initialSlice+1:finalSlice]

        index = initialSlice + 1
        for layer in locallayer
            seg = LineSegment(face, layer)
            if seg != Nothing
                push!(segmentlist[index].segments, seg)
            end
            index = index + 1
        end
    end
    return segmentlist
end

################################################################################
#
# PolygonMesh:
#   bounds
#   faces
#
# outer constructors:
#   PolygonMesh(m::IOStream)
#       Create a mesh from an STL file IOStream
#
################################################################################

function PolygonMesh(path::String)
    # create a mesh representation 
    file = open(path, "r")

    faces = Face[]
    bounds = Bounds(0,0,0,Inf,Inf,Inf)

    # Discover file type
    if endswith(path, ".stl")
        header = ascii(readbytes(file, 5))
        if lowercase(header) == "solid"
            s = :ascii_stl
        else
            readbytes(file, 75) # throw out header
            s = :binary_stl
            read(file, Uint32) # throwout triangle count
        end
    end

    # Construct mesh
    while !eof(file)
        f = Face(file, s)
        if f != Nothing
            push!(faces, f)
            update!(bounds, f)
        end
    end

    close(file)
    return PolygonMesh(bounds, faces)
end


################################################################################
#
# Face:
#   vertices : [[x, y, z], ...]
#       An array of the vertices in the face
#   normal : [x::Float64, y::Float64, x::Float64]
#
#
# outer constructors:
#   Face(f::IOStream)
#       Pulls a face from an STL file IOStream
#
################################################################################

function Face(m::IOStream, s::Symbol)
    #  facet normal -1 0 0
    #    outer loop
    #      vertex 0 0 10
    #      vertex 0 10 10
    #      vertex 0 0 0
    #    endloop
    #  endfacet
    vertices = [zeros(3) for i = 1:3]
    normal = zeros(3)
    if s == :ascii_stl
        line = split(lowercase(readline(m)))
        if line[1] == "facet"
            normal = float64(line[3:5])
            readline(m) # Throw away outerloop
            for i = 1:3 # Get vertices
                line = split(lowercase(readline(m)))
                vertices[i] = float64(line[2:4])
            end
            return Face(vertices, normal)
        end

    elseif s == :binary_stl
        normal = [read(m, Float32) for i = 1:3]
        vertices = [[read(m, Float32) for i = 1:3] for j = 1:3]
        read(m, Uint16) # throwout attribute
        return Face(vertices, normal)
    end

    # If we can't find anything, return nothing
    return Nothing
end

function (==)(a::Face, b::Face)
    return (a.vertices == b.vertices &&
            a.normal == b.normal)
end

################################################################################
#
# LineSegment:
#   start : [x::Float64, y::Float64]
#   finish : [x::Float64, y::Float64]
#   slope :
#       slope in slice plane, computed automatically by inner constructor
#
# outer constructors:
#   LineSegment(f::Face, z::Number)
#   LineSegment(p0, p1, p2, z::Number)
#       p0, p1, p2 are expected to be Arrays of size 3 containing numbers
#
################################################################################

function LineSegment(f::Face, z::Float64)

    p0 = f.vertices[1]
    p1 = f.vertices[2]
    p2 = f.vertices[3]

    if p0[3] < z && p1[3] >= z && p2[3] >= z
        return LineSegment(p0, p2, p1, z, f.normal)
    elseif p0[3] > z && p1[3] < z && p2[3] < z
        return LineSegment(p0, p1, p2, z, f.normal)
    elseif p1[3] < z && p0[3] >= z && p2[3] >= z
        return LineSegment(p1, p0, p2, z, f.normal)
    elseif p1[3] > z && p0[3] < z && p2[3] < z
        return LineSegment(p1, p2, p0, z, f.normal)
    elseif p2[3] < z && p1[3] >= z && p0[3] >= z
        return LineSegment(p2, p1, p0, z, f.normal)
    elseif p2[3] > z && p1[3] < z && p0[3] < z
        return LineSegment(p2, p0, p1, z, f.normal)
    else
        return Nothing
    end

end

function LineSegment(p0::Array{Float64}, p1::Array{Float64}, p2::Array{Float64}, z::Float64, normal::Array{Float64})
    start = zeros(2)
    finish = zeros(2)
    start[1] = p0[1] + (p1[1] - p0[1]) * (z - p0[3]) / (p1[3] - p0[3]);
    start[2] = p0[2] + (p1[2] - p0[2]) * (z - p0[3]) / (p1[3] - p0[3]);
    finish[1] = p0[1] + (p2[1] - p0[1]) * (z - p0[3]) / (p2[3] - p0[3]);
    finish[2] = p0[2] + (p2[2] - p0[2]) * (z - p0[3]) / (p2[3] - p0[3]);
    return LineSegment(start, finish, normal);
end


function (==)(a::LineSegment, b::LineSegment)
    return (a.start == b.start &&
            a.finish == b.finish)
end

################################################################################
#
# Bounds:
#   xmax
#   ymax
#   zmax
#   xmin
#   ymin
#   zmin
#
################################################################################


function update!(box::Bounds, face::Face)
    x = sort(face.vertices, by=x->x[1]) # Sort by x
    y = sort(face.vertices, by=x->x[2]) # Sort by y

    box.xmin = min(x[1][1], box.xmin)
    box.ymin = min(y[1][2], box.ymin)
    box.zmin = min(face.vertices[1][3], box.zmin)
    box.xmax = max(x[3][1], box.xmax)
    box.ymax = max(y[3][2], box.ymax)
    box.zmax = max(face.vertices[3][3], box.zmax)
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
