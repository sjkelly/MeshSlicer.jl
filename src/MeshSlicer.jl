module MeshSlicer

type PolygonSlice
    segments
    layer
end

function slice(path::String, thickness)
    file = open(path, "r")

    mesh = PolygonMesh(file)

    close(file)

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

        initialSlice = convert(Int64, floor((face.vertices[face.order[1]][3] - mesh.bounds.zmin)/sliceheight))
        finalSlice = convert(Int64, floor((face.vertices[face.order[3]][3] - mesh.bounds.zmin)/sliceheight))

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

function findorder(vertices, index)
    # findorder
    # Given an array of vectors, return an ordered list of their maximum values.
    order = [1,1,1]
    for i = 1:3
        if vertices[i][index] < vertices[order[1]][index] #min
            order[2] = order[1]
            order[1] = i
        elseif vertices[i][index] > vertices[order[3]][index] #max
            order[2] = order[3]
            order[3] = i
        elseif vertices[i][index] == vertices[order[1]][index] #same as min
            order[2] = order[1]
        elseif vertices[i][index] == vertices[order[3]][index] #same as max
            order[2] = order[3]
        else
            order[2] = i
        end
    end
    return order
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

type PolygonMesh
    bounds
    faces
end

function PolygonMesh(m::IOStream)
    # create a mesh representation 
    
    faces = Face[]
    bounds = Bounds(0,0,0,Inf,Inf,Inf)
    
    while !eof(m)
        f = Face(m)
        if f != Nothing
            push!(faces, f)
            update!(bounds, f)
        end
    end
    
    return PolygonMesh(bounds, faces)
end


################################################################################
#
# Face:
#   vertices : [[x, y, z], ...]
#       An array of the vertices in the face
#   normal : [x::Float64, y::Float64, x::Float64]
#   order : [min, middle, max]
#       The order of the faces based on their ordinance against the slicng plane
#   angle : Float64
#       Angle in radians with the slicing plane. Ex: Normal of [0,0,1] = 0,
#       Normal of [1,0,0] = pi/2, Normal of [0,0,-1] = pi
#
#
# outer constructors:
#   Face(f::IOStream)
#       Pulls a face from an STL file IOStream
#
################################################################################

type Face
    vertices
    normal
    order # indices of min, middle, max z height
    angle # angle (in radians) the face makes with the slice plane
end

function Face(m::IOStream)
    #  facet normal -1 0 0
    #    outer loop
    #      vertex 0 0 10
    #      vertex 0 10 10
    #      vertex 0 0 0
    #    endloop
    #  endfacet
    vertices = Array[]
    normal = zeros(3)
    line = split(lowercase(readline(m)))
    if line[1] == "facet"
        normal = float64(line[3:5])
        normal = normal/norm(normal) # make sure normal is actually normal
        angle = acos(dot(normal, [0,0,1])) # find angle against plane
        readline(m) # Throw away outerloop
        for i = 1:3 # Get vertices
            line = split(lowercase(readline(m)))
            push!(vertices, float64(line[2:4]))
        end
        # Find Height ordering
        order = findorder(vertices, 3)

        return Face(vertices, normal, order, angle)
    else
        return Nothing
    end
end

function (==)(a::Face, b::Face)
    return (a.vertices == b.vertices &&
            a.normal == b.normal &&
            a.order == b.order &&
            a.angle == b.angle)
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

type LineSegment
    start
    finish
    slope

    LineSegment(x,y) = new(x,y,(x[1]-y[1])/(x[2]-y[2])) # compute slope
end

function LineSegment(f::Face, z)

    p0 = f.vertices[1]
    p1 = f.vertices[2]
    p2 = f.vertices[3]

    if p0[3] < z && p1[3] >= z && p2[3] >= z
        return LineSegment(p0, p2, p1, z)
    elseif p0[3] > z && p1[3] < z && p2[3] < z
        return LineSegment(p0, p1, p2, z)
    elseif p1[3] < z && p0[3] >= z && p2[3] >= z
        return LineSegment(p1, p0, p2, z)
    elseif p1[3] > z && p0[3] < z && p2[3] < z
        return LineSegment(p1, p2, p0, z)
    elseif p2[3] < z && p1[3] >= z && p0[3] >= z
        return LineSegment(p2, p1, p0, z)
    elseif p2[3] > z && p1[3] < z && p0[3] < z
        return LineSegment(p2, p0, p1, z)
    else
        return Nothing
    end

end

function LineSegment(p0::Array, p1::Array, p2::Array, z)
    start = zeros(2)
    finish = zeros(2)
    start[1] = p0[1] + (p1[1] - p0[1]) * (z - p0[3]) / (p1[3] - p0[3]);
    start[2] = p0[2] + (p1[2] - p0[2]) * (z - p0[3]) / (p1[3] - p0[3]);
    finish[1] = p0[1] + (p2[1] - p0[1]) * (z - p0[3]) / (p2[3] - p0[3]);
    finish[2] = p0[2] + (p2[2] - p0[2]) * (z - p0[3]) / (p2[3] - p0[3]);
    return LineSegment(start, finish);
end


function (==)(a::LineSegment, b::LineSegment)
    return (a.start == b.start &&
            a.finish == b.finish &&
            a.slope == b.slope)
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

type Bounds
    xmax
    ymax
    zmax
    xmin
    ymin
    zmin
end

function update!(box::Bounds, face::Face)
    xordering = findorder(face.vertices, 1)
    yordering = findorder(face.vertices, 2)
    zordering = face.order

    box.xmin = min(face.vertices[xordering[1]][1], box.xmin)
    box.ymin = min(face.vertices[yordering[1]][2], box.ymin)
    box.zmin = min(face.vertices[zordering[1]][3], box.zmin)
    box.xmax = max(face.vertices[xordering[3]][1], box.xmax)
    box.ymax = max(face.vertices[yordering[3]][2], box.ymax)
    box.zmax = max(face.vertices[zordering[3]][3], box.zmax)
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
