module MeshSlicer

# Setup types
type Face
    vertices
    normal
    heights # indices of min, middle, max z height
end

type Bounds
    xmax
    ymax
    zmax
    xmin
    ymin
    zmin
end

type Volume
    bounds
    faces
end

type UnstitchedPolygon
    segments
    layer
end


# End Type Setup


function slice(path::String, thickness)
    file = open(path, "r")

    volume = Volume(file)

    startZ = volume.bounds.zmin

    #We can only print an integer number of layers
    layercount = round((volume.bounds.zmax - volume.bounds.zmin)/thickness)

    #Adjust sliceheight
    sliceheight = (volume.bounds.zmax - volume.bounds.zmin)/layercount

    layers = [volume.bounds.zmin:sliceheight:volume.bounds.zmax]

    segmentlist = Array(UnstitchedPolygon,convert(Int64,layercount))

    for i = 1:layercount
        segmentlist[i] = UnstitchedPolygon(LineSegment[],layers[i])
    end

    #println(segmentlist)
    for face in volume.faces

        initialSlice = convert(Int64, floor((face.vertices[face.heights[1]][3] - volume.bounds.zmin)/sliceheight))
        finalSlice = convert(Int64, floor((face.vertices[face.heights[3]][3] - volume.bounds.zmin)/sliceheight))

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
    return (segmentlist)
end

function Volume(m::IOStream)
    # create a volume representation 
    
    faces = Face[]
    bounds = Bounds(0,0,0,Inf,Inf,Inf)
    
    while !eof(m)
        f = Face(m)
        if f != Nothing
            push!(faces, f)
            updatebounds!(bounds, f)
        end
    end
    
    return Volume(bounds, faces)
end

function updatebounds!(box::Bounds, face)
    xordering = findorder(face.vertices, 1)
    yordering = findorder(face.vertices, 2)
    zordering = face.heights

    box.xmin = min(face.vertices[xordering[1]][1], box.xmin)
    box.ymin = min(face.vertices[yordering[1]][2], box.ymin)
    box.zmin = min(face.vertices[zordering[1]][3], box.zmin)
    box.xmax = max(face.vertices[xordering[3]][1], box.xmax)
    box.ymax = max(face.vertices[yordering[3]][2], box.ymax)
    box.zmax = max(face.vertices[zordering[3]][3], box.zmax)
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
        readline(m) # Throw away outerloop
        for i = 1:3 # Get vertices
            line = split(lowercase(readline(m)))
            push!(vertices, float64(line[2:4]))
        end
        # Find Height ordering
        heights = findorder(vertices, 3)

        return Face(vertices, normal, heights)
    else
        return Nothing
    end
end

function findorder(vertices, index)
    # findheights
    # Given an array of vectors, return an ordered list of their maximum values.
    heights = [1,1,1]
    for i = 1:3
        if vertices[i][index] < vertices[heights[1]][index] #min
            heights[2] = heights[1]
            heights[1] = i
        elseif vertices[i][index] > vertices[heights[3]][index] #max
            heights[2] = heights[3]
            heights[3] = i
        elseif vertices[i][index] == vertices[heights[1]][index] #same as min
            heights[2] = heights[1]
        elseif vertices[i][index] == vertices[heights[3]][index] #same as max
            heights[2] = heights[3]
        else
            heights[2] = i
        end
    end
    return heights
end

################################################################################
#
# LineSegment:
#   start : [x::Float64, y::Float64]
#   finish : [x::Float64, y::Float64]
#   layer : z::Number
#
# constructors:
#   LineSegment(f::Face, z::Number)
#   LineSegment(p0, p1, p2, z::Number)
#       p0, p1, p2 are expected to be Arrays of size 3 containing numbers
#
################################################################################


type LineSegment
    start
    finish
    layer
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
    return LineSegment(start, finish, z);
end

end # module
