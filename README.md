# MeshSlicer
This is a library for slicing [polygon mesh](http://en.wikipedia.org/wiki/Polygon_mesh) structures into flat polygons.

![](./img/sliced_cylinder.png)

## Install
This package is not yet in the Julia package repository. For now, you can call ```Pkg.clone(https://github.com/sjkelly/MeshSlicer.jl.git)``` in the Julia REPL.

## Supported File Formats
* STL (ASCII and Binary)


## API Reference

### Composite Types

#### Bounds
```
type Bounds{T<:Number}
    xmax::T
    ymax::T
    zmax::T
    xmin::T
    ymin::T
    zmin::T
end
```

**Constructors:**
* ```Bounds()```

**Operators:**
* ``` == ```


#### Face
```
type Face
    vertices::Array{Vector3{Float64}}
    normal::Vector3{Float64}
    next::Union(Face,Nothing)

    Face(v, n) = new(v, n, nothing)
end
```
**Constructor:**
* ```Face(m::IOStream, s::Symbol)```

**Operators:**
* ``` == ```


#### PolygonMesh
```
type PolygonMesh
    bounds::Bounds
    faces::Union(Face,Nothing)
end
```
**Constructors:**
* ```PolygonMesh() = PolygonMesh(Bounds(), nothing)```
* ```PolygonMesh(path::String)```

**Operators:**
* Iteration over faces. Ex: ```for face in mesh::PolygonMesh```

#### LineSegment
```
type LineSegment
    start::Vector2{Float64}
    finish::Vector2{Float64}
    normal::Vector3{Float64}
end
```
**Constructors:**
* ```LineSegment(f::Face, z::Float64)```
* ```LineSegment(p0::Vector3, p1::Vector3, p2::Vector3, z::Float64, normal::Vector3)```


**Operators:**
* ``` == ```

#### PolygonSlice
```
type PolygonSlice
    segments::Array{LineSegment}
    layer::Float64
end
```

**Constructors:**

* ```PolygonSlice(mesh::PolygonMesh, height::Float64)```
* ```PolygonSlice(mesh::PolygonMesh, heights::Array{Float64})```



### Functions
* update!(box::Bounds, face::Face)
* rotate!(mesh::PolygonMesh, angle::Float64, axis::Array{Float64}, through::Array{Float64})
* rotate(x, y, z, a, b, c, u, v, w, angle)
* push!(mesh::PolygonMesh, f::Face)


## Examples
* http://nbviewer.ipython.org/github/sjkelly/MeshSlicer.jl/blob/master/examples/slice_tree.ipynb

## Build Status
[![Build Status](https://travis-ci.org/sjkelly/MeshSlicer.jl.svg)](https://travis-ci.org/sjkelly/MeshSlicer.jl)
[![Coverage Status](https://img.shields.io/coveralls/sjkelly/MeshSlicer.jl.svg)](https://coveralls.io/r/sjkelly/MeshSlicer.jl)

This package is developed under the latest [development verion of Julia](https://github.com/julialang/julia).

## License
The MeshSlicer.jl package is licensed under the MIT "Expat" License. See [LICENSE.md](./LICENSE.md).
