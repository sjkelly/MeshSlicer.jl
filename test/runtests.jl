#! /usr/bin/env julia

using MeshSlicer
using Base.Test
using ImmutableArrays

# Make pwd MeshSlicer/test
data_path = Pkg.dir("MeshSlicer")*"/test/data/"

# getmesh/getface test
println("Testing PolygonMesh...")

# get the first face from cube.stl
testFace = Face(Vector3[Vector3(10.0,0.0,10.0),
                      Vector3(10.0,0.0,0.0),
                      Vector3(10.0,10.0,10.0)],
                      Vector3(1.0,0.0,0.0))

testBounds = Bounds3{Float64}(10.0,10.0,10.0,0,0,0)

resultmesh = PolygonMesh(data_path*"cube.stl")


@test resultmesh.faces[end].vertices == testFace.vertices
@test resultmesh.faces[end].normal == testFace.normal


println("Testing Bounds...")
@test resultmesh.bounds.xmax == testBounds.xmax
@test resultmesh.bounds.ymax == testBounds.ymax
@test resultmesh.bounds.zmax == testBounds.zmax
@test resultmesh.bounds.xmin == testBounds.xmin
@test resultmesh.bounds.ymin == testBounds.ymin
@test resultmesh.bounds.zmin == testBounds.zmin
@test resultmesh.bounds == testBounds


# Test bounds on translated cube
# translate([1,2,3])cube([4,5,6]);
testBounds2 = Bounds3(4.0,6.0,9.0,1.0,2.0,4.0)
updateFace = Face(Vector3[Vector3(5.0,3.0,4.0),
                      Vector3(2.0,7.0,4.0),
                      Vector3(2.0,3.0,3.0)],
                      Vector3(-1.0,0.0,0.0)) # contains correct bounds for cube
update!(testBounds2, updateFace)
resultmesh2 = PolygonMesh(data_path*"translated_cube.stl")

@test resultmesh2.bounds.xmax == testBounds2.xmax
@test resultmesh2.bounds.ymax == testBounds2.ymax
@test resultmesh2.bounds.zmax == testBounds2.zmax
@test resultmesh2.bounds.xmin == testBounds2.xmin
@test resultmesh2.bounds.ymin == testBounds2.ymin
@test resultmesh2.bounds.zmin == testBounds2.zmin
@test resultmesh2.bounds == testBounds2


# Test patio detection
println("Tesing patio detection...")

@test resultmesh2.patios == [3.0, 9.0]
@test resultmesh.patios == [0.0, 10.0]


# Test LineSegment
testpoint1 = Vector3[Vector3(0.0,0.0,0.0),Vector3(0.0,1.0,0.0),Vector3(0.0,1.0,1.0)]
normal = Vector3(1.0, 1.0, 1.0)
l1 = LineSegment(testpoint1[3],testpoint1[1],testpoint1[2], 0.5, normal)
l2 = LineSegment(testpoint1[3],testpoint1[2],testpoint1[1], 0.5, normal)
l3 = LineSegment(testpoint1[3],testpoint1[2],testpoint1[1], 0.0, normal)

testpoint2 = Vector3[Vector3(0.0,0.0,0.0),Vector3(1.0,1.0,1.0),Vector3(0.0,1.0,0.0)]
l4 = LineSegment(testpoint2[2],testpoint2[1],testpoint2[3], 0.5, normal)
l5 = LineSegment(testpoint2[2],testpoint2[3],testpoint2[1], 0.5, normal)

l6 = LineSegment(testFace, 5.0)
l7 = LineSegment(testpoint2[2],testpoint2[3],testpoint2[1], 1.0, normal)

testpoint3 = Vector3[Vector3(0.0,0.0,0.0),Vector3(0.0,1.0,1.0),Vector3(0.0,0.0,2.0)]
testface3 = Face(testpoint3, normal)
l8 = LineSegment(testface3, 0.9)
l9 = LineSegment(Face(Vector3[testpoint3[3], testpoint3[1], testpoint3[2]], normal), 1.5)

@test l1.start == Vector2(0.0, 0.5)
@test l1.finish == Vector2(0.0, 1.0)
@test l2.start == Vector2(0.0, 1.0)
@test l2.finish == Vector2(0.0, 0.5)
@test l3.start == Vector2(0.0, 1.0)
@test l3.finish == Vector2(0.0, 0.0)
@test l4.start == Vector2(0.5, 0.5)
@test l4.finish == Vector2(0.5, 1.0)
@test l5.start == Vector2(0.5, 1.0)
@test l6.start == Vector2(10.0, 0.0)
@test l6.finish == Vector2(10.0, 5.0)
@test l6.normal == Vector3(1.0,0.0,0.0)
@test l7.start == Vector2(1.0, 1.0)
@test l7.finish == Vector2(1.0, 1.0)
@test l8.start == Vector2(0.0, 0.0)
@test l8.finish == Vector2(0.0, 0.9)
@test l9.start == Vector2(0.0, 0.0)
@test l9.finish == Vector2(0.0, 0.5)


# test binary STL
println("Testing Binary STL import...")
binarymesh = PolygonMesh(data_path*"cube_binary.stl")
@test binarymesh.bounds == Bounds3(10.0,10.0,10.0,0.0,0.0,0.0)

println("Testing slicing...")
slice1 = MeshSlice(binarymesh, [0.0,5.0,10.0])
for slice in slice1
    for poly in slice.polygons
        n = length(poly.segments)
        for i = 1:n-1
            @test poly.segments[i].finish == poly.segments[i+1].start
        end
    end
end

println("Tesing Polygon construction...")
ep = 0.1
l1 = LineSegment(Vector2(0.0+ep,0.0), Vector2(0.0,1.0+ep), Vector3(0.0,0.0,0.0))
l2 = LineSegment(Vector2(0.0+ep,1.0), Vector2(0.5-ep,1.0), Vector3(0.0,0.0,0.0))
l3 = LineSegment(Vector2(0.5+ep,1.0), Vector2(0.0,0.0+ep), Vector3(0.0,0.0,0.0))
@test length(Polygon([l1,l2,l3], eps=ep*2, autoeps=false)) == 1
@test length(Polygon([l1,l2,l3], eps=ep/2, autoeps=false)) == 0

# test mesh rotation
println("Testing Mesh Rotation...")
r1 = rotate(1, 1, 0, 0, 0, 0, 1, 0, 0, pi/2)
r2 = rotate(1, 1, 0, 0, 0, 0, 0, 1, 0, pi/2)
r3 = rotate(1, 1, 0, 0, 0, 0, 0, 0, 1, pi/2)
rr1 = [1.0, 0, 1.0]
rr2 = [0.0, 1.0, -1.0]
rr3 = [-1.0, 1.0, 0.0]
for i = 1:3
    @test isapprox(r1[i] ,rr1[i])
    @test isapprox(r2[i] ,rr2[i])
    @test isapprox(r3[i] ,rr3[i])
end

rotate!(binarymesh, 45.0, [1.0,0.0,0.0], [0.0,0.0,0.0])
@test binarymesh.bounds == Bounds3(10.0,10.0,13.762255133518483,0.0,-8.509035245341185,0.0)

rotate!(binarymesh, 45.0, [0.0,1.0,0.0], [0.0,0.0,0.0])
@test binarymesh.bounds == Bounds3(16.96357128682594,10.0,13.762255133518483,0.0,-8.509035245341185,-8.509035245341185)

# Issue 21
mesh = PolygonMesh(data_path*"cube.stl")
a = MeshSlice(mesh, [eps()])
@test length(a) == 1
@test length(a[1].polygons) == 1
@test length(a[1].polygons[1].segments) == 4

