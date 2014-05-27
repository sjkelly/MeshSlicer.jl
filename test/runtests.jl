using MeshSlicer
using Base.Test
using ImmutableArrays

# Since we use relative paths, you need to test from MeshSlicer/test

# setup test handlers
test_handler(r::Test.Success) = println("Test success: $(r.expr)")
test_handler(r::Test.Failure) = error("Test fail: $(r.expr)")
test_handler(r::Test.Error) = rethrow(r)

# getmesh/getface test
println("Testing PolygonMesh...")

# get the first face from cube.stl
testFace = Face(Vector3[Vector3(0.0,0.0,0.0),
                      Vector3(0.0,0.0,10.0),
                      Vector3(0.0,10.0,10.0)],
                      Vector3(-1.0,0.0,0.0))

testBounds = Bounds(10.0,10.0,10.0,0,0,0)

resultmesh = PolygonMesh("./data/cube.stl")

Test.with_handler(test_handler) do
    @test resultmesh.faces[1].vertices == testFace.vertices
    @test resultmesh.faces[1].normal == testFace.normal
    @test resultmesh.faces[1] == testFace
end

println("Testing Bounds...")
Test.with_handler(test_handler) do
    @test resultmesh.bounds.xmax == testBounds.xmax
    @test resultmesh.bounds.ymax == testBounds.ymax
    @test resultmesh.bounds.zmax == testBounds.zmax
    @test resultmesh.bounds.xmin == testBounds.xmin
    @test resultmesh.bounds.ymin == testBounds.ymin
    @test resultmesh.bounds.zmin == testBounds.zmin
    @test resultmesh.bounds == testBounds
end


# Test bounds on translated cube
# translate([1,2,3])cube([4,5,6]);
testBounds2 = Bounds(4.0,6.0,9.0,1.0,2.0,4.0)
updateFace = Face(Vector3[[5.0,3.0,4.0],
                      [2.0,7.0,4.0],
                      [2.0,3.0,3.0]],
                      Vector3(-1.0,0.0,0.0)) # contains correct bounds for cube
update!(testBounds2, updateFace)
resultmesh2 = PolygonMesh("./data/translated_cube.stl")
Test.with_handler(test_handler) do
    @test resultmesh2.bounds.xmax == testBounds2.xmax
    @test resultmesh2.bounds.ymax == testBounds2.ymax
    @test resultmesh2.bounds.zmax == testBounds2.zmax
    @test resultmesh2.bounds.xmin == testBounds2.xmin
    @test resultmesh2.bounds.ymin == testBounds2.ymin
    @test resultmesh2.bounds.zmin == testBounds2.zmin
    @test resultmesh2.bounds == testBounds2
end


# Test LineSegment
testpoint1 = Vector3[[0.0,0.0,0.0],[0.0,1.0,0.0],[0.0,1.0,1.0]]
normal = Vector3(1.0, 1.0, 1.0)
l1 = LineSegment(testpoint1[3],testpoint1[1],testpoint1[2], 0.5, normal)
l2 = LineSegment(testpoint1[3],testpoint1[2],testpoint1[1], 0.5, normal)
l3 = LineSegment(testpoint1[3],testpoint1[2],testpoint1[1], 0.0, normal)

testpoint2 = Vector3[[0.0,0.0,0.0],[1.0,1.0,1.0],[0.0,1.0,0.0]]
l4 = LineSegment(testpoint2[2],testpoint2[1],testpoint2[3], 0.5, normal)
l5 = LineSegment(testpoint2[2],testpoint2[3],testpoint2[1], 0.5, normal)

Test.with_handler(test_handler) do
    @test l1.start == Vector2(0.0, 0.5)
    @test l1.finish == Vector2(0.0, 1.0)
    @test l2.start == Vector2(0.0, 1.0)
    @test l2.finish == Vector2(0.0, 0.5)
    @test l3.start == Vector2(0.0, 1.0)
    @test l3.finish == Vector2(0.0, 0.0)
    @test l4.start == Vector2(0.5, 0.5)
    @test l4.finish == Vector2(0.5, 1.0)
    @test l5.start == Vector2(0.5, 1.0)
end

l6 = LineSegment(testFace, 5.0)

Test.with_handler(test_handler) do
    @test l6.start == Vector2(0.0, 5.0)
    @test l6.finish == Vector2(0.0, 0.0)
    @test l6.normal == Vector3(-1.0,0.0,0.0)
end

# test binary STL
println("Testing Binary STL import...")
binarymesh = PolygonMesh("./data/cube_binary.stl")
Test.with_handler(test_handler) do
    @test binarymesh.bounds == Bounds(10.0,10.0,10.0,0.0,0.0,0.0)
end

# test mesh rotation
println("Testing Mesh Rotation...")
r1 = rotate(1, 1, 0, 0, 0, 0, 1, 0, 0, pi/2)
r2 = rotate(1, 1, 0, 0, 0, 0, 0, 1, 0, pi/2)
r3 = rotate(1, 1, 0, 0, 0, 0, 0, 0, 1, pi/2)
rr1 = [1.0, 0, 1.0]
rr2 = [0.0, 1.0, -1.0]
rr3 = [-1.0, 1.0, 0.0]
Test.with_handler(test_handler) do
    for i = 1:3
        @test isapprox(r1[i] ,rr1[i])
        @test isapprox(r2[i] ,rr2[i])
        @test isapprox(r3[i] ,rr3[i])
    end
end

rotate!(binarymesh, 45.0, [1.0,0.0,0.0], [0.0,0.0,0.0])
Test.with_handler(test_handler) do
    @test binarymesh.bounds == Bounds(10.0,10.0,13.762255133518483,0.0,-8.509035245341185,0.0)
end
rotate!(binarymesh, 45.0, [0.0,1.0,0.0], [0.0,0.0,0.0])
Test.with_handler(test_handler) do
    @test binarymesh.bounds == Bounds(16.96357128682594,10.0,13.762255133518483,0.0,-8.509035245341185,-8.509035245341185)
end
