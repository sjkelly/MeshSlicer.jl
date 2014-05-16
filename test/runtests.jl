using MeshSlicer
using Base.Test

# write your own tests here

# setup test handlers
test_handler(r::Test.Success) = println("Test success: $(r.expr)")
test_handler(r::Test.Failure) = error("Test fail: $(r.expr)")
test_handler(r::Test.Error) = rethrow(r)




# findorder tests
println("Testing findorder...")
a = Array[[1.0,2.0,3.0],[2.0,3.0,1.0],[3.0,1.0,2.0]]
b = Array[[1.0,2.0,3.0],[2.0,2.0,3.0],[2.0,1.0,3.0]]

Test.with_handler(test_handler) do
    @test MeshSlicer.findorder(a, 1) == [1,2,3]
    @test MeshSlicer.findorder(a, 2) == [3,1,2]
    @test MeshSlicer.findorder(a, 3) == [2,3,1]
    @test MeshSlicer.findorder(b, 1) == [1,2,2]
    @test MeshSlicer.findorder(b, 2) == [3,1,1]
    @test MeshSlicer.findorder(b, 3) == [1,1,1]
end


# getmesh/getface test
println("Testing PolygonMesh...")

# get the first face from cube.stl
f = open("./data/cube.stl","r")

testFace = MeshSlicer.Face(Array[[0.0,0.0,10.0],
                                 [0.0,10.0,10.0],
                                 [0.0,0.0,0.0]],
                           [-1.0,0.0,0.0],
                           [3,1,1])

testBounds = MeshSlicer.Bounds(10.0,10.0,10.0,0,0,0)

resultmesh = MeshSlicer.PolygonMesh(f)

Test.with_handler(test_handler) do
    @test resultmesh.faces[1].vertices == testFace.vertices
    @test resultmesh.faces[1].normal == testFace.normal
    @test resultmesh.faces[1].order == testFace.order
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

close(f)

# Test bounds on translated cube
# translate([1,2,3])cube([4,5,6]);
f = open("./data/translated_cube.stl","r")
testBounds2 = MeshSlicer.Bounds(5.0,7.0,9.0,1.0,2.0,3.0)
resultmesh2 = MeshSlicer.PolygonMesh(f)
Test.with_handler(test_handler) do
    @test resultmesh2.bounds.xmax == testBounds2.xmax
    @test resultmesh2.bounds.ymax == testBounds2.ymax
    @test resultmesh2.bounds.zmax == testBounds2.zmax
    @test resultmesh2.bounds.xmin == testBounds2.xmin
    @test resultmesh2.bounds.ymin == testBounds2.ymin
    @test resultmesh2.bounds.zmin == testBounds2.zmin
    @test resultmesh2.bounds == testBounds2
end
close(f)


# Test getlinesegment
testpoints = Array[[0,0,0],[0,1,0],[0,1,1]]
println(MeshSlicer.LineSegment(testpoints[1],testpoints[2],testpoints[3], 0.5))
println(MeshSlicer.LineSegment(testpoints[1],testpoints[3],testpoints[2], 0.5))
println(MeshSlicer.LineSegment(testpoints[3],testpoints[1],testpoints[2], 0.5))
println(MeshSlicer.LineSegment(testpoints[3],testpoints[2],testpoints[1], 0.5))
println(MeshSlicer.LineSegment(testpoints[2],testpoints[3],testpoints[1], 0.5))
println(MeshSlicer.LineSegment(testpoints[2],testpoints[1],testpoints[3], 0.5))

