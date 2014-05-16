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


# getvolume/getface test
println("Testing Volume...")

# get the first face from cube.stl
f = open("./data/cube.stl","r")

testFace = MeshSlicer.Face(Array[[0.0,0.0,10.0],
                                 [0.0,10.0,10.0],
                                 [0.0,0.0,0.0]],
                           [-1.0,0.0,0.0],
                           [3,1,1])

testBounds = MeshSlicer.Bounds(10.0,10.0,10.0,0,0,0)

resultVolume = MeshSlicer.Volume(f)

Test.with_handler(test_handler) do
    @test resultVolume.faces[1].vertices == testFace.vertices
    @test resultVolume.faces[1].normal == testFace.normal
    @test resultVolume.faces[1].heights == testFace.heights
end

println("Testing Bounds...")
Test.with_handler(test_handler) do
    @test resultVolume.bounds.xmax == testBounds.xmax
    @test resultVolume.bounds.ymax == testBounds.ymax
    @test resultVolume.bounds.zmax == testBounds.zmax
    @test resultVolume.bounds.xmin == testBounds.xmin
    @test resultVolume.bounds.ymin == testBounds.ymin
    @test resultVolume.bounds.zmin == testBounds.zmin
end

close(f)

# Test bounds on translated cube
# translate([1,2,3])cube([4,5,6]);
f = open("./data/translated_cube.stl","r")
testBounds2 = MeshSlicer.Bounds(5.0,7.0,9.0,1.0,2.0,3.0)
resultVolume2 = MeshSlicer.Volume(f)
Test.with_handler(test_handler) do
    @test resultVolume2.bounds.xmax == testBounds2.xmax
    @test resultVolume2.bounds.ymax == testBounds2.ymax
    @test resultVolume2.bounds.zmax == testBounds2.zmax
    @test resultVolume2.bounds.xmin == testBounds2.xmin
    @test resultVolume2.bounds.ymin == testBounds2.ymin
    @test resultVolume2.bounds.zmin == testBounds2.zmin
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

