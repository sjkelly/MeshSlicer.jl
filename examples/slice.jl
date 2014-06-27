using MeshSlicer
using PyPlot

println("First mesh load")
@time mesh = PolygonMesh("cell_lamp.stl");
println("First rotation")
@time rotate!(mesh, 45.0, [1.0,0.0,0.0],[0.0,0.0,0.0])
println("Second mesh load")
@time mesh = PolygonMesh("cell_lamp.stl");
println("Second rotation")
@time rotate!(mesh, 45.0, [1.0,0.0,0.0],[0.0,0.0,0.0])
layers = 100;

# find slice heights
step = (mesh.bounds.zmax - mesh.bounds.zmin)/layers

heights = [mesh.bounds.zmin:step:mesh.bounds.zmax]
# create an array of slices
println("First Slice")
@time a = MeshSlice(mesh,heights)
println("Second Slice")
@time b = MeshSlice(mesh,heights)

for slice in b
    print(length(slice.polygons), " ")
    for poly in slice.polygons
        for line in poly.segments
            plot3D([line.start[1],line.finish[1]], [line.start[2],line.finish[2]], slice.layer)
        end
   end
end


