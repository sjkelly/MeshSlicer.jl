#! /usr/bin/env julia

using Benchmark
using MeshSlicer
using DataFrames

outfile = Pkg.dir("MeshSlicer")*"/perf/slice1_benchmarks.csv"
example_dir = Pkg.dir("MeshSlicer")*"/examples/"

function slice1()
    mesh = PolygonMesh(example_dir*"cell_lamp.stl");
    layers = 1000;

    # find slice heights
    step = (mesh.bounds.zmax - mesh.bounds.zmin)/layers

    heights = [mesh.bounds.zmin:step:mesh.bounds.zmax]
    a = MeshSlice(mesh,heights)
end

out = benchmark(slice1, "Slicing", "cell lamp", 10)

if isfile(outfile)
    previous = readtable(outfile)
    out = vcat(previous, out)
end
writetable(outfile, out)

