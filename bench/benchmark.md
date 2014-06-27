# Before Optimization

steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 2.045935488 seconds (833052088 bytes allocated)
elapsed time: 0.86993292 seconds (351221616 bytes allocated)
elapsed time: 0.527344829 seconds (64522600 bytes allocated)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 2.687707098 seconds (884174624 bytes allocated)
elapsed time: 0.899445194 seconds (351215552 bytes allocated)
elapsed time: 0.527980913 seconds (64513028 bytes allocated)

# With Immutable bounds, not good.

steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 2.716703447 seconds (899433748 bytes allocated)
elapsed time: 0.913785208 seconds (360521920 bytes allocated)
elapsed time: 0.536434056 seconds (64513028 bytes allocated)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 2.740043702 seconds (900698276 bytes allocated)
elapsed time: 0.945200606 seconds (360521920 bytes allocated)
elapsed time: 0.593300534 seconds (64513028 bytes allocated)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 2.66593846 seconds (899349700 bytes allocated)
elapsed time: 0.899419326 seconds (360521920 bytes allocated)
elapsed time: 0.576144777 seconds (64513028 bytes allocated)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ 


# with faces as linked list

steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 2.036814304 seconds (830264312 bytes allocated)
elapsed time: 0.880437183 seconds (359156448 bytes allocated)
elapsed time: 0.490653605 seconds (69616416 bytes allocated)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 2.805084688 seconds (882807912 bytes allocated)
elapsed time: 0.876525983 seconds (359150384 bytes allocated)
elapsed time: 0.536929522 seconds (69607004 bytes allocated)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 2.747059389 seconds (881580412 bytes allocated)
elapsed time: 0.894625001 seconds (359150384 bytes allocated)
elapsed time: 0.488832929 seconds (69607004 bytes allocated)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ 

# remove sorting of vertices

steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 1.218973974 seconds (448063936 bytes allocated)
elapsed time: 0.266686171 seconds (51653052 bytes allocated)
elapsed time: 0.372984835 seconds (25322384 bytes allocated)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 1.208036712 seconds (448024988 bytes allocated)
elapsed time: 0.268237677 seconds (51653052 bytes allocated)
elapsed time: 0.373860517 seconds (25322384 bytes allocated)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 1.222956401 seconds (448092096 bytes allocated)
elapsed time: 0.294331718 seconds (51653052 bytes allocated)
elapsed time: 0.383517176 seconds (25322384 bytes allocated)

# Break early in height list
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 0.919651744 seconds (415934104 bytes allocated)
elapsed time: 0.251398088 seconds (50586416 bytes allocated)
elapsed time: 0.260027178 seconds (82725620 bytes allocated)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 0.966058951 seconds (415936568 bytes allocated)
elapsed time: 0.254233233 seconds (50586416 bytes allocated)
elapsed time: 0.264334856 seconds (82725620 bytes allocated)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 0.924111572 seconds (415936168 bytes allocated)
elapsed time: 0.253502532 seconds (50586416 bytes allocated)
elapsed time: 0.260585732 seconds (82725620 bytes allocated)

# using DataStructures.jl

steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 0.904120562 seconds (418402728 bytes allocated)
elapsed time: 0.280351949 seconds (50606828 bytes allocated)
elapsed time: 0.329599369 seconds (82709456 bytes allocated)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 0.917570881 seconds (418403080 bytes allocated)
elapsed time: 0.285277793 seconds (50606828 bytes allocated)
elapsed time: 0.342525914 seconds (82709456 bytes allocated)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia-git slice_tree.jl 
elapsed time: 0.899747286 seconds (418407304 bytes allocated)
elapsed time: 0.286881589 seconds (50606828 bytes allocated)
elapsed time: 0.336618854 seconds (82709456 bytes allocated)

# With Polygon splicing
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia slice_tree.jl 
INFO: Loading help data...
First mesh load
elapsed time: 1.090027502 seconds (425745616 bytes allocated, 32.56% gc time)
First rotation
elapsed time: 0.301161317 seconds (71385580 bytes allocated, 17.39% gc time)
Second mesh load
elapsed time: 1.012526352 seconds (407552960 bytes allocated, 55.04% gc time)
Second rotation
elapsed time: 0.209350957 seconds (66227520 bytes allocated, 28.22% gc time)
First Slice
elapsed time: 0.345537349 seconds (30581936 bytes allocated)
Second Slice
elapsed time: 0.290897034 seconds (25900720 bytes allocated, 20.59% gc time)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia slice_tree.jl 
First mesh load
elapsed time: 1.01764528 seconds (425456808 bytes allocated, 26.55% gc time)
First rotation
elapsed time: 0.338779086 seconds (71248392 bytes allocated, 26.63% gc time)
Second mesh load
elapsed time: 0.950085899 seconds (407552960 bytes allocated, 53.94% gc time)
Second rotation
elapsed time: 0.212357985 seconds (66227520 bytes allocated, 26.03% gc time)
First Slice
elapsed time: 0.392493393 seconds (30361584 bytes allocated, 13.49% gc time)
Second Slice
elapsed time: 0.23300473 seconds (25900720 bytes allocated)

# With Arrays instead of LinkedList
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia slice_tree.jl 
First mesh load
elapsed time: 0.965414743 seconds (425867484 bytes allocated, 23.73% gc time)
First rotation
elapsed time: 0.294804637 seconds (66166284 bytes allocated, 23.38% gc time)
Second mesh load
elapsed time: 0.797142521 seconds (407927472 bytes allocated, 44.25% gc time)
Second rotation
elapsed time: 0.173766517 seconds (61133120 bytes allocated, 24.52% gc time)
First Slice
elapsed time: 0.325421391 seconds (25262964 bytes allocated)
Second Slice
elapsed time: 0.260986107 seconds (20807272 bytes allocated, 16.69% gc time)
steve@sjkellyT420:~/.julia/MeshSlicer/examples$ julia slice_tree.jl 
First mesh load
elapsed time: 0.938301709 seconds (425876704 bytes allocated, 24.10% gc time)
First rotation
elapsed time: 0.291845671 seconds (66166284 bytes allocated, 23.09% gc time)
Second mesh load
elapsed time: 0.743311234 seconds (407927472 bytes allocated, 40.87% gc time)
Second rotation
elapsed time: 0.206246307 seconds (61133120 bytes allocated, 38.10% gc time)
First Slice
elapsed time: 0.326519483 seconds (25262964 bytes allocated)
Second Slice
elapsed time: 0.261493267 seconds (20807272 bytes allocated, 16.69% gc time)

