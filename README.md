# Mesh Improvement via Sampling:

Implementation and results of "A Constrained Resampling Strategy for Mesh Improvement", Proceedings of 
The Eurographics Symposium on Geometry Processing (SGP) 2017.

# Description:

The code includes the implementation for curved surface meshes. This includes the following applications:
1) _Non-obtuse retriangulation_:
   Eliminating obtuse triangles from input mesh while preserving the minimum angle bound, smoothness, and Delaunay property. 

<img src="https://github.com/Ahdhn/MeshImp/blob/master/data/input/gargoyle/input.png" height="435" width="435"><img src="https://github.com/Ahdhn/MeshImp/blob/master/data/NonObtuse/gargoyle/output.png" height="435" width="435">

```
	Input (Obtuse triangles in red)        					   Output
```


2) _Mesh Simplification_:
   Achieving Level of Detail while preserving the angle bound, and triangle quality.

<img src="https://github.com/Ahdhn/MeshImp/blob/master/data/input/cvt/davidhead/input.png" height="215" width="215"><img src="https://github.com/Ahdhn/MeshImp/blob/master/data/DelaunaySifting/cvt/davidhead/output.png" height="215" width="215"><img src="https://github.com/Ahdhn/MeshImp/blob/master/data/MeshSimplification/cvt/davidhead/level_1.png" height="215" width="215"><img src="https://github.com/Ahdhn/MeshImp/blob/master/data/MeshSimplification/cvt/davidhead/level_2.png" height="215" width="215">

```
	30K triangles                20K triangles                12K triangles                  1.5K triangles
```


3) _Delauany Sifting_:
   Reducing the number of Steiner points while preserving angle bound, triangle quality, smoothness, and Delaunay property. 

<img src="https://github.com/Ahdhn/MeshImp/blob/master/data/input/dr/bimba/input.png" height="435" width="435"><img src="https://github.com/Ahdhn/MeshImp/blob/master/data/DelaunaySifting/dr/bimba/output.png" height="435" width="435">

```
		      Input         				     	      Output
```

The results as presented in the paper in Table 2, Table 3, Table 5 are included as well under `/data`. 




# Installation:

The code has no prerequisites except CMake (3.6 or higher). The code was developed and tested under Microsoft Visual 12.0 C++11 on Windows 7 (x64).

- Windows:
	* On CMake GUI, set the source code directory to where the CMakeLists.txt is and set the build directory to where you would like the project to be built. Click `Configure` followed by `Generate` followed by `Open Project`. We specified the generator as "Visual Studio 12 2013 Win64". Other generators should work too.

	* On Visual Studio, Set the Solution Configuration to Release. Under Build menu, 	click Build Solution. The executable should then be found under `/BUILD_DIR/Release`.

# Usage:
The code relies on switches to specify the right application and objectives to improve as follows:

#### Command Line Syntax is:
```
mesh_imp.exe -APP [-tar -sam -smooth -dih -ring -del -minang -maxang] INPUT.obj
```
`-APP` could be `-sim` for Mesh Simplification or `-obt` for Non-obtuse Retriangulation. The code accepts .obj files that contains only raw vertices (preceded by `v`) and faces (preceded by `f`) in addition to comments (preceded by `#`).

#### Command Line Switches are:
```
   -sim      Invokes the Mesh Simplification algorithm on the input mesh
             where the complexity of the input mesh is reduced while
             respecting the specified mesh quality.
             
   -tar      The number followed represents the target number of samples desired
             in the simplified output mesh. If this number is not set, 
             The code will attempt to reduce the samples count as much as possible.

   -obt      Invokes the non-obtuse retriangulation algorithm on the input
             mesh where we attempt to eliminate the obtuse triangles from 
             the input mesh.

   -sam      The number followed represents the number of successive successful
             darts before termination.
	     During the sampling process, we can sample more than one successful
	     sample and pick the best one. The best one can have different objective.
	     In the current implementation, the `best' here is the one with maximum.
	     minimum apex angle across different algorihtms.
             The default is 10 samples.
			 
   -ring     The number followed represents the number of rings (layers)
             taken as background grid for sampling. The default is 3 rings.

   -del      If set, the Delaunay property will be preserved.
             The default is false.

   -minang   The number followed represents the minimum angle preserved.
             The default is set to the minimum angle of the input mesh.
			 
   -maxang   The number followed represents the maximum angle preserved.
             The default is set to the maximum angle of the input mesh.   

   -smooth   If set, the smoothness will be preserved.
             The default is false.
			 
   -dih      The maximum dihedral angle to which the smoothness will be
             preserved.
			 
   -v        Display various statistics throughout the execution.
   
   -h        Display the use message, and quit.
  ```

### Examples

Non-obtuse Retriangulation:

```
mesh_imp.exe -obt -sam10 -smooth -dih170 -del -minang  INPUT.obj
```

Delaunay Sifting:
```
mesh_imp.exe -sim -sam10 -smooth -dih170 -del -minang -maxang INPUT.obj
```

Mesh Simplification:
```
mesh_imp.exe -sim -sam10 -minang -maxang INPUT.obj
```


# Reporting problems:

To report bugs or compilation issue, please file an issue [here](https://github.com/Ahdhn/MeshImp/issues).
