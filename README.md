# Mesh Improvement via Sampling:

Implementation and results of "A Constrained Resampling Strategy for Mesh Improvement" to appear in Proceedings of 
The Eurographics Symposium on Geometry Processing (SGP) 2017.

# Description:

The code includes the implementation for curved surface meshes. This includes the following applications:
- Non-obtuse retriangulation 
- Mesh Simplification 
- Delauany Sifting

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
             of the simplified in the output mesh. If this number is not set,
             The code will attempt to reduce the samples count as much as possible.

   -obt      Invokes the non-obtuse retriangulation algorithm on the input
             mesh where we attempt to eliminate the obtuse triangles from
             the input mesh.

   -sam      The number followed represents the number of successful samples.
             The default is 10 samples.

   -ring     The number followed represents the number of rings (layers)
             taken as background grid for sampling. The default is 3 rings.

   -del      If set, the Delaunay property will be preserved.
             The default is false.

   -minang   The number followed represents the minimum angle preserved.
             The default is set to the minimum angle of the input mesh.
   -maxang   The number followed represents the maximum angle preserved.
             The default is set to the maximum angle of the input mesh.

   -minedge  The number followed represents the minimum edge preserved.
   -maxedge  The number followed represents the maximum edge preserved.

   -smooth   If set, the smoothness will be preserved.
             The default is false.
   -dih      The maximum dihedral angle to which the smoothness will be
             preserved.
   -v        Display various statistics throughout the execution.
   -h        Display the use message.
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