# openVDB_mayaGUI
openVDB_mayaGUI is a Maya plugin that allows you to write openVDB files from Maya Fluids. It can export a frame range, temperature, velocity, but also the density texture and final pixel color. It can also bake the transform into the vdb file. When exporting a range, you can also specify a timestretch allowing you to do slowmotion. It also has a filter system, allowing you to modify the vdb file.

openVDB_mayaGUI works as a command, allowing it to be intergrated into custom scripts and tools. It also has a UI for easier use.

# Usage
There are two ways to use openVDB_mayaGUI. By using the command and by using the UI.

To use the UI, use the command `ovdbUI`. This opens up the UI. You can easily save that command as a shelf item. The UI has annotations on the controls to give you information about them.

Using the command is also quite straight forward.
```
ovdbWrite [-r] [-rs INT] [-re INT] [-rts FLOAT] [-rsn INT] [-f STRING] [-t] [-e] [-c] [-v] [-tf] [-vb] [-p] <-o PATH> [FLUID OBJECTS]
```
| Flag  | Long              | Usage         | Description
| ----- | ----------------- | ------------- | -----------------
| -r    | -range            | -r            | Enable the export of a frame range. -rs and -re need to be specified then!
| -rs   | -rangestart       | -rs INT       | Specify the starting frame. Needed when using the -r flag.
| -re   | -rangeend         | -re INT       | Specify the ending frame, excluding the frame. Needed when using the -r flag.
| -rts  | -rangetimestretch | -rts FLOAT    | Specify the timestretch. 0.5 is 2x slower and 2.0 is 2x faster. The frame number to add in the filename is divided by the timestretch to produce whole numbers.
| -rsn  | -rangestartnumber | -rsn INT      | If timestretch is 0.25 and the starting frame is 2, then the filename starts at frame number 8. This is not always what you would want, so with this you can specify at what frame number the filename should start at. Use -1 to use the starting frame, like normal.
| -f    | -filter           | -f STRING     | Add a filter that should be executed on the grids. You can use this multiple times to add multiple filters. The first filter executes first.
| -e    | -emission         | -e            | Export the temperature grid of the fluid. This is almost always used as an emission, thus the name.
| -v    | -velocity         | -v            | Export the velocity grid of the fluid.
| -t    | -textures         | -t            | Use a different method that samples the fluid, allowing it to export textures. When using this option it will also bake the textures and other shading effects, that affect the density, into the density grid.
| -c    | -color            | -c            | If the flag -t is used, then using this flag it will also export the final pixel color of the fluid. If you have the 'Real Lights' option enabled, it probably also uses those lights. Do be careful with this!.
| -vb   | -verbose          | -vb           | If this is enabled, it will print out all the progress information.
| -p    | -progress         | -p            | If this is enabled, it will try to register a progress window and use that. If there is already one, it will ignore this flag and thus not mess with the existing progress window.
| -o    | -output           | -o PATH       | Specify the output path to write the file(s) to. If you want to have the frame number be inserted into the filename, like you should when exporing a range, you use `$F` in the name. `$F` will be replaced with a six digit whole number, which is the frame number.

## Filters
There are a few filters that you can use.

| Filter            | Description
| ----------------- | -----------------
| voxelize          | This voxelizes the grid. It makes each voxel look like a cube. Great for a pixelated look. This does increase the filesize greatly.
| voxelize_small    | This is the same as voxelize, but it creates a smaller filesize. The quality of the effect is lesser, but often enough.
| resample_{factor} | Change the grid resolution by the specified factor. It uses a quadratic sampler for the resampling.
| blur_{size}       | Blur the grid using a gaussian blur with the specified size.

## Grid names
When exporting the grids, these are how they are names.

| Grid name         | Description
| ----------------- | -----------------
| density           | This is the standard density grid.
| emission          | This is the emission / temperature grid.
| color_{r/g/b}     | This is the color grid. Each channel is its own float grid.
| velocity_{x/y/z}  | This is the velocity grid. Each component is its own float grid.

# Installing
I already have the plugin compiled for Maya 2018 on Windows. If you use a different version or operating system, then you need to recompile it and replace the binary files. Other than that, installing it should be pretty easy.

In the reposity, you have a folder called `maya_plugin`. In here, you find the module file and a files folder with all the files needed. A script that manages the UI and a bin folder with the actual plugin file and the libraries needed. These files are what you need.

1. Make sure that Maya is not running.
2. Create a folder somewhere where you know that it will not be removed. For example on Windows `C:/Program Files/openVDB_mayaGUI`.
3. Place the contents of the folder `maya_plugin` (so not the folder itself), into the just created folder.
4. If you need to recompile the plugin, you place the new plugin files in `CREATED_FOLDER/files/bin`, thus replacing the old ones.
5. Navigate to your Maya environment file, and open it. On Windows you can find it in `Documents/maya/[VERSION]/Maya.env`.
6. Add the following line. ```MAYA_MODULE_PATH = [PATH TO FOLDER CONTAINING THE modFile.mod FILE]```
   If there is already a line with `MAYA_MODULE_PATH`, just append the path to the end by placing a `;` between the path and what was already there. ```MAYA_MODULE_PATH = C:/some/stuff;[PATH TO FOLDER CONTAINING THE modFile.mod FILE]```
7. Done!

# Building
The plugin source code is only one file, which makes it easy to compile it using any tool, since you only need to work with one file.

In order to compile, you need the Maya SDK and the OpenVDB library together with all the libraries that it needs.

When compiling, do make sure that SDL checks is turned off, or the warning level is turned up. There are some warnings which are seen as errors by some compilers.

In the includes you need
* zlib
* ilmbase
* boost
* openvdb-3.2.0 (If you use other versions, there is a good chance that it cannot be read by most render engines.)
* tbb
* Maya SDK

In the linker you need
* OpenMaya.lib                              (Maya)
* OpenMayaFX.lib                            (Maya)
* OpenMayaAnim.lib                          (Maya)
* OpenMayaRender.lib                        (Maya)
* OpenMayaUI.lib                            (Maya)
* Foundation.lib                            (Maya)
* openvdb.lib                               (OpenVDB)
* Half.lib                                  (ILM Base)
* libboost_filesystem_vc141-mt-x64-1_67.lib (BOOST) (Naming will probably be different based on what you use.)

I have found that compiling OpenVDB on windows is very difficult. You do have OpenVDB on vcpkg, but that is OpenVDB 5.0.0, which does not work with most render engines. I ended up creating a Visual Studio project and copying over the files and compiled it from there.
