/*
Copyright 2018 Bram Stout

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the 
Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH 
THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

struct IUnknown;	// Workaround for "combaseapi.h(229): error C2187: syntax error: 'identifier' was unexpected here" when using /permissive-
                    // From https://developercommunity.visualstudio.com/content/problem/185399/error-c2760-in-combaseapih-with-windows-sdk-81-and.html
					// GridTransformer.h indirectly includes combaseapi.h and thus creates that error.

#include <maya/MSimple.h>
#include <maya/MIOStream.h>
#include <maya/MGlobal.h>
#include <maya/MSyntax.h>
#include <maya/MArgDatabase.h>
#include <maya/MSelectionList.h>
#include <maya/MFnFluid.h>
#include <maya/MDagPath.h>
#include <maya/MFileObject.h>
#include <maya/MAnimControl.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MFloatVector.h>
#include <maya/MFloatMatrix.h>
#include <maya/MFloatPointArray.h>
#include <maya/MFloatPoint.h>
#include <maya/MRenderUtil.h>
#include <maya/MDGContextGuard.h>
#include <maya/MCommandResult.h>
#include <maya/MFloatArray.h>
#include <maya/MProgressWindow.h>
#include <maya/MMatrix.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/Dense.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Filter.h>
#include <boost/filesystem.hpp>
#include <random>
#include <string>
#include <cstdlib>
#include <regex>

/*######################################*/
/*										*/
/*		STANDARD BOILERPLATE CODE		*/
/*										*/
/*######################################*/
class ovdbWriteCmd : public MPxCommand{
public:
	ovdbWriteCmd() {};
	~ovdbWriteCmd() override;
	MStatus			doIt(const MArgList& args) override;

	static void*	creator();

	static MSyntax newSyntax();

};

ovdbWriteCmd::~ovdbWriteCmd() {}

void* ovdbWriteCmd::creator(){
	return new ovdbWriteCmd();
}

MStatus initializePlugin(MObject obj) {
	MFnPlugin plugin(obj);
	plugin.registerCommand("ovdbWrite",
		ovdbWriteCmd::creator, ovdbWriteCmd::newSyntax);

	MGlobal::executeCommand("source openvdb_mayagui.mel"); // This registers the script that handles the UI

	return MStatus::kSuccess;
}

MStatus uninitializePlugin(MObject obj) {

	MFnPlugin plugin(obj);
	plugin.deregisterCommand("ovdbWrite");

	return MStatus::kSuccess;
}

/*######################################*/
/*										*/
/*			 Actual code				*/
/*										*/
/*######################################*/

MSyntax ovdbWriteCmd::newSyntax() { // Standard syntax code
	MSyntax syntax;

	syntax.addFlag("-r", "-range");									// Enable the export of a range of frame.
	syntax.addFlag("-rs", "-rangestart", MSyntax::kLong);			// Specify the start of the range.
	syntax.addFlag("-re", "-rangeend", MSyntax::kLong);				// Specify the end of the range.
	syntax.addFlag("-rts", "-rangetimestretch", MSyntax::kDouble);	// Specify the time stretch when exporting a range of frames.
	syntax.addFlag("-rsn", "-rangestartnumber", MSyntax::kLong);	// The frame number in the vdb file's name starts from this number. -1 to disable it.
	syntax.addFlag("-f", "-filter", MSyntax::kString);				// Add a filter to the filter stack. Add this for every filter.
	syntax.makeFlagMultiUse("-f");
	syntax.addFlag("-t", "-textures");								// Enable textures when exporting. This is needed when exporting colors.
	syntax.addFlag("-e", "-emission");								// Export emission which is the temperature channel.
	syntax.addFlag("-c", "-color");									// Export the scatter color.
	syntax.addFlag("-v", "-velocity");								// Export the velocity.
	syntax.addFlag("-tf", "-transform");							// Export the fluid's transform.
	syntax.addFlag("-vb", "-verbose");								// Print out the progress.
	syntax.addFlag("-p", "-progress");								// Show the progress in a progress window.
	syntax.addFlag("-o", "-output", MSyntax::kString);				// Where to save the file. If it is not absolute, the workspace is prefixxed.
	syntax.setObjectType(MSyntax::kSelectionList);
	syntax.useSelectionAsDefault(true);

	return syntax;
}

MStatus ovdbWriteCmd::doIt(const MArgList& args){
	
	// Get all the flags
	MArgDatabase argData(syntax(), args);

	bool verbose = argData.isFlagSet("-vb");
	bool progress = argData.isFlagSet("-p");

	if(progress) // Setup and show the progress window if it is enabled.
		if (!MProgressWindow::reserve())
			progress = false; // If it already exists, then we cannot use the progress window, so disable it
	if (progress) { // If everything still went right, then continue setting up the progress window
		MProgressWindow::setProgressRange(0, 100);
		MProgressWindow::setTitle("ovdbWriter");
		MProgressWindow::setProgress(0);
		MProgressWindow::setInterruptable(true);
		MProgressWindow::startProgress();
	}

	//Continue setting the flags
	bool isRanged = argData.isFlagSet("-r");
	int startFrame = 0;
	int endFrame = 0;
	int startNumber = -1;
	double timeStretch = 1.0;
	if (isRanged) {
		argData.getFlagArgument("-rs", 0, startFrame);
		argData.getFlagArgument("-re", 0, endFrame);
		if(argData.isFlagSet("-rsn"))
			argData.getFlagArgument("-rsn", 0, startNumber);
		if (argData.isFlagSet("-rts")) {
			argData.getFlagArgument("-rts", 0, timeStretch);
			if (timeStretch <= 0.000001) // Prevent division by 0. It can only go a millionth times slower. Should be enough.
				timeStretch = 0.000001;
		}
	}
	startFrame /= timeStretch; // We scale these right now so that we increase the amount of frames to generate.
	endFrame /= timeStretch;   // Later on we multiply by timeStretch when setting the time, to make sure that that stays in the range.
	if (startNumber < 0)       // If the startNumber was not set (thus -1) then set it to the startFrame.
		startNumber = startFrame;

	unsigned int amtFilters = argData.numberOfFlagUses("-f");
	std::vector<std::string> filters(amtFilters);
	for (unsigned int i = 0; i < amtFilters; i++) {
		MString fString;
		argData.getFlagArgument("-f", i, fString);
		filters[i] = std::string(fString.asChar());
	}

	bool isTextured = argData.isFlagSet("-t");
	bool hasEmission = argData.isFlagSet("-e");
	bool hasColor = argData.isFlagSet("-c");
	bool hasVelocity = argData.isFlagSet("-v");

	bool exportTransform = argData.isFlagSet("-tf");

	if (!argData.isFlagSet("-o")) {
		MGlobal::displayError("No output location specified! Please specify the path of the file(s) to write with the flag '-o'.");
		if (progress)  // Make sure to end the progress!
			MProgressWindow::endProgress();
		return MStatus::kFailure;
	}
	MString oString;
	argData.getFlagArgument("-o", 0, oString);
	MString workspaceDir;
	MGlobal::executeCommand("workspace -q -rd", workspaceDir); // Get the project folder
	MFileObject directory;
	directory.setRawFullName(workspaceDir);
	MFileObject absoluteFile;
	absoluteFile.setRawFullName(oString);
	if (!MFileObject::isAbsolutePath(oString)) { // If the path specified was not absolute, add it to the project folder
		MString absoluteFileName = directory.resolvedFullName() + "/" + oString;
		absoluteFile.setRawFullName(absoluteFileName);
	}
	std::string outFile = std::string(absoluteFile.resolvedFullName().asChar());

	MSelectionList objects;
	argData.getObjects(objects); // Get the objects passed to the command
	bool foundFluid = false;
	MFnFluid fluid;
	MDagPath fluidDagPath;
	for (unsigned int i = 0; i < objects.length(); i++) { // Go over every object and use the first fluid
		objects.getDagPath(i, fluidDagPath);
		fluidDagPath.extendToShape(); // If it is a transform, get its child. If it is already a shape, it does nothing
		if (fluidDagPath.apiType() == MFn::kFluid) {
			fluid.setObject(fluidDagPath);
			foundFluid = true;
			break; // We will only allow one fluid at a time
		}
	}
	if (!foundFluid) { // If no fluid was found, then there is no reason to continue.
		MGlobal::displayError("Could not find fluid shape node! Please select a fluid node or append the its name to the command.");
		if (progress)
			MProgressWindow::endProgress();
		return MStatus::kFailure;
	}


	openvdb::initialize(); // This is needed to be called for openvdb to work.
	MTime oldTime = MAnimControl::currentTime(); // Get the current time so that if we use a range, we can set it back at the end

	//If frame ranges is disabled, then startFrame and endFrame will be 0, and the for loop will only execute one. Thus no if statements needed.
	for (int currentFrame = startFrame; currentFrame <= endFrame; currentFrame++) {
		if (verbose)
			MGlobal::displayInfo((std::string("ovdbWrite: frame ") + 
				std::to_string(currentFrame - startFrame) + "/" + std::to_string(endFrame - startFrame)).c_str());

		if (progress) {
			MProgressWindow::setProgressStatus((std::string("Frame ") +
				std::to_string(currentFrame - startFrame) + "/" + std::to_string(endFrame - startFrame)).c_str());
			if (MProgressWindow::isCancelled()) {
				MProgressWindow::endProgress();
				return MStatus::kSuccess;
			}
		}

		if(isRanged)
			MGlobal::viewFrame(((double)currentFrame) * timeStretch); // Set the frame. Make sure to first cast current frame to double when interpolating.

		std::vector<openvdb::FloatGrid::Ptr> gridsTemp; // This is a list of temporary grids, we save them to the right place after the filters
		openvdb::GridPtrVec grids; // The actual list of grids to write
		openvdb::MetaMap meta; // A meta map with data. If you want to add your own metadata, then you can do that here
		meta.insertMeta("creator", openvdb::StringMetadata("openVDB_mayaGUI")); // Just add in that this plugin was used to create the file

		unsigned int xRes = 32;
		unsigned int yRes = 32;
		unsigned int zRes = 32;
		fluid.getResolution(xRes, yRes, zRes); // Get the resolution. There is not really a reason for xRes to be 32, but why not

		double xDim = 32;
		double yDim = 32;
		double zDim = 32;
		fluid.getDimensions(xDim, yDim, zDim); // Get the dimensions

		double xOff = -xDim * 0.5;
		double yOff = -yDim * 0.5;
		double zOff = -zDim * 0.5; // Calculate the offset to center the grids around the origin like with a fluid

		xDim /= (double)xRes;
		yDim /= (double)yRes;
		zDim /= (double)zRes; // We want x/y/zDim to be the size of a voxel

		openvdb::Vec3d voxelSize(xDim, yDim, zDim);
		openvdb::Vec3d containerOffset(xOff, yOff, zOff);
		openvdb::Mat4R mat(voxelSize[0], 0.0, 0.0, 0.0,
			0.0, voxelSize[1], 0.0, 0.0,
			0.0, 0.0, voxelSize[2], 0.0,
			containerOffset[0], containerOffset[1], containerOffset[2], 1.0); // Create the transform for the grids to have it be the same as a fluid container

		if (exportTransform) { // If we want to also bake in the fluid's transformation, then we get it and multiply it with the transformation we already have
			MObject parent = fluid.parent(0);
			MFnDagNode parentNode(parent);
			MMatrix mm = parentNode.transformationMatrix();
			openvdb::Mat4R mat2(mm(0, 0), mm(0, 1), mm(0, 2), mm(0, 3),
				mm(1, 0), mm(1, 1), mm(1, 2), mm(1, 3),
				mm(2, 0), mm(2, 1), mm(2, 2), mm(2, 3),
				mm(3, 0), mm(3, 1), mm(3, 2), mm(3, 3));

			mat *= mat2;
		}

		openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(mat); // Create the transform to pass to the grids

		openvdb::CoordBBox bbox(openvdb::Coord(0), openvdb::Coord(xRes - 1, yRes - 1, zRes - 1)); // Bounding box that shows the voxel coordinate range

		float* data = nullptr; // All the grid data pointers
		bool isDataAllocated = false; // Did we allocate data ourselves and we need to free it?
		float* dataEmission = nullptr;
		bool isDataEmissionAllocated = false;
		float* dataColorR = nullptr;
		float* dataColorG = nullptr;
		float* dataColorB = nullptr;
		bool isDataColorAllocated = false;
		float* dataVelocityX = nullptr;
		float* dataVelocityY = nullptr;
		float* dataVelocityZ = nullptr;
		bool isDataVelocityAllocated = false;


		if (isTextured) { // If we want to export textures we need to use a different method to get the fluid data
			if (hasEmission) // We cannot get the temperature / emission from sampling the fluid as a shading network
				dataEmission = fluid.temperature();

			data = new float[xRes * yRes * zRes]; // We allocate our own data array and fill it in
			isDataAllocated = true;

			if (hasColor) { // Same stuff for color
				dataColorR = new float[xRes * yRes * zRes];
				dataColorG = new float[xRes * yRes * zRes];
				dataColorB = new float[xRes * yRes * zRes];
				isDataColorAllocated = true;
			}
			if (hasVelocity) // Just like temperature / emission we just have to get it as the grid data
				fluid.getVelocity(dataVelocityX, dataVelocityY, dataVelocityZ);
			
			{ // MDGContextGuard works per scope, so we just create a new scope. It is not needed but I would rather have it
				MDGContext context; // The context allows us to change stuff and have it be resetted at the end
				MDGContextGuard tempContext(context);

				MPlug plugFluidPointNear = fluid.findPlug("pointObj"); // Get all the plugs to the attributes
				MPlug plugFluidPointFar = fluid.findPlug("farPointObj");
				MPlug plugFluidOutColor = fluid.findPlug("outColor");
				MPlug plugFluidOutTransparency = fluid.findPlug("outTransparency");
				MPlug plugFluidFilterSize = fluid.findPlug("filterSize");
				MDataHandle dhPointNear, dhPointFar, dhOutColor, dhOutTransparency, dhFilterSize; // Use data handles to access it
				plugFluidPointNear.getValue(dhPointNear); // If we want to set a datahandle, we first have to initialize it like this
				plugFluidPointFar.getValue(dhPointFar);
				plugFluidFilterSize.getValue(dhFilterSize);

				dhFilterSize.set3Float(1.0 / xRes, 1.0 / yRes, 1.0 / zRes); // Set the filter size
				plugFluidFilterSize.setValue(dhFilterSize); // Make sure to update the plug

				float px; // Some variables to be used in the loop. This is the sampling point
				float py;
				float pz;
				float pxFar; // Sampling a fluid works by sampling between to points, so we offset it on the x axis
				float voxelSize = 1.0f / ((float)xRes); // The size of a voxel and also the amount to offset the sampling point for the other point
				double transparency; // A temporary variable to store the transparency

				for (unsigned int i = 0; i < xRes; ++i) { // Three loops to loop over each voxel
					if (verbose)
						MGlobal::displayInfo((std::string("ovdbWrite: frame ") +
							std::to_string(currentFrame - startFrame) + "/" + std::to_string(endFrame - startFrame)
							+ " | Sampling fluid: " + std::to_string(
							(int)((((float) i) / ((float) xRes)) * 100.0f)
							) + "%").c_str());

					if (progress) {
						MProgressWindow::setProgressStatus((std::string("Frame ") +
							std::to_string(currentFrame - startFrame) + "/" + std::to_string(endFrame - startFrame)
							+ " | Sampling fluid").c_str());
						MProgressWindow::setProgress((int)((((float)i) / ((float)xRes)) * 100.0f));
						if (MProgressWindow::isCancelled()) {
							MProgressWindow::endProgress();
							if (isDataAllocated)
								delete[] data;
							if (isDataEmissionAllocated)
								delete[] dataEmission;
							if (isDataVelocityAllocated) {
								delete[] dataVelocityX;
								delete[] dataVelocityY;
								delete[] dataVelocityZ;
							}
							return MStatus::kSuccess;
						}
					}

					px = ((((float)i) + 0.5f) / ((float)xRes)) * 2.0f - 1.0f; // Calculate the sampling point's x value. It is object space which
					pxFar = px + voxelSize;									  // for some reason ignores the fluid's dimension and it from -1.0 to 1.0

					for (unsigned int j = 0; j < yRes; ++j) {
						py = ((((float)j) + 0.5f) / ((float)yRes)) * 2.0f - 1.0f; // Calculate y value

						for (unsigned int k = 0; k < zRes; ++k) {
							pz = (((float)k) / ((float)zRes)) * 2.0f - 1.0f; // Calculate z value

							dhPointNear.set3Float(px, py, pz); // Set the near sampling point
							plugFluidPointNear.setValue(dhPointNear); // Update its plug
							dhPointFar.set3Float(pxFar, py, pz); // Set the far sampling point
							plugFluidPointFar.setValue(dhPointFar); // Update its plug

							plugFluidOutTransparency.getValue(dhOutTransparency); // Retrieve the transparency value
							transparency = dhOutTransparency.asFloat3()[0]; // Set the variable

							data[i + j * xRes + k * xRes * yRes] = 1.0 - transparency; // Save it in our data array

							if (hasColor) { // If we have color, get retrieve the color value and save that in our data array
								plugFluidOutColor.getValue(dhOutColor);
								dataColorR[i + j * xRes + k * xRes * yRes] = dhOutColor.asFloat3()[0];
								dataColorG[i + j * xRes + k * xRes * yRes] = dhOutColor.asFloat3()[1];
								dataColorB[i + j * xRes + k * xRes * yRes] = dhOutColor.asFloat3()[2];
							}
						}
					}
				}
			}
		}
		else { // If we do not need textures, we just use our normal method
			data = fluid.density(); // Retrieve just the grid data, quite simple
			if (hasEmission)
				dataEmission = fluid.temperature();
			if (hasVelocity) 
				fluid.getVelocity(dataVelocityX, dataVelocityY, dataVelocityZ); // Velocity is done like that, since it consists out of three grids
		}
		if (data == nullptr) { // If we for some reason could not actually get the fluid data, we cancle
			MGlobal::displayError("Could not retrieve fluid data!");
			if (progress)
				MProgressWindow::endProgress();
			return MStatus::kFailure;
		}

		if (verbose)
			MGlobal::displayInfo((std::string("ovdbWrite: frame ") +
				std::to_string(currentFrame - startFrame) + "/" + std::to_string(endFrame - startFrame)).c_str());

		{	// Create the grid with everything
			openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0f); // Create the grid with a background value of 0.0

			openvdb::tools::Dense<float, openvdb::tools::LayoutXYZ> dense(bbox, data); // Use the Dense tool to convert a normal grid to a sparse one
			openvdb::tools::copyFromDense(dense, grid->tree(), 1e-7f); // Copy it over, the 1e-7f is a threshold, everything below is just 0.0

			grid->setTransform(transform); // Set the transform to the grid
			grid->setGridClass(openvdb::GRID_FOG_VOLUME); // We only use this plugin for volumetrics, so tell the reader that it is a volume
			grid->setName("density"); // Set the name

			grid->treePtr()->prune(); // Prune it to make the file size smaller

			gridsTemp.push_back(grid); // Save it to the temporary grid list. This list is read by the filters and the filters save it to the actual list
		}
		if (dataEmission != nullptr) {	// Create the grid, but for emission
			openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0f);

			openvdb::tools::Dense<float, openvdb::tools::LayoutXYZ> dense(bbox, dataEmission);
			openvdb::tools::copyFromDense(dense, grid->tree(), 1e-7f);

			grid->setTransform(transform);
			grid->setGridClass(openvdb::GRID_FOG_VOLUME);
			grid->setName("emission");

			grid->treePtr()->prune();

			gridsTemp.push_back(grid);
		}
		if (dataColorR != nullptr) {
			openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0f);

			openvdb::tools::Dense<float, openvdb::tools::LayoutXYZ> dense(bbox, dataColorR);
			openvdb::tools::copyFromDense(dense, grid->tree(), 1e-7f);

			grid->setTransform(transform);
			grid->setGridClass(openvdb::GRID_FOG_VOLUME);
			grid->setName("color_r");

			grid->treePtr()->prune();

			gridsTemp.push_back(grid);
		}
		if (dataColorG != nullptr) {
			openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0f);

			openvdb::tools::Dense<float, openvdb::tools::LayoutXYZ> dense(bbox, dataColorG);
			openvdb::tools::copyFromDense(dense, grid->tree(), 1e-7f);

			grid->setTransform(transform);
			grid->setGridClass(openvdb::GRID_FOG_VOLUME);
			grid->setName("color_g");

			grid->treePtr()->prune();

			gridsTemp.push_back(grid);
		}
		if (dataColorB != nullptr) {
			openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0f);

			openvdb::tools::Dense<float, openvdb::tools::LayoutXYZ> dense(bbox, dataColorB);
			openvdb::tools::copyFromDense(dense, grid->tree(), 1e-7f);

			grid->setTransform(transform);
			grid->setGridClass(openvdb::GRID_FOG_VOLUME);
			grid->setName("color_b");

			grid->treePtr()->prune();

			gridsTemp.push_back(grid);
		}
		if (dataVelocityX != nullptr) {
			openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0);

			openvdb::tools::Dense<float, openvdb::tools::LayoutXYZ> dense(bbox, dataVelocityX);
			openvdb::tools::copyFromDense(dense, grid->tree(), 1e-7f);

			grid->setTransform(transform);
			grid->setGridClass(openvdb::GRID_FOG_VOLUME);
			grid->setName("velocity_x");

			grid->treePtr()->prune();

			gridsTemp.push_back(grid);
		}
		if (dataVelocityY != nullptr) {
			openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0);

			openvdb::tools::Dense<float, openvdb::tools::LayoutXYZ> dense(bbox, dataVelocityY);
			openvdb::tools::copyFromDense(dense, grid->tree(), 1e-7f);

			grid->setTransform(transform);
			grid->setGridClass(openvdb::GRID_FOG_VOLUME);
			grid->setName("velocity_y");

			grid->treePtr()->prune();

			gridsTemp.push_back(grid);
		}
		if (dataVelocityZ != nullptr) {
			openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0);

			openvdb::tools::Dense<float, openvdb::tools::LayoutXYZ> dense(bbox, dataVelocityZ);
			openvdb::tools::copyFromDense(dense, grid->tree(), 1e-7f);

			grid->setTransform(transform);
			grid->setGridClass(openvdb::GRID_FOG_VOLUME);
			grid->setName("velocity_z");

			grid->treePtr()->prune();

			gridsTemp.push_back(grid);
		}

		//Filters
		for (unsigned int j = 0; j < gridsTemp.size(); j++) {
			openvdb::FloatGrid::Ptr grid = gridsTemp[j]; // Retrieve the grid

			for (unsigned int i = 0; i < amtFilters; i++) { // Go through each filter, from the first filter to the last one
				if (verbose)
					MGlobal::displayInfo((std::string("ovdbWrite: frame ") +
						std::to_string(currentFrame - startFrame) + "/" + std::to_string(endFrame - startFrame)
						+ " | Filters: " + std::to_string(
						(int) ((((float)(j * gridsTemp.size() + i)) / ((float) (gridsTemp.size() * amtFilters))) * 100.0f)
						)).c_str());

				if (progress) {
					MProgressWindow::setProgressStatus((std::string("Frame ") +
						std::to_string(currentFrame - startFrame) + "/" + std::to_string(endFrame - startFrame)
						+ " | Filters").c_str());
					MProgressWindow::setProgress((int)((((float)(j * gridsTemp.size() + i)) / ((float)(gridsTemp.size() * amtFilters))) * 100.0f));
					// Filters are fast, so we do not need to check for a cancle
				}

				if (filters[i] == "voxelize") { // If it is the voxelize filter, we create a newer smaller grid and copy the data over using a PointSampler
					const double scaleFactor = 0.125; // The scale factor.

					openvdb::FloatGrid::Ptr gridDest = openvdb::FloatGrid::create(0.0f);

					openvdb::Mat4R matDest(voxelSize[0] * scaleFactor, 0.0, 0.0, 0.0,
						0.0, voxelSize[1] * scaleFactor, 0.0, 0.0,
						0.0, 0.0, voxelSize[2] * scaleFactor, 0.0,
						containerOffset[0], containerOffset[1], containerOffset[2], 1.0); // Create the transform but each voxel is smaller
					openvdb::math::Transform::Ptr transformDest = openvdb::math::Transform::createLinearTransform(matDest);
					openvdb::Mat4R xformVoxelize = grid->transformPtr()->baseMap()->getAffineMap()->getMat4() *
						transformDest->baseMap()->getAffineMap()->getMat4().inverse(); // This is a transform that maps one grid to the other

					openvdb::tools::GridTransformer transformer(xformVoxelize); // Create the instance that copies over the grid data
					transformer.transformGrid<openvdb::tools::PointSampler, openvdb::FloatGrid>(*grid, *gridDest); // Copy it over
					gridDest->setTransform(transformDest->copy()); // Copy over the transfor
					gridDest->setGridClass(grid->getGridClass()); // Copy over the grid class
					gridDest->setName(grid->getName()); // Copy over the name

					gridDest->treePtr()->prune(); // We have a new tree, so prune it for smaller file sizes

					grid = gridDest; // Set the grid pointer to the right grid, so that the next filter can use it
				}else if (filters[i] == "voxelize_small") { // Same thing as before, but the scaleFactor is larger
					const double scaleFactor = 0.25;

					openvdb::FloatGrid::Ptr gridDest = openvdb::FloatGrid::create(0.0f);

					openvdb::Mat4R matDest(voxelSize[0] * scaleFactor, 0.0, 0.0, 0.0,
						0.0, voxelSize[1] * scaleFactor, 0.0, 0.0,
						0.0, 0.0, voxelSize[2] * scaleFactor, 0.0,
						containerOffset[0], containerOffset[1], containerOffset[2], 1.0);
					openvdb::math::Transform::Ptr transformDest = openvdb::math::Transform::createLinearTransform(matDest);
					openvdb::Mat4R xformVoxelize = grid->transformPtr()->baseMap()->getAffineMap()->getMat4() *
						transformDest->baseMap()->getAffineMap()->getMat4().inverse();

					openvdb::tools::GridTransformer transformer(xformVoxelize);
					transformer.transformGrid<openvdb::tools::PointSampler, openvdb::FloatGrid>(*grid, *gridDest);
					gridDest->setTransform(transformDest->copy());
					gridDest->setGridClass(grid->getGridClass());
					gridDest->setName(grid->getName());

					gridDest->treePtr()->prune();

					grid = gridDest;
				}
				else if (filters[i].compare(0, 8, "resample") == 0) { // If we need to resample, we do the same thing as voxelize, but using a QuadraticSampler
					std::string sfactor = filters[i].substr(9); // The scale factor is specified by the user
					const double scaleFactor = std::atof(sfactor.c_str());

					openvdb::FloatGrid::Ptr gridDest = openvdb::FloatGrid::create(0.0f);

					openvdb::Mat4R matDest(voxelSize[0] * scaleFactor, 0.0, 0.0, 0.0,
						0.0, voxelSize[1] * scaleFactor, 0.0, 0.0,
						0.0, 0.0, voxelSize[2] * scaleFactor, 0.0,
						containerOffset[0], containerOffset[1], containerOffset[2], 1.0);
					openvdb::math::Transform::Ptr transformDest = openvdb::math::Transform::createLinearTransform(matDest);
					openvdb::Mat4R xformVoxelize = grid->transformPtr()->baseMap()->getAffineMap()->getMat4() *
						transformDest->baseMap()->getAffineMap()->getMat4().inverse();

					openvdb::tools::GridTransformer transformer(xformVoxelize);
					transformer.transformGrid<openvdb::tools::QuadraticSampler, openvdb::FloatGrid>(*grid, *gridDest);
					gridDest->setTransform(transformDest->copy());
					gridDest->setGridClass(grid->getGridClass());
					gridDest->setName(grid->getName());

					gridDest->treePtr()->prune();

					grid = gridDest;
				}
				else if (filters[i].compare(0, 4, "blur") == 0) { // Blur it. Get the blur size and use the standard filter tool to blur it
					std::string sfactor = filters[i].substr(5);
					const double sizeFactor = std::atof(sfactor.c_str());

					openvdb::tools::Filter<openvdb::FloatGrid> filter(*grid);
					filter.gaussian((int) sizeFactor);
				}
			}

			grids.push_back(grid); // All the filters have been applied to this grid, so just save it to the list of grids to write to the file
		}

		if (verbose)
			MGlobal::displayInfo((std::string("ovdbWrite: frame ") +
				std::to_string(currentFrame - startFrame) + "/" + std::to_string(endFrame - startFrame)
				+ " | Saving...").c_str());

		if (progress) {
			MProgressWindow::setProgressStatus((std::string("Frame ") +
				std::to_string(currentFrame - startFrame) + "/" + std::to_string(endFrame - startFrame)
				+ " | Saving...").c_str());
		}

		std::string fileLocation = outFile; // Replace '$F' with the frame number. If none is found, none is replaced
		std::string frameString = std::to_string(startNumber++); //The start number is also incremented here
		std::string paddedFrameString = std::string(6 - frameString.length(), '0').append(frameString);
		fileLocation = std::regex_replace(fileLocation, std::regex("\\$F"), paddedFrameString);

		// Check if the folders exist. If not, then make them
		boost::filesystem::path dir(fileLocation);
		dir.remove_filename();
		if (!boost::filesystem::exists(dir))
			boost::filesystem::create_directories(dir);
		
		// Create the file object
		openvdb::io::File file(fileLocation);
		file.write(grids, meta); // Write the grid data and the metadata
		file.close(); // Close the file

		// If 'float* data' was created by using the 'new' operator, make sure to free the memeory.
		if (isDataAllocated)
			delete[] data;
		if (isDataEmissionAllocated)
			delete[] dataEmission;
		if (isDataVelocityAllocated) {
			delete[] dataVelocityX;
			delete[] dataVelocityY;
			delete[] dataVelocityZ;
		}

		if (verbose)
			MGlobal::displayInfo((std::string("ovdbWrite: frame ") +
				std::to_string(currentFrame - startFrame) + "/" + std::to_string(endFrame - startFrame)
				+ " | Done").c_str());
	}

	MGlobal::viewFrame(oldTime); // View the frame we used to be on. If ranged was off, we go to the frame we are already one and nothign changes

	if (verbose)
		MGlobal::displayInfo("ovdbWrite: done. ");

	if (progress) {
		MProgressWindow::endProgress();
	}

	return MStatus::kSuccess; // Return
}