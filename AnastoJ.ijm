// Parameters to check before running the macro
minBranchArea = 10;



//setBatchMode(true);
// Clear ROI manager and close all opened images
roiManager("reset");
close("*");

// Ask for input directory
inputDir = getDirectory("Please select a directory containing images to analyze");

// Check poreTemplate.tif file is in input directory
if(!File.exists(inputDir + "poreTemplate.tif")){
	exit("Please download poreTemplate.tif file and drop it into your images directory.");
}

// Get all files in input directory
imgList = getFileList(inputDir);


// Create a global output directory
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
outDir = inputDir+"Results_"+(month+1)+"-"+dayOfMonth+"_"+hour+"-"+minute+File.separator();
if (!File.isDirectory(outDir)) {
	File.makeDirectory(outDir);
}

 
// Create a results csv file
resultsFile = File.open(outDir+ "results_skeleton.csv");
print(resultsFile,"Image name, Pore_ID , Branch Length, Total branch length , Euclidean branch length , Total Euclidean branch length , Branches nb, Junctions nb, Endoints nb , Nb of rejected branches < "+minBranchArea+" \n");
File.close(resultsFile);

frontiersResultsFile = File.open(outDir+"results_frontiers.csv");
print(frontiersResultsFile,"Image name, Anastomosis area (μm2), Pore_ID 1, Pore_ID 2\n");
File.close(frontiersResultsFile);
// Set foreground and background colors
Color.setForeground("white");
Color.setBackground("black");

// Loop through all files with .nd2 extension in input directory
for (i = 0; i < imgList.length; i++) {
	if (endsWith(imgList[i], ".nd2")) {
		print("Analyzing image " + imgList[i] + "...");
		
    	// Create a specific output directory based on image name
		imgName = replace(imgList[i], ".nd2", "");
		
		run("Bio-Formats Importer", "open=["+inputDir+"poreTemplate.tif] color_mode=Default view=Hyperstack stack_order=XYCZT");
    	getDimensions(poreWidth, poreHeight, channels, slices, frames);
		
    	// Open image
    	run("Bio-Formats Importer", "open=["+inputDir+imgList[i]+"] color_mode=Default view=Hyperstack stack_order=XYCZT");
    	

    	
    	// Z-project
    	run("Z Project...", "projection=[Max Intensity]");
    	Stack.setChannel(3);
    	
    	
    	close(imgList[i]);
    	
    	run("Set Measurements...", "area mean redirect=None decimal=2");
 	 	// Split channels
    	run("Split Channels");
    	close("C1-MAX_"+imgList[i]);
    	
    	/*
    	 * GFP segmentation
    	 */
    	 
    	selectImage("C2-MAX_"+imgList[i]);
    	getPixelSize(unit, pixelWidth, pixelHeight);
    	getDimensions(widthImage, heightImage, channels, slices, frames);
    	// Median Filter 
    	run("Median...", "radius=3");
    	// automatic Mean threshold
		setAutoThreshold("Mean dark");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		// fill Holes betwenn 0 and 50 pixels
		fillHoles(0, 50);
		// Keep only the big pores > 1000 pixels and remove error in segmentation
    	run("Analyze Particles...", "size=1000-Infinity show=Masks");
    	// get the segmentatin mask as a variable
		run("Invert LUT");
		rename("C2-Mask");
		titre_mask = getTitle();
		
  ////////////////////////////////
 //   Template-Matching Part   //
////////////////////////////////

		// Use the template matching method to find pores in Brightfield image
    	run("Template Matching Image", "template=poreTemplate.tif image=[C3-MAX_"+imgList[i]+"] flip_template_vertically flip_template_horizontally rotate=[] matching_method=[Normalised 0-mean cross-correlation] number_of_objects=30 score_threshold=0.4 maximal_overlap=0 add_roi ");
    	selectImage("C3-MAX_"+imgList[i]);
    	Width = getWidth();
    	run("Set Measurements...", "mean standard redirect=None decimal=2");
    	// Create a dark image for later
    	newImage("test", "8-bit black", Width, Width, 1);
    	// For each pore check if the standard deviation is > 10000, ensure the presence of a pore 
    	for (r = 0; r < roiManager("count"); r++) {
    		selectImage("C3-MAX_"+imgList[i]);
    		roiManager("select", r);
    		
    		List.setMeasurements;
    		std_dev =List.getValue("StdDev");
    		 
    		selectImage(titre_mask);
    		
    		 List.setMeasurements;
    		mean_GFP = List.getValue("Mean");
    		if (std_dev < 10000 || mean_GFP == 0) {
    			roiManager("delete");
    			r=r-1;
    		} else {
    			// If there is a pore, create an ellipse at the center of it and fill it in white in the drak image we create before
    			// this will be used for the watershed-Split algo as seeds for each pore.
    			selectImage(titre_mask);
    			roiManager("select", r);
    			Roi.getBounds(x, y, width, height);
    			centerX = x+width/2;
    			centerY = y+height/2;
    			
				selectImage("test");
				makeEllipse(centerX-20, centerY-20, centerX+20, centerY+20, 1);
				// Fill it in white
				setColor(255);
				run("Fill");
				run("Select None");
    		}
    	}
    	
    	nb_pores = roiManager("count");
    	

    	
  //////////////////////////
 //   Anastomosis Part   //
//////////////////////////

    	// Use the watershed split algo allowing us to separate segmented pore in GFP channels
    	run("3D Watershed Split", "binary="+titre_mask+" seeds=test radius=2");
    	pore_mask = getTitle();
    	run("Duplicate...", " ");
    	setThreshold(1, 65000);
    	setOption("BlackBackground", true);
		run("Convert to Mask");
		
    	run("Outline");
    	title_outline = getTitle();
    	// Multyply both the mask and the outline to obtain the frontier between pores
    	imageCalculator("Multiply create", title_outline,titre_mask);
    	// Extend these frontiers by one pixel to calculate the value of each associated pore
    	run("Outline");
    	run("Fill Holes");
    	// Exclude frontier with the edges of the image
    	run("Analyze Particles...", "add exclude");
    	run("Set Measurements...", "area redirect=None decimal=2");
    	run("Clear Results");
    	
    	for (a = nb_pores; a < roiManager("count"); a++){
    		roiManager("select", a);
    		 List.setMeasurements ;
    		area_frontier = List.getValue("Area");
    		run("Clear Results");
    		// Return all the pixels composing the frontier 
    		Roi.getContainedPoints(xpoints, ypoints);
    		
    		value_array = newArray();
    		res=newArray();
    		selectImage("Split");
    		// Loop throught each pixel and get his value from the watershed segmented image
    		for (f = 0; f < xpoints.length; f++) {
    			value_array[f] = getValue(xpoints[f], ypoints[f]);
    		}
    		// We have now x pixels in the array but only 3 values, next step will create a new array with only 3 values and sort it as res[0] = 0 (background) res[1] = pore_ID 1 and res[2] = pore_ID 2
    		res = getUnique(value_array);
    		// Write results 
    		if (res.length > 2) {
    		frontier_params = imgName+","+(area_frontier*pixelWidth)+","+res[1]+","+res[2]; 
			File.append(frontier_params, outDir + "results_frontiers.csv");	
    		}
    	}
    		// Loop throught roi to delete the frontier ROI
    	while (roiManager("count")!= nb_pores) {
    	
    		roiManager("select", nb_pores);
			roiManager("delete");
    	}
    	
    	
  ///////////////////////
 //   Skeleton Part   //
///////////////////////

		
    	nb_init_pores = roiManager("count");
    	Euc_d = 0;
    	total_distance = 0;
    	// Loop thought each pores 
		for (pore_ID = 1; pore_ID <= nb_init_pores; pore_ID++){
			selectImage(pore_mask);
			run("Select None");
			run("Duplicate...", " ");
			// Set specific Threshold to segment pore one by one
			setThreshold(pore_ID, pore_ID, "raw");
			setOption("BlackBackground", true);
			run("Convert to Mask");
			run("Skeletonize");
			skelImgName = getTitle();
			// Return the skeleton without junctions so we can filter out short branch
			getSkeletonWithoutJunctions(true);
	
			run("Set Measurements...", "area min redirect=[Tagged skeleton] decimal=2");
			
			run("Analyze Particles...", "size=0-Infinity show=Nothing display add record");
			
			selectWindow(skelImgName);
			setBackgroundColor(0, 0, 0);
			d=0;
			total_Ed = 0;
			branch_before = roiManager("count")-nb_init_pores;
			
			// Loop throught each branch and remove the one < than the min branch length in um
			for (m = 0 ; m < branch_before; m++) {		
				if(getResult("Area", m) < (minBranchArea*pixelWidth) && getResult("Min", m) == 30) {
					roiManager("select", nb_init_pores+m);
					run("Clear", "slice");
				}
			}
			
			run("Skeletonize");
	
			close(skelImgName+"NoJct");
			close("Tagged skeleton");
			nbResults = roiManager("count")-nb_init_pores;
			
			for (n = 0 ; n < nbResults; n++) {		
				roiManager("select", nb_init_pores);
				roiManager("delete");	
			}
			// Analyse the skeleton and get results
			run("Skeletonize");
			run("Analyze Skeleton (2D/3D)", "prune=none show");
			if (branch_before != 0) {
				selectWindow("Branch information");
				c = Table.getColumn("Branch length");
				Ed = Table.getColumn("Euclidean distance");
				for (b=0 ; b < c.length; b++) {
					d = d+(c[b]*pixelWidth);
					total_distance = total_distance + (c[b]*pixelWidth);
				
					total_Ed = total_Ed+(Ed[b]*pixelWidth);
					Euc_d = Euc_d + (Ed[b]*pixelWidth);
			}
			total_branches = getResult("# Branches", 0);
			total_junctions = getResult("# Junctions", 0);
			total_endpoints = getResult("# End-point voxels", 0);
			diff_branches = branch_before - total_branches;
			branch_before = 0;
			
			params = imgName+","+pore_ID+","+d+","+total_distance+","+total_Ed+","+Euc_d+","+total_branches+","+total_junctions+","+total_endpoints+","+diff_branches; 
						File.append(params, outDir + "results_skeleton.csv");
			}
			close();
			close("Split-*");
			}
			selectImage(pore_mask);
			Overlay.clear;
			roiManager("reset");
			saveAs("tiff", outDir+imgName);
			
			close("*");
			roiManager("reset");
		}
}
		
		

function getUnique(a) {
    // tableau résultat et compteur
    unique = newArray();
    o = 0;

    for (i = 0; i < lengthOf(a); i++) {
        val = a[i];
        found = false;
        for (j = 0; j < o; j++) {
            if (unique[j] == val) {
                found = true;
                break;
            }
        }
        if (!found) {
            // agrandir manuellement unique
            tmp = newArray(o + 1);
            for (k = 0; k < o; k++)
                tmp[k] = unique[k];
            tmp[o] = val;
            unique = tmp;
            o++;
        }
    }

    // Si on veut trier (optionnel)
    if (o > 1) Array.sort(unique);

    // Stocke le tableau dans une variable globale fiable
    GETUNIQUE_RESULT = unique;

    // Retourne une chaîne sûre pour affichage / debugging
    return unique;
}

function getSkeletonWithoutJunctions(keepTaggedSkel) {
	outputImgName = getTitle()+"NoJct";
	run("Duplicate...", "title="+outputImgName);
	
	run("Analyze Skeleton (2D/3D)", "prune=none");
	run("Duplicate...", "title=jct");
	setThreshold(60, 80, "raw");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	imageCalculator("Subtract", outputImgName, "jct");
	
	selectWindow(outputImgName);
	close("jct");
	close("Results");
	
	if(!keepTaggedSkel) {
		close("Tagged skeleton");
	}
}  	

// Fill holes whose area is between minHoleArea and maxHoleArea
function fillHoles(minHoleArea, maxHoleArea) {
	run("Invert");
	run("Analyze Particles...", "size=" + minHoleArea + "-" + maxHoleArea + " add");
	run("Invert");
	if (roiManager("count") > 0) {
		roiManager("deselect");
		roiManager("combine");
		run("Fill", "slice");
		roiManager("reset");
		run("Select None");
	}
}