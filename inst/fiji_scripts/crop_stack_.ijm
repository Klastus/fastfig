// Create a Stack.tif based on a Roi 
// 488, 555, 647 and DAPI (fusion experiments)

// Use R to iterate through appropriate folders 
// Each roi will be subjected to the following Fiji macro.


LenOfArgs = 4;
args = getArgument();
args = split(args, ";");

if(args.length != LenOfArgs && args.length < LenOfArgs){
	print("Not enough arguments");
	exit;
}

input = args[0];
output = args[1];
DyeList = split(args[2], ",");
SpecificRoi = args[3];

/*
input = "D:/Piotrek/Experiments/PEG/raw_images/2019-12-19-PT116_inter/2019-12-19_000/Well F10/";
output = "D:/Piotrek/Experiments/ICF/IFNG_fusion/images/15min_10ng/cell_test";
DyeList = newArray("Alexa 488 - Confocal_529-24 - n000000", 
"Alexa 555 Confocal fusion - n000000", 
"Alexa 647 - Confocal fusion - n000000",
"DAPI Confocal - n000000"); //order determines colors!
SpecificRoi = "PT116_centered_F10.roi";
*/

setBatchMode(true);

for( image=0; image < DyeList.length; image++){
    open(input + DyeList[image]);
}

if(DyeList.length > 1){
 run("Images to Stack", "name=Stack title=[] use");
}


roiManager("Open", output + "/" + SpecificRoi);
roiManager("Select", 0);
run("Crop");

OutputPath = output + "/Stack" + ".tif";
saveAs("Tiff", OutputPath);

run("Close All");

run("Quit");

exit();

setBatchMode(false);
