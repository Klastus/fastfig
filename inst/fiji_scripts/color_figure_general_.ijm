// Create a figure with jpeg images next to the other 
// based on the images from stack

// Use R to iterate through appropriate folders 
// Each stack will be subjected to the following Fiji macro.

LenOfArgs = 8;
args = getArgument();
args = split(args, ";");

if(args.length != LenOfArgs && args.length < LenOfArgs){
	print("Not enough arguments");
	exit;
}


input = args[0];
output = args[1];
StackName = args[2];
min = split(args[3], ",");
max = split(args[4], ",");
DyeList = split(args[5], ",");
NameList = split(args[6], ",");
ScaleBarWidth = args[7];
ScaleBarOnSeparate = args[8];


//input = "D:/Piotrek/Experiments/ICF/TNF-immuno_fusion/images/15min_10ng/same/cell_test/";
//output = "D:/Piotrek/Experiments/ICF/TNF-immuno_fusion/images/15min_10ng/same/cell_test/";
//DyeList = newArray("Alexa 488 - Confocal_529-24 - n000000", 
//"Alexa 555 Confocal fusion - n000000", 
//"Alexa 647 - Confocal fusion - n000000",
//"DAPI Confocal - n000000"); //order determines colors!
//max = newArray(1000, 600, 250, 700);
//NameList = newArray("488", "green", "red", "blue")

setBatchMode(true);
open(input + StackName);

run("Stack to Images");

format = ".jpg";
for( image=0; image < DyeList.length; image++){
    selectWindow(DyeList[image]);
    setMinAndMax(min[image], max[image]);
}

MergeCommand = "";
/*ImageJColors = newArray("red", "green", "blue", "grey");
for(image = 0; image < DyeList.length; image++){
    UserColorName = NameList[image];
    for(color = 0; color < ImageJColors.length; color++){
        if(ImageJColors[color] == UserColorName){
            UserColorName = color+1;
            break;
        }
    }
    MergeCommand = MergeCommand + "c" + UserColorName + "=[" + DyeList[image] + "] ";
}
*/
ImageJColors = newArray("red", "green", "blue");
image = 0;
for(name = 0; name < NameList.length; name++){
    UserColorName = NameList[name];
    for(color = 0; color < ImageJColors.length; color++){
            if(ImageJColors[color] == UserColorName){
                UserColorName = color+1;
                MergeCommand = MergeCommand + "c" + UserColorName + "=[" + DyeList[image] + "] ";
                image++;
            }
        }
}
    
    
MergeCommand = MergeCommand + "create";
//MergeCommand = "c1=[" + DyeList[0] + "] c2=[" + DyeList[1] + "] c3=[" + DyeList[3] + "] create";
run("Merge Channels...", MergeCommand);


for(name = 0; name < NameList.length; name++){
    if(NameList[name] == "Composite"){
        OutputPath = output + "Composite" + format;
        saveAs("Jpeg", OutputPath);
    }
}

run("Set Scale...", "distance=1018 known=0.3 pixel=1 unit=mm");
run("Stack to Images");

WhichGrey = 0;
for(name = 0; name < NameList.length; name++){
    if(NameList[name] == "grey"){
        selectWindow(DyeList[WhichGrey]);
        rename("Composite-000" + DyeList.length);
    }
    if(NameList[name] != "Composite"){
        WhichGrey++;
    }
}

ImageJColors = newArray("red", "green", "blue", "grey");
NameColorValid = 0;
for(color=0; color < ImageJColors.length; color++){ 
    for(name = 0; name < NameList.length-1; name++){//NameList.length-1 because "montage_name" should be avoided and is the last element
        if(NameList[name] == ImageJColors[color]){
            CompositeNumber = NameColorValid +1;
            CompositeName = "Composite-000" + CompositeNumber;
            selectWindow(CompositeName);
            OutputPath = output + NameList[name] + format;
            
            if(ScaleBarOnSeparate == 1){
                run("Scale Bar...", "width=" + ScaleBarWidth + " height=10 font=14 color=White background=None location=[Lower Right] hide");
            }
            saveAs("Jpeg", OutputPath);
            NameColorValid++;
        }
    }
}

run("Close All");


for( image=0; image < (NameList.length-1); image++){
	open(input + NameList[image] + ".jpg");
}


run("Images to Stack", "name=Stack title=[] use");
run("Set Scale...", "distance=1018 known=0.3 pixel=1 unit=mm");
setSlice(NameList.length-1);
//ScaleCommand = "width=" + ScaleBarWidth + " height=4 font=14 color=White background=None location=[Lower Right] hide";
run("Scale Bar...", "width=" + ScaleBarWidth + " height=10 font=14 color=White background=None location=[Lower Right] hide");

MontageCommand = "columns=" + NameList.length-1 + " rows=1 scale=1 border=5";
run("Make Montage...", MontageCommand);

OutputPath = output + NameList[NameList.length-1] + format;
saveAs("Jpeg", OutputPath);
run("Close All");
run("Quit");

exit();


