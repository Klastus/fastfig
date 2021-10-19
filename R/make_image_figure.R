#### Compose a color image from multiple microscopic images, side by side ####
## It is suitable for preparation of the publication or thesis figure ##


crop_location <- function(roi.input,
                          pattern.roi = "\\[\\].*.roi",
                          path.to.fiji = "D:/Piotrek/programs/Fiji.app",
                          overwrite = TRUE,
                          dyes = list("Alexa 555 Confocal fusion - n000000",
                                      "Alexa 647 - Confocal fusion - n000000",
                                      "DAPI Confocal - n000000",
                                      "Alexa 488 - Confocal_529-24 - n000000"),
                          all.exper.input = "D:/Piotrek/Experiments/PEG/raw_images",
                          stack.name = "Stack.tif"){

  # Important! This function is set up for windows
  # and needs "crop_stack.ijm" file to be present in the ImageJ folder

  # the function to crop an image according to ROI file and creating the Stack.tif.
  # Stack.tif is then used to compose the figure
  # The function searches all .roi files in the provided "roi.input" directory

  # the roi filename should follow this convention:
  # name_of_the_experiment[]wellname.roi
  # [] is a separator. The experiment has to be in all.exper.input folder.
  # The name_of_the_eperiment has to explicitly indicate the experiment
  # in the all.exper.input

  # - dyes; the list of image names taken to stack.tif
  # - overwrite: should the already created Stacks be overwritten?
  #   FALSE is usefull when new roi are created
  # - stack.name: obvious

  setwd(path.to.fiji)
  command.base <- "ImageJ-win64.exe --console -macro crop_stack_.ijm"
  dirs <- dirname(list.files(roi.input,
                             pattern = pattern.roi,
                             full.names = TRUE,
                             recursive = TRUE))
  if(length(dirs) > 0 ){

    for(i in 1:length(dirs)){
      stack.exists <- file.exists(paste(dirs[i], "/", stack.name, sep = ""))
      if(stack.exists == FALSE | overwrite == TRUE){
        roi <- list.files(dirs[i], pattern = pattern.roi, full.names = FALSE, recursive = TRUE)
        metadata <- strsplit(roi, "\\[\\]|\\.")[[1]]

        experiments <- list.dirs(all.exper.input, recursive = FALSE)
        experiment <- experiments[ grepl(paste(metadata[1], "$", sep = ""),
                                         experiments) ]
        exper.path <- list.dirs(experiment, recursive = FALSE)

        imageJ.input <- paste(exper.path, "/Well ", metadata[2], "/",  sep = "")

        if(FALSE %in% file.exists(paste(imageJ.input, dyes, sep = ""))){
          return(paste("not all imagenames provided exist in a wellname directory.",
                       "Inspect it!"))
        }

        imageJ.output <- dirs[i]

        dye.string <- paste(dyes, collapse = ',')

        command <- paste(command.base, ' "', imageJ.input, ";", imageJ.output, ";",
                         dye.string, ";", roi, ";", sep = '')
        system(command)

      } else {
        print("the stack already exists, change overwrite to TRUE")
      }

    }

  } else {
    return("no ROI found, inspect the roi.input")
  }
}

make_image_figure <- function(stack.input,
                                 dyes = list("Alexa 555 Confocal fusion - n000000",
                                             "Alexa 647 - Confocal fusion - n000000",
                                             "DAPI Confocal - n000000",
                                             "Alexa 488 - Confocal_529-24 - n000000"),
                                 color.names = list("green",
                                              "red",
                                              "blue",
                                              "Composite",
                                              "grey",
                                              "Montage"),
                                 stack.name = "Stack.tif",
                                 contrast.max = list(600, 250, 700, 1000),
                                 contrast.min = list(600, 250, 700, 1000),
                                 overwrite = TRUE,
                                 path.to.fiji = "D:/Piotrek/programs/Fiji.app",
                                 scale.bar.mm = 0,
                                 keep.separate.images = FALSE,
                                 scale.bar.on.separate.images = 0,
                              collect.montages.folder = FALSE) {
  # A foo for combining the microscopic images side by side.
  # A stack (from ImageJ) of images is needed.

  # - stack.input; path where the stack.name occurs,
  #   can be higher directory, as the functions look for recursively
  # - dyes; list of image names present in stack, without an extension,
  # - color.names; list of names of colors you want to use.
  #   Allowed are: "red", "green", "blue", "Composite", "any_name_for_montage".
  #   The order of names will indicate the look of the montage.
  #   the order of dyes should correspond to the order of names,
  #   with composite inserted in any place, but not influencing the order. e.g.,
  #   dyes = list("DAPI", "Alexa 555")
  #   names = list("blue", "Composite", "red", "Montage")
  #   is ok. But the "Montage" has to be the last one! And "Montage" can be
  #   any name
  # - stack.name; what is the name of the stack;
  # - contrast.max/contrast.min; values for adjusting the brightness/contrast
  #   The order follows the order of dyes
  # - owerwrite = it's a precaution for unwanted ovewriting
  #   and if you want to add new rois
  # - scale.bar.mm; either 0 or the length in mm.
  #   As of now the scale bar option is available only for pathway images,
  #   so if the image comes from different place,
  #   use 0 and provide the bar manually.
  # - scale.bar.on.separate.images; should the scale bar be on each of separate images?
  #   with keep.sepparate.images = FALSE it make no sense to put this to 1
  # - collect.montages; TRUE if you want to have the montages in one folder
  #   for convenient opening and eye inspection

  dirs <- dirname(list.files(stack.input,
                             recursive = TRUE,
                             full.names = TRUE,
                             patter = stack.name))
  for(input in dirs){
    montage.exists <- file.exists(paste(input, "/Montage.jpg", sep = ""))
    if(montage.exists == FALSE | overwrite == TRUE){
      cell.input <- paste(input, "/", sep="")
      cell.output <- cell.input

      command.base <-
        "ImageJ-win64.exe --console --headless -macro color_figure_general_.ijm"

      setwd(path.to.fiji)
      dye.string <- paste(dyes, collapse = ',')
      min.string <- paste(contrast.min, collapse = ',')
      max.string <- paste(contrast.max, collapse = ',')
      name.string <- paste(color.names, collapse = ',')

      command <- paste(command.base, ' "', cell.input, ";", cell.output, ";",
                       stack.name, ";",
                       min.string, ";", max.string, ";", dye.string, ";",
                       name.string, ";", scale.bar.mm, ";",
                       scale.bar.on.separate.images, ";", sep = '')
      system(command)

      if(keep.separate.images == FALSE){
        for(name in color.names){
          if(name %in% c("Montage")){ next }
          file.remove(list = paste(cell.output, name, ".jpg", sep = ""))
        }
      }
    } else {
      cat(paste("Montage.jpg in \n", input, "\nexists, skipping to next stack"))
    }
  }

  collect_montages <- function(input.path,
                               collecting.folder,
                               montage.name){
    # a small foo for collecting the montages to one folder
    # for convenient opening and eye inspection

    push.dir <- function(folder.name){
      if(!dir.exists(folder.name)){
        dir.create(folder.name, recursive = TRUE)
      }
    }

    montage.dirs <- list.files(input.path,
                                       pattern = montage.name,
                                       recursive = TRUE,
                                       full.names = TRUE)

    for(montage.dir in montage.dirs){

      path.from <- montage.dir

      path.to <- paste(dirname(dirname(montage.dir)), "/",
                       collecting.folder, "/",
                       basename(dirname(montage.dir)), ".jpg", sep = "")
      push.dir(dirname(path.to))

      file.copy(from = path.from,
                to = path.to,
                overwrite = TRUE)
    }
  }
  if(collect.montages.folder != FALSE){

  collect_montages(input.path = stack.input,
                   collecting.folder = collect.montages.folder,
                   montage.name = paste(color.names[[length(color.names)]], ".jpg", sep = ""))
  }
}




