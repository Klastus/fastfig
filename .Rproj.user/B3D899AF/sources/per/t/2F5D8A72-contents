package.list <- list("ggplot2", 
                     "gridExtra",
                     "data.table",
                     "reshape2",
                     "deamer",
                     "foreach",
                     "doParallel",
                     "dplyr",
                     "SysBioSigTheme",
                     "scales",
                     "RColorBrewer",
                     "ggforce",
                     "rlist",
                     "extrafont",
                     "flowCore",
                     "cowplot")
try({package.list
  package.load <- sapply(package.list, function(package.name){
    package.exist <- require(package.name, character.only = TRUE)
    if(!package.exist){
      install.packages(package.name)
      return(library(package.name, character.only = TRUE))
    }
    return(package.exist)
  })
})

normalize_data <- function(data,
                           normalize_factor = 65535){
  data.intensity_colnames <- grepl("Intensity", colnames(data)) & 
    !grepl("Location", colnames(data))
  data[, data.intensity_colnames] <- data[, data.intensity_colnames] * normalize_factor
  return(list(data = data))
}


cc.assign <- function(dane, cc.borders, DAPI.column){
  c.phase <- c()
  for (i in dane[[DAPI.column]]) {
    if (i < cc.borders[1]){
      c.phase <- c(c.phase, "outliers")
    } else if (i < cc.borders[2]){
      c.phase <- c(c.phase, "G1")
    } else if (i < cc.borders[3]){
      c.phase <- c(c.phase, "S")
    } else if (i < cc.borders[4]){
      c.phase <- c(c.phase, "G2/M")
    } else { 
      c.phase <- c(c.phase, "outliers")
    }
  }
  return(c.phase)
}

discriminate <- function(data, borders, var.column, factor.list){
  c.phase <- c()
  for (i in data[[var.column]]) {
    if (i <= borders[[1]]){
      c.phase <- c(c.phase, factor.list[[1]])
    } else if (i <= borders[[2]]){
      c.phase <- c(c.phase, factor.list[[2]])
    } else if (i <= borders[[3]]){
      c.phase <- c(c.phase, factor.list[[3]])
    } else if (i <= borders[[4]]){
      c.phase <- c(c.phase, factor.list[[4]])
    } else { 
      c.phase <- c(c.phase, factor.list[[1]])
    }
  }
  return(factor(c.phase, levels = paste(factor.list)))
}

distinguish <- function(data, 
                        borders, 
                        var.column, 
                        factor.list = list("outliers", "G1", "S", "G2/M")){
  c.phase <- c()
  for (i in data[, colnames(data) %in% var.column]) {
    for(discriminated in 1:length(borders)){
      if (i <= borders[[discriminated]]){
        c.phase <- c(c.phase, factor.list[[discriminated]])
        break
      } else if(discriminated == length(borders)){
        c.phase <- c(c.phase, factor.list[[1]])
      }
    }
  }
  return(factor(c.phase, levels = paste(factor.list)))
}

distinguish.cc <- function(data = a, 
                           borders, 
                           var.column, 
                           factor.list = list("outliers", "G1", "S", "G2/M")){
  
  c.phase <- ifelse(data[[var.column]] <= borders[[1]], factor.list[[1]],
                    ifelse(data[[var.column]] <= borders[[2]], factor.list[[2]], 
                           ifelse(data[[var.column]] <= borders[[3]], factor.list[[3]],
                                  ifelse(data[[var.column]] <= borders[[4]], factor.list[[4]],
                                         factor.list[[1]]))))
  
  return(factor(c.phase, levels = paste(factor.list)))
}

change.dapi.name <- function(x){
  dapi.name <- grep("DAPI", colnames(x), value = TRUE)
  x <- x %>%
    dplyr::rename(Integrated_DAPI_epi = all_of(dapi.name))
  return(x)
}

dna.histogram <- function(data, channel, bins = 100, 
                          border.list, title = "random title", 
                          xlab = "DAPI", ylab = "count", xlim,
                          thickness = 1, 
                          aspect.ratio = 1,
                          linetype = 2){
  g <- ggplot(data, aes_string(x=channel))+
    geom_histogram(bins=bins)+
    geom_vline(xintercept = border.list[[1]], col="red", 
               size = thickness, linetype = linetype)+
    geom_vline(xintercept = border.list[[2]], col="orange", 
               size = thickness, linetype = linetype)+
    geom_vline(xintercept = border.list[[3]], col="green", 
               size = thickness, linetype = linetype)+
    geom_vline(xintercept = border.list[[4]], col="blue", 
               size = thickness, linetype = linetype)+
    ggtitle(title)+
    scale_x_continuous(expand = c(0, 0), name = xlab, limits = xlim)+
    scale_y_continuous(expand = c(0, 0), name = ylab)+
    theme_trajectories(aspect.ratio = aspect.ratio)+
    ylab(ylab)
  return(g)
}

cc.to.factor <- function(phases=c("G1","S","G2/M"), df){
  df.subset <- df[df$phase %in% phases, ]
  df.subset$phase <- factor(df.subset$phase, levels = phases)
  return(df.subset)
}

# variable.subset <- function(data, columns, new.columns = 0){
#   data.2 <- data[, colnames(data) %in% columns]
#   if(new.columns[1]!=0){
#     colnames(data.2) <- new.columns
#   }
#   return(data.2)
# }

variable.subset <- function(data, 
                            columns,
                            new.columns = 0){
  data <- data %>%
    select(all_of(columns))
  if(new.columns[1]!=0){
    colnames(data) <- new.columns
  }
  return(data)
}


randomRows <- function(df, n){
  return(df[sample(nrow(df), n), ])
}


CV <- function(x){
  return(sd(x)/mean(x))
}

push.dir <- function(folder.name){
  if(!dir.exists(folder.name)){
    dir.create(folder.name, recursive = TRUE)
  }
}

zero.time.multiply <- function(data.zero, data.no.zero, 
                               time.of.zero = 30,
                               stimulation.of.zero = 0){
  data.zero.unchanged <- data.zero
  data.zero.unchanged$duplication <- "no"
  data.zero.unchanged$time.1.1 <- time.of.zero
  data.zero.unchanged$stimulation.1.1 <- stimulation.of.zero
  data.zeros <- data.frame()
  
  stims <- unique(data.no.zero$stimulation.1.1)
  
  for(stim in stims){
    data.zero$stimulation.1.1 <- stim
    data.zero$time.1.1 <- 0
    data.zeros <- rbind(data.zeros, data.zero)
  }
  data.zeros$duplication <- "yes"
  data.no.zero$duplication <- "no"
  
  if(stimulation.of.zero == -1){
    return(rbind(data.no.zero, data.zeros))
  } else {
    return(rbind(data.zero.unchanged, data.no.zero, data.zeros))
  }
}



zero.multiply <- function(data.zero, 
                          data.no.zero, 
                          time.of.zero = 30){
  # data.zero.unchanged <- data.zero
  # data.zero.unchanged$duplication <- "no"
  # data.zero.unchanged$time.1.1 <- time.of.zero
  # data.zero.unchanged$stimulation.1.1 <- stimulation.of.zero
  data.zeros <- data.frame()
  
  stims <- unique(data.no.zero$stimulation.1.1)
  
  for(stim in stims){
    data.zero$stimulation.1.1 <- stim
    data.zero$time.1.1 <- 0
    data.zeros <- rbind(data.zeros, data.zero)
  }
  
  for(time.point in time.of.zero){
    data.zero$stimulation.1.1 <- 0
    data.zero$time.1.1 <- time.point
    data.zeros <- rbind(data.zeros, data.zero)
  }
  
  data.zeros$duplication <- "yes"
  data.no.zero$duplication <- "no"
  
  return(rbind(data.no.zero, data.zeros))
  
}




copy.IPIQA <- function(ipiqa.path = "D:/IPIQA/analysis_output/2019-08-27/",
                       core.destination = "D:/Piotrek/Experiments/ICF/",
                       experiment.full.name,
                       csv.name.list = c("ShrinkedNuclei.csv"),
                       project,
                       shift.value = FALSE) {
  
  string <- tail(strsplit(experiment.full.name, "-")[[1]], 1)
  exper.last.name <- regmatches(string, regexpr("(PT|JW)\\d+", string))
  
  exper.path = list.dirs(paste(ipiqa.path, experiment.full.name, sep = ""),
                         recursive = FALSE, full.names = TRUE) [1]
  
  normalizations <- list.dirs(exper.path, recursive = FALSE, full.names = FALSE)
  for(normaliz in normalizations){
    for(csv.name in csv.name.list){
      csv.path <- paste(exper.path, 
                        normaliz, 
                        "data_quantify", 
                        csv.name,
                        sep = "/")
      destination.path <- paste(core.destination, 
                                exper.last.name, 
                                "/R/", 
                                project, 
                                "/input/", 
                                normaliz, 
                                "/",
                                if(shift.value != FALSE){
                                  shift.value},
                                sep = "")
      push.dir(destination.path)
      file.copy(from = csv.path,
                to = destination.path,
                overwrite = TRUE)
    }
  }
}

theme_trajectories <- function(border.thickness = 0.5,
                               axis.num.size = 10, 
                               axis.name.size = 11,
                               aspect.ratio = FALSE){
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size=border.thickness, fill=NA),
        panel.background = element_blank(),
        axis.ticks = element_line(colour = "black", size = border.thickness),
        axis.text.x = element_text(colour = "black", size = axis.num.size),
        axis.text.y = element_text(colour = "black", size = axis.num.size),
        strip.background = element_blank(),
        axis.title = element_text(face="plain", size = axis.name.size),
        plot.title = element_text(hjust = 0.5),
        legend.key=element_blank())+
    if(aspect.ratio != FALSE){
      theme(aspect.ratio = aspect.ratio)
    } else {theme()}
}

add_letter <- function(letter = "A",
                       size = 20, 
                       face = "plain", 
                       color = "black",
                       family = "sans"){
  list(labs(tag = letter),
       theme(plot.tag = element_text(size = size, 
                                     face = face, 
                                     color = color,
                                     family = family)))
} 

theme.itrc <- function(){
  theme_itrc <-theme(axis.title = element_text(size= 11,
                                               face="plain",
                                               vjust = 0.5),
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     axis.line.x = element_blank(),
                     axis.line.y = element_blank(),
                     panel.background = element_blank(),
                     strip.background = element_blank(),
                     strip.text =
                       element_text(size= 10,
                                    face="plain",
                                    vjust = 0.5#,lineheight = theme.text_size*3
                       ))
}

theme_thesis <- function(border.thickness = 0.5,
                         axis.num.size = 10, 
                         axis.name.size = 11,
                         aspect.ratio = 1,
                         family = "Poppins"){
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size=border.thickness, fill=NA),
        panel.background = element_blank(),
        axis.ticks = element_line(colour = "black", size = border.thickness),
        axis.text.x = element_text(colour = "black", size = axis.num.size),
        axis.text.y = element_text(colour = "black", size = axis.num.size),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 12),
        axis.title = element_text(face="plain", size = axis.name.size),
        plot.title = element_text(hjust = 0.5),
        legend.key=element_blank(),
        legend.background = element_blank(),
        text = element_text(family = family))+
    if(aspect.ratio != FALSE){
      theme(aspect.ratio = aspect.ratio)
    } else {theme()}
}

nice.colors <- c("forestgreen", "steelblue", "coral3", "cadetblue", "mediumorchid4")
cc.pal <- c("coral3", "darkorange3", "skyblue3")
# IFN.pal <- c(brewer.pal(5, "PuBu"))[c(1,2,3,5)]
IFN.pal <- c(brewer.pal(8, "PuBu"))[c(3,5,7,8)]
TNF.pal <- c(brewer.pal(5, "OrRd"))[c(1,2,3,5)]
OSM.pal <- c("#cbdeda", "#63bda3", "#218845", "#114320")
overlaps.pal <- c(brewer.pal(5, "BuGn"))[c(1,2,3,5)]

pdf.push <- function(file, ...){
  if(!dir.exists(dirname(file))){
    dir.create(dirname(file), recursive = TRUE)
  }
  pdf(file, ...)
}

extract.MK <- function(exp.name){ 
  base.path <- "D:/Piotrek/Experiments/ICF/"
  rds.path <- paste(base.path, "MK_all/R/input/data_ffc_filtered_oranized.RDS", sep = "")
  poster.data.list <- readRDS(file = rds.path)
  
  F1=poster.data.list[[exp.name]]
  
  F1$time <- factor(F1$time)
  F1$Intensity_MeanIntensity_Alexa488 <- F1$Intensity_MeanIntensity_Alexa
  F1$stimulation <- as.factor(F1$stimulation)
  return(F1)
}

pbs.reference <- 19.6605
pbs.normalize <- function(data, 
                          pbs.camcor, 
                          pbs.reference = 19.6605){
  data.intensity_colnames <- grepl("^Intensity_[A-z]*", colnames(data))
  data[,data.intensity_colnames] <- data[,data.intensity_colnames]*pbs.reference/pbs.camcor
  return(data)
}

load.bio <- function(bio.path,
                     columns = c("Intensity_IntegratedIntensity_Alexa488",
                                 "Intensity_MeanIntensity_Alexa488")){
  # a foo for loading the experiment .csv.
  
  
  bio.nuc <- read.table(bio.path, header=TRUE, sep=",")
  
  bio.nuc <- normalize_data(bio.nuc)$data
  
  from <- which(colnames(bio.nuc) == grep("\\.1\\.1|\\.1\\.0", colnames(bio.nuc), value = TRUE)[1])
  to <- length(bio.nuc[1 ,])
  
  columns.basic <- c("AreaShape_Area",
                     "AreaShape_Center_X",
                     "AreaShape_Center_Y")
  
  if(columns[1] == "all"){
    columns.all <- colnames(bio.nuc)
  } else {
    columns.all <- c(columns.basic, columns, 
                     colnames(bio.nuc)[from:to])
  }
  bio.nuc <- bio.nuc[, colnames(bio.nuc) %in% columns.all]
  
  return(bio.nuc)
}

load_fcs <- function(file, sample.type){
  fcs <- read.FCS(filename = file,
                  alter.names = TRUE)
  epp <- data.frame(exprs(fcs))
  epp$sample_type <- sample.type
  return(epp)
}

options("max.print" = 100)

add.cls.well.name.fast <- function(data){
  library(dplyr)
  data.well <- data.frame("Row" = rep(1:8, each = 12),
                          "Column" = rep(1:12, times = 8))
  data.well$well.name <- paste(LETTERS[data.well$Row], 
                               formatC(data.well$Column,
                                       width = 2,
                                       flag = "0"), 
                               sep = "")
  
  return(data %>%
           left_join(data.well, by = c("Row", "Column")))
}

extract.bck <- function(path.to.bck.csv,
                        wells,
                        bck.column = "Intensity_MeanIntensity_bck"){
  # a foo to extract the single value of background from bck csv file, 
  #   e.g., Image.csv of cellprofiler#
  # you need to load normalize_data.R first!
  
  # - path.to.bck.csv: a path to csv file contating bck values, 
  #   normally from Image.csv file 
  # - wells: vector of wells from which bck value will be extracted
  # - bck.column: should be in .csv file
  if(!file.exists(path.to.bck.csv)){return("No bsk file found")}
  
  bck <- normalize_data(read.table(path.to.bck.csv, 
                                   header=TRUE, sep=","))$data
  bck <- bck[bck$well.name %in% wells, ]
  # bck.column %in% colnames(bck)
  mean.bck <- mean(bck[[bck.column]])
  return(mean.bck)
}


normalize.to.global.avg <- function(df,
                                    min.bck = 0,
                                    avg.zero.bck.subtraction = TRUE,
                                    avg.zero.push = FALSE,
                                    fluorescence.colname,
                                    zero.wells,
                                    remove.negative.cells = TRUE){
  ## a foo to make a normalized fluorescnce column ##
  
  # df, fluo.column, attributes of mean wells, min. value (with false if you want to take minimum)
  
  # - min.bck: the value of bck to be subtracted. 
  #   If FALSE, the minimum value of the fluorescence.colname in df is taken
  # - avg.zer.push: a number which you want to provide as the mean value. 
  #   FALSE, if you want to use the wells. 
  #   Note: mean.zero.push has the priority over the "zero.wells" argument
  #   Note2: if a number, it should already be bck-subtracted!
  # - avg.zero.bck.subtraction: TRUE if the bck should be subtracted 
  #   from the provided mean.zer.push value
  # - zero.wells: reference wells, considered as 1 value of fluorescence. Should 
  #   be identified outside of this function 
  # - fluorescence.colname: the column name to be normalized
  
  if(!fluorescence.colname %in% colnames(df)){
    return("Fluorescence colname not in df!")
  }
  min.zero <- min.bck
  if (min.bck == FALSE){
    min.zero <- min(df[[fluorescence.colname]])
  } 
  
  if(avg.zero.push == FALSE){
    avg.zero <- 
      mean(df[df$well.name %in% zero.wells, ][[fluorescence.colname]]) - min.zero
  } else {
    if(avg.zero.bck.subtraction == TRUE){
      avg.zero <- avg.zero.push - avg.zero
      print("Information: avg.zero.push is bck-subtracted")
    } else if(avg.zero.bck.subtraction == FALSE){
      avg.zero <- avg.zero.push
      print("Information: avg.zero.push is not bck-subtracted")
    }
  }
  
  df[[paste(fluorescence.colname, "_fold", sep = "")]] <- 
    (df[[fluorescence.colname]] - min.zero)/avg.zero
  if(remove.negative.cells == TRUE){
    df <- df[df[[paste(fluorescence.colname, "_fold", sep = "")]] > 0, ]
  }
  return(df)
}

# how to install and use fonts:
# 1. download the font from any webpage (ttf format) and install it (right click on windows)
# 2. load font in R: windowsFonts("LM Roman 10"=windowsFont("LM Roman 10"))
# 3. e.g. family = "LM Roman 10" in theme() or geom_xxx
# 4. save cairo_pdf()
# in case of troubles try the package extrafont

# on mac: run loadfonts() before
# loadfonts()
# font_import(pattern = "08635755*") 
# but it doesn't work anyway


# second way of font installing:
# 1. download the font from any webpage (ttf format) and install it for all users (right click on windows)
 # It is crucial to have the font installed in C:Windows/Fonts, instead of C:Users.... !!
# 2. extrafont::font_import()
# 3. extrafont::loadfonts(device="win")

cbind.compartments <- function(path,
                               compartment.files,
                               fluo.columns.to.stay,
                               new.column.names,
                               first.columns.to.stay = 0){
  # normalize_data() is needed!
  # library(dplyr)
  # this is function to Pathway data only!
  # - path: path to folder with .csv files
  # - compartment.files: name of csv filet to be merged, with extention!
  # - fluo.columns.to.stay: character vector of fluorescent columns to stay. 
  #   Note: area, x,y and other first columns can stay 
  #   according to first.columns.to.stay
  # - new.column.names: should follow the order of fluo.columns.to.stay
  # - first.columns.to.stay: this number of first columns stay. Usually their are x,y, area etc,  
  
  df <- list()
  possible.compartments <- list(list("ucl",
                                     "cell",
                                     "cyto|plasm"),
                                list("nucleus",
                                     "cell",
                                     "cytoplasm"))
  for(compartment.file in compartment.files){
    for(compartment.guess in 1:length(possible.compartments[[1]])){
      comp.matched <- grepl(pattern = possible.compartments[[1]][[compartment.guess]], 
                            x = compartment.file, 
                            ignore.case = TRUE)
      if(comp.matched){
        df[[possible.compartments[[2]][[compartment.guess]]]] <- 
          normalize_data(read.table(paste(path, compartment.file, sep=''), 
                                    header=TRUE, sep=","))$data
        break 
      }
      if(compartment.guess == 3){
        return(paste(compartment.file, "not found!"))
      }
    }
  }
  
  
  df.for.columns <- df[[1]] 
  for(name in names(df)){
    
    df[[name]] <- variable.subset(df[[name]], 
                                  columns = fluo.columns.to.stay,
                                  new.columns = paste0(new.column.names, "_",
                                                       name))
    
  }
  # compartment.files <- c("CellsFiltered647.csv",
  #                        "ShrinkedNucleiMasked.csv")
  
  
  from <- which(colnames(df.for.columns) == grep("\\.1\\.1|\\.1\\.0", colnames(df.for.columns), 
                                                 value = TRUE)[1])
  to <- length(df.for.columns[1 ,])
  
  
  
  return(cbind(df.for.columns[, head(colnames(df.for.columns), first.columns.to.stay)],
               dplyr::bind_cols(df),
               df.for.columns[, colnames(df.for.columns)[from : to]]))
}