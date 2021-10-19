\name{crop_location}
\alias{crop_location}

\title{
crop the location
}
\description{
\code{crop.location} is a function to crop an image according to ROI file
(region of interest, a file from ImageJ). The function creats the Stack.tif.
  Stack.tif can be used to make the figure using make.image.figure function.
}
\usage{
# the example with minimum number of arguments:

crop_location(roi.input,
dyes = list("DAPI Confocal - n000000",
"Alexa 488 - Confocal_529-24 - n000000"),
path.to.fiji = "D:/Piotrek/programs/Fiji.app",
all.exper.input = "D:/Piotrek/Experiments/PEG/raw_images")
}
\arguments{

  \item{roi.input}{
the path to an existing roi file. All roi files (recursively) in this path will be processed.

}
  \item{all.exper.input}{
the path where the experiments folders are found
}
  \item{pattern.roi}{
regular expression of pattern to find the roi. The deafult follows the convention:
"name_of_the_experiment[]wellname.roi", where "[]" is a separator.
The experiment has to be in all.exper.input folder.
  The name_of_the_eperiment has to explicitly indicate the experiment
  in the all.exper.input
}
  \item{path.to.fiji}{
the path to fiji folder, for example
"D:/Piotrek/programs/Fiji.app"
}
  \item{overwrite}{
should the already created Stacks be overwritten?
  FALSE is usefull when new roi are created
}
  \item{dyes}{
the list of image names taken to create a stack.
}
  \item{stack.name}{
a name of a stack, default: "Stack.tif"
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a tiff file with the \item{stack.name}{} name
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Piotr Topolewski
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
crop_location(roi.input = yourinputpath,
              dyes = list("Alexa 488 - Confocal_529-24 - n000000.tif",
                          "DAPI_PT-nconf cycle - n000000.tif",
                          "Alexa 647 - Confocal - n000000.tif"))
}
