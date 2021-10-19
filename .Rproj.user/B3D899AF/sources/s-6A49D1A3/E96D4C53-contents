\name{make_image_figure}
\alias{make_image_figure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
put images next to the other creating a nice publish-ready figure

}
\description{
A function for combining the microscopic images side by side.
A stack (from ImageJ) of images is needed.
\cr
\cr
Important! Please, place each Stack.tiff file into a separate folder
(eg. area_1/Stack.tiff; area_2/Stack.tiff etc.).
This would enable collecting the montage.jpg files into one folder
(see collect.montages.folder argument)
}
\usage{
make_image_figure(stack.input,
dyes = list("Alexa 555 Confocal fusion - n000000", "Alexa 647 - Confocal fusion - n000000",),
names = list("red", "yellow", "Composite", "Montage"),
stack.name = "Stack.tif",
contrast.max = list(600, 250),
contrast.min = list(0, 0),
path.to.fiji = "D:/Piotrek/programs/Fiji.app",
scale.bar.mm = 0.02)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stack.input}{
path where the stack.name occurs,
  can be higher directory, as the functions look for recursively
}
  \item{dyes}{
a list of image names present in stack, without an extension
}
  \item{color.names}{
  a list of names of colors you want to use.
    Allowed are: "red", "green", "blue", "Composite", "any_name_for_montage".
    The order of names will indicate the look of the montage.
    the order of dyes should correspond to the order of names,
    with composite inserted in any place, but not influencing the order. e.g.,
    dyes = list("DAPI", "Alexa 555")
    names = list("blue", "Composite", "red", "Montage")
    is ok. But the "Montage" has to be the last one! And "Montage" can be
    any name
}
  \item{stack.name}{
a name of a stack, default: "Stack.tif"
}
  \item{contrast.max}{
values for adjusting the brightness/contrast
The order follows the order of dyes
}
  \item{contrast.min}{
values for adjusting the brightness/contrast
The order follows the order of dyes
}
  \item{overwrite}{
it's a precaution for unwanted ovewriting
and if you want to add new rois
}
  \item{path.to.fiji}{
the path to fiji folder, for example
"D:/Piotrek/programs/Fiji.app"
}
  \item{scale.bar.mm}{
either 0 or the length in mm.
As of now the scale bar option is available only for pathway images,
so if the image comes from different place,
use 0 and provide the bar manually.
}
  \item{keep.separate.images}{
during the function run, each image is saved as a separate .jpg file.
As a default, those images are removed after the montage is created.
You can keep them changeing this argument to TRUE
}
  \item{scale.bar.on.separate.images}{
should the scale bar be on each of separate images? Default is 0.
with keep.sepparate.images = FALSE it make no sense to put this to 1
}
\item{collect.montages.folder}{
FALSE, the default, does nothing. Put any character string (e.g. "collected")
if you want to have the montages in that new folder
for convenient opening and eye inspection
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
