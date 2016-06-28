\name{jglfit}
\alias{jglfit}
\docType{data}
\title{
Example data set for generating simulated data sets.
}
\description{
This is the output of applying New-JGL to the gene expression data set used in Gao et al. (2016).
}
\usage{data(jglfit)}
\format{
  A data frame with 0 observations on the following 2 variables.
  \describe{
    \item{\code{x}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list with the following elements:
\item{pie}{A vector of the mixing proportion of each component in the model.} \cr
\item{mu}{A matrix that contains the mean of each component in the model in each row.} \cr
\item{covinv}{A list that contains the precision matrices for each component in the model.} \cr
\item{membership}{The class assignment for each observation.} \cr
\item{par}{The optimal value of lambda_1 and lambda_2.} \cr

}

%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
# load data set

data(jglfit)

# mixing proportion
jglfit$pie

# dimension of an estimated precision matrix
dim(jglfit$covinv[[1]])

}
\keyword{datasets}
