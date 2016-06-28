\name{pGMM}
\alias{pGMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating multiple Gaussian graphical models (GGM) in penalized Gaussian mixture models (GMM).%%  ~~function to do ... ~~
}
\description{
Fit the Gaussian mixture model using either New-SP or New-JGL method. Sparsity and fusion penalty are imposed on the component-specific precision matrices with either non-convex (New-SP) or convex (New-JGL) function. The output contains the class assignment, the mixing proportion, the means and the component-specific precision matrices and the optimal tuning parameter (if a vector of tuning parameter is supplied).
}
\usage{
pGMM(Y, k, method, lambda1, lambda2, ncluster, tau=0.01, threshold=1e-5, 
  MAX_iter=100, seed=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  The input data as an \code{n} by \code{p} matrix, where \code{n} is the number of observations and \code{p} is the number of variables.
%%     ~~Describe \code{Y} here~~
}
  \item{k}{
  The fold for cross-validation. Default is 5.
%%     ~~Describe \code{k} here~~
}
  \item{method}{
Either "New-SP" or "New-JGL". "New-SP" uses the non-convex truncated lasso penalty (TLP). "New-JGL" uses the a convex lasso penalty. Default is "New-SP".
}
  \item{lambda1}{
  Tuning parameter for the sparsity penalty, either a single value or a vector. If the input is a vector, cross-validation is applied to select the optimal value.
%%     ~~Describe \code{lambda1} here~~
}
  \item{lambda2}{
  Tuning parameter for the fusion penalty, either a single value or a vector. If the input is a vector, cross-validation is applied to select the optimal value.
%%     ~~Describe \code{lambda2} here~~
}
  \item{ncluster}{
  Number of components in the Gaussian mixture.
%%     ~~Describe \code{ncluster} here~~
}
  \item{tau}{
  The tuning paramter in the non-convex penalty function. Default is 0.01.
%%     ~~Describe \code{tau} here~~
}
  \item{threshold}{
  Threshold for convergence. Default is 1e-5.
%%     ~~Describe \code{threshold} here~~
}
  \item{MAX_iter}{
  Maximum number of iterations. Default is 100.
%%     ~~Describe \code{MAX_iter} here~~
}
  \item{seed}{
  The seed used when splitting data for cross-validation. Default is 1.
%%     ~~Describe \code{seed} here~~
}
}
\details{
It aims to estimate multiple networks in the presence of sample
heterogeneity, where the independent samples (i.e. observations) may come
from different and unknown populations or distributions. Specifically, we
consider penalized estimation of multiple precision matrices in the framework
of a Gaussian mixture model. A major innovation is to take advantage
of the commonalities across the multiple precision matrices through possibly
nonconvex fusion regularization, which for example makes it possible
to achieve simultaneous discovery of unknown disease subtypes and detection
of differential gene (dys)regulations in functional genomics. We embed
in the EM algorithm one of two recently proposed methods for estimating
multiple precision matrices in Gaussian graphical models: one is based
on the convex Lasso penalty (Danaher et al., 2014), and the other is based
on the con-convex TLP penalty (Zhu et al., 2014); the former is faster,
while the latter is less biased and will perform better with a relatively larger sample size.

The majority of the existing literature on mixture modeling focus on regularizing only the mean parameters with diagonal covariance matrices, though some have started considering regularization of the covariance parameters too, all of which, however, do not touch on the key issue of identifying both common and varying substructures of the precision matrices across the components of a mixture model. Since these methods always give different networks for different components unless a common network is assumed, they do not address the question of interest here: which parts of the networks change with the components.
Since a fusion penalty is used to  shrink multiple networks towards
each other, the proposed methods not only are statistically more efficient with information borrowing, but also facilitate interpretation in identifying differential network
substructures. In particular, due to the use of a non-convex penalty, the adopted method
of Zhu et al. (2014) strives to uncover the commonalities among multiple
networks while maintaining their unique substructures too.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

\item{pie}{A vector of the mixing proportion of each component in the model.} \cr
\item{mu}{A matrix that contains the mean of each component in the model in each row.} \cr
\item{covinv}{A list that contains the precision matrices for each component in the model.} \cr
\item{membership}{The class assignment for each observation.} \cr
\item{par}{The optimal value of lambda_1 and lambda_2.} \cr

}
\references{
% ~~ Literature or other references for background information ~~
Gao, C., hu, Y., Shen, X., and Pan, W. (2016). Estimation of multiple networks in Gaussian mixture models, \emph{Electronic Journal of Statistics}, \bold{10}(1), 1133--1154. \cr
\url{http://projecteuclid.org/download/pdfview_1/euclid.ejs/1462192266} \cr
Zhu, Y., Shen, X., and Pan, W. (2014). Structural pursuit over multiple
undirected graphs, \emph{Journal of the American Statistical Association}, \bold{109}(508), 1683--1696. \cr
\url{http://www.tandfonline.com/doi/pdf/10.1080/01621459.2014.921182} \cr
Danaher, P., Wang, P., and Witten, D. M. (2014). The joint graphical lasso
for inverse covariance estimation across multiple classes, \emph{Journal of the Royal
Statistical Society, Series B}, \bold{76}(2), 373--397. \cr
\url{http://onlinelibrary.wiley.com/doi/10.1111/rssb.12033/epdf} \cr
}
\author{
Chen Gao, Yunzhang Zhu
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

data(jglfit)

# simulate data

set.seed(12345)



Y1 = rmvnorm(1*sum(jglfit$membership==1), jglfit$mu[1,], solve(round(jglfit$covinv[[1]], digits=3)))
Y2 = rmvnorm(1*sum(jglfit$membership==2), jglfit$mu[2,], solve(round(jglfit$covinv[[2]], digits=3)))
Y3 = rmvnorm(1*sum(jglfit$membership==3), jglfit$mu[3,], solve(round(jglfit$covinv[[3]], digits=3)))
Y4 = rmvnorm(1*sum(jglfit$membership==4), jglfit$mu[4,], solve(round(jglfit$covinv[[4]], digits=3)))

Y=rbind(Y1,Y2,Y3,Y4)
p=ncol(Y)

# test

test = pGMM(Y, k=2, method="New-SP", lambda1=1, lambda2=1, ncluster=4)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
