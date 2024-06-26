\name{ProSpect-package}
\alias{ProSpect-package}
\alias{ProSpect}
\docType{package}
\title{
\packageTitle{ProSpect}
}
\description{
\packageDescription{ProSpect}
}
\details{
\tabular{ll}{
Package: \tab ProSpect\cr
Type: \tab Package\cr
Version: \tab 1.4.1\cr
Date: \tab 2024-02-26\cr
License: \tab LGPL-3\cr
Depends: \tab R (>= 3.0)\cr
Imports:  \tab celestial, NISTunits, pracma\cr
Suggests: \tab shiny, magicaxis, knitr, LaplacesDemon, foreach, data.table, hyper.fit, doParallel\cr
}
}

\author{
\packageAuthor{ProSpect}

Maintainer: \packageMaintainer{ProSpect}
}
\references{
Akhlaghi et al, 2015, ApJSS, 220, 1 \cr
Bertelli et al., 1994, A&AS, 106, 275 \cr
Brammer, van Dokkum & Coppi, 2008, ApJ, 686, 1503 \cr
Bruzual and Charlot, 2003, MNRAS, 344, 1000 \cr
Chabrier, 2003, PASP, 115, 763 \cr
Condon & Yin, 1990. ApJ, 357, 97 \cr
da Cunha et al, 2008, MNRAS, 388, 1595 \cr
Dale et al., 2014, ApJ, 784, 11 \cr
Draine et al, 2007, ApJ, 663, 866 \cr
Feltre et al, 2012, MNRAS, 426, 120 \cr
Fritz et al, 2006, MNRAS, 366, 767 \cr
Kennicutt et al, 1994, ApJ, 435, 22 \cr
Kennicutt et al, 1998, ARAA, 36, 189 \cr
Kroupa, 2001, MNRAS, 322, 231 \cr
Levesque, Kewley & Larson, 2010, AJ, 139, 712 \cr
Marvil et al, 2015, AJ, 149, 32 \cr
Robotham et al, 2020, MNRAS, 495, 905 \cr
Songaila, 2004, AJ, 127, 2598 \cr
Stanway and Eldridge, 2018, MNRAS, 479, 75 \cr
Vazdekis et al., 2016, MNRAS, 463, 3409
}
\examples{
SFHdemo_dust=SFHfunc()
plot(SFHdemo_dust$out[,2:3], log='x', ylim=c(26,17), xlab=BC03lr$Labels$Wavelab,
ylab='Mag', col=rev(rainbow(26, end=2/3)), pch=16)
}
\concept{SED}
