\name{Inoue_IGM}
\alias{tau_IGM_LSLAF}
\alias{tau_IGM_LSDLA}
\alias{tau_IGM_LCLAF}
\alias{tau_IGM_LCDLA}
\alias{tau_IGM_Tot}
\alias{Inoue14_IGM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
IGM absorption model from Inoue (2014). 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Compute the optical depth and attenuation as a function of observed wavelengths using the Inoue (2014) model.
}
\usage{
tau_IGM_LSLAF(l_obs, z)
tau_IGM_LSDLA(l_obs, z)
tau_IGM_LCLAF(l_obs, z)
tau_IGM_LCDLA(l_obs, z)
tau_IGM_Tot(l_obs, z)
Inoue14_IGM(l_obs, z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l_obs}{
%%     ~~Describe \code{x} here~~
  Numeric vector; observed wavelengths in Angstrom.
  }
  \item{z}{
  Numeric scalar; source redshift. 
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The light emitted at the source redshift is attenuated by intervening neutral hydrogen clouds between the source and the observer. Inoue (2014) connects the optical depth of the IGM to the distribution of absorbers as a function of redshift and hydrogen column density. The model includes opacity from Lyman continuum and Lyman series line absorption in the Lyman alpha forest (LAF, log(N_H / cm^-2) < 17.2), Lyman limit systems (LLS, 17.2 < log(N_H / cm^-2) <= 20.3) and damped Lyman systems (DLA, log(N_H / cm^-2) > 20.3). 
Eqs. 21-29 in Inoue (2014). 

The model is only valid for observed wavelengths longer than 911.75 Angstrom. Shorter observed wavelengths are hence not attenuated. 

The observed spectrum will be attenuated as F_obs = F_emit * Inoue14_IGM(l_obs, z). 
}
\value{
For \code{tau_IGM_LSLAF}, \code{tau_IGM_LSDLA}, \code{tau_IGM_LCLAF}, \code{tau_IGM_LCDLA}, \code{tau_IGM_Tot}: 
Vector of optical depth, tau, corresponding to Lyman series/continuum in LAF/LLS/DLA or total for each observed wavelength.

For \code{Inoue_IGM}: 
Vector of exp(-1 * tau_IGM_tot) for each observed wavelength.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Inoue, A. K., Shimizu, I., Iwata, I., & Tanaka, M. 2014, Monthly Notices of the Royal Astronomical Society, 442, 1805
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Aaron Robotham
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{ProSpectSED}}
}
\examples{
l_obs = seq(1000, 50000, 100)
z = 6

plot(
  l_obs, 
  tau_IGM_Tot(l_obs,z), 
  log = "xy", 
  type = "l", 
  xlab = "Observed wavelength [Angstrom]", 
  ylab = "tau"
)
lines(
  l_obs, 
  tau_IGM_LSLAF(l_obs,z), 
  col = "blue"
)
lines(
  l_obs, 
  tau_IGM_LSDLA(l_obs,z), 
  col = "cornflowerblue"
)
lines(
  l_obs, 
  tau_IGM_LCLAF(l_obs,z), 
  col = "red"
)
lines(
  l_obs, 
  tau_IGM_LCDLA(l_obs,z), 
  col = "purple"
)
legend(
  "topleft", 
  lty = 1,
  col = c("black", "blue", "cornflowerblue", "red", "purple"),
  c("Total", 
    "Lyman series LAF", 
    "Lyman series LLS/DLA", 
    "Lyman continuum LAF", 
    "Lyman continuum LLS/DLA")
)

plot(
  l_obs, 
  Inoue14_IGM(l_obs, 6),
  log = "x", 
  type = "l", 
  col = "blue", 
  xlab = "Observed wavelength [Angstrom]",
  ylab = "exp(-1*tau)"
)
lines(
  l_obs, 
  Inoue14_IGM(l_obs, 2), 
  col = "red"
)
legend(
  "topleft", 
  lty = 1,
  col = c("blue", "red"),
  c("z=6", "z=2")
)
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
