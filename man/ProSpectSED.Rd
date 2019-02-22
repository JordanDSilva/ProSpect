\name{ProSpectSED}
\alias{ProSpectSED}
\alias{ProSpectSEDlike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Top Level ProSpect SED Functions
}
\description{
\code{ProSpectSED} computes a fully attenuated and re-emitted SED model incompassing star light (for a number of different star formation history models), dust and AGN. \code{ProSpectSEDlike} is an interface for computing likelihoods, and by extension is the workhorse code that allows model fitting using a number of techniques.
}
\usage{
ProSpectSED(SFH, z = 0.1, tau_birth = 1, tau_screen = 0.3, tau_AGN = 1,
pow_birth = -0.7, pow_screen = -0.7, pow_AGN = -0.7, alpha_SF_birth = 1,
alpha_SF_screen = 3, alpha_SF_AGN = 0, AGNlum = 0, sparse = 5, speclib = NULL,
Dale = NULL, AGN = NULL, filtout = NULL, Dale_M2L_func = NULL, returnall = TRUE,
H0 = 67.8, OmegaM = 0.308, OmegaL = 1-OmegaM, waveout = seq(2,7,by=0.01), ref, ...)
ProSpectSEDlike(parm = c(8,9,10,10,0,-0.5,0.2), Data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SFH}{
Function; the star formation history to use. Must be one of \code{\link{SFHp4}}, \code{\link{SFHp5}} and \code{\link{SFHfunc}}.
}
  \item{z}{
Numeric scalar; cosmological redshift of the source. If the provided \option{z}=0 then the returned fluxes will be absolute, i.e. absolute magnitudes if \option{outtype}='mag'. If \option{z}<0 then only luminosities are computed (flux parts are skipped).
}
  \item{tau_birth}{
Numeric scalar; dust tau of birth clouds. Associated with sub 10 Myr star formation. See \code{\link{CF_birth}}.
}
  \item{tau_screen}{
Numeric scalar; dust tau of the dust screen. See \code{\link{CF_screen}}.
}
  \item{tau_AGN}{
Numeric scalar; dust tau of the dust screen. See \code{\link{CF}}.
}
  \item{pow_birth}{
Numeric scalar; power to further raise the dust attenuation of birth clouds. See \code{\link{CF_birth}}.
}
  \item{pow_screen}{
Numeric scalar; power to further raise the dust attenuation of the dust screen. See \code{\link{CF_screen}}.
}
  \item{pow_AGN}{
Numeric scalar; power to further raise the dust attenuation of the AGN taurus. See \code{\link{CF_screen}}.
}
  \item{alpha_SF_birth}{
Numeric scalar; alpha slope of the birth cloud dust. Lower values mean hotter dust. See \code{\link{Dale}}.
}
  \item{alpha_SF_screen}{
Numeric scalar; alpha slope of the screen dust. Lower values mean hotter dust. See \code{\link{Dale}}.
}
  \item{alpha_SF_AGN}{
Numeric scalar; alpha slope of the AGN taurus dust. Lower values mean hotter dust. See \code{\link{Dale}}.
}
  \item{AGNlum}{
%%     ~~Describe \code{AGNlum} here~~
}
  \item{sparse}{
%%     ~~Describe \code{sparse} here~~
}
  \item{speclib}{
Object; required. Pass in the spectral library directly. Must be one of \code{\link{BC03lr}}, \code{\link{BC03hr}}, \code{\link{EMILES}}.
}
  \item{Dale}{
Object; required. Pass in the Dale dust library \code{\link{Dale_NormTot}} directly (otherwise will run, but not work correctly).
}
  \item{AGN}{
Object; required. Pass in the AGN template to use. Must be one of \code{\link{AGN_UnOb}} or \code{\link{AGN_UnOb_Sparse}}.
}
  \item{filtout}{
Object; required. Pass in the photometric filters directly. If \code{filters} is a vector of filter names, this can be achieved efficiently with a command like: filtout=foreach(i = filters)\%do\%{getfilt(i)} (see \code{\link{getfilt}}).
}
  \item{Dale_M2L_func}{
Function; if dust masses are desired then the user should supply a function to convert \option{alpha_SF} to mass-to-light, e.g. \code{\link{Dale_M2L_func}}
}
  \item{returnall}{
Logical; if true then returns all the various sub SEDs. See Details. If FALSE then code just returns \option{photom_out} (photometry in Jansky).
}
  \item{H0}{
Numeric scalar; Hubble constant as defined at z=0 (default is H0=67.8 (km/s)/Mpc).
}
  \item{OmegaM}{
Numeric scalar; Omega Matter today (default is 0.308).
}
  \item{OmegaL}{
Numeric scalar; Omega Lambda today (default is for a flat Universe with OmegaL = 1 - OmegaM = 0.692).
}
  \item{waveout}{
Numeric vector; desired output log10 wavelength grid to use in Ang. Default covers typical range of galaxy SED data.
}
  \item{ref}{
The name of a reference cosmology to use, one of 137 / 737 / Planck / Planck13 / Planck15 / Planck18 / WMAP / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. Planck = Planck18 and WMAP = WMAP9. The usage is case insensitive, so wmap9 is an allowed input. This overrides any other settings for H0, OmegaM and OmegaL. If OmegaR is missing from the reference set then it is inherited from the function input (0 by default). See \code{\link{cosref}} for details.
}
  \item{parm}{
Numeric vector; arguments to be passed to \code{ProSpectSED}.
}
  \item{Data}{
Object; data object containing everything else needs. See Details.  
}
  \item{\dots}{
All further arguments are passed directly onto \option{SFH}, so by this mechanism you can set function specific arguments that control the star formation history etc.
}
}
\details{
This function allows full modelling of quite sophisticated star formation histories, dust attentuation and re-emission, as well as correctly attenuated AGN. It is useful directly, but it is also the workhorse function for the \code{\link{ProSpectSEDlike}} function, that allows us to compute likelihoods against target data and fit models in a generic manner.

Some more info on the required \option{Data} list inputs:

\describe{
  \item{Data$flux}{
Object; required. Data.frame(4 columns, filter name, central wave, flux (Jy), fluxerr (Jy))
}
  \item{Data$SFH}{
Object; required. The star formation history function to use. Must be one of \code{\link{SFHp4}}, \code{\link{SFHp5}} and \code{\link{SFHfunc}}.
}
  \item{Data$arglist}{
List; optional list of named parameters and values that are passed directly into \code{ProSpectSED}. Arguments that are not matched in \code{ProSpectSED} are passed directly onto \option{SFH}, so by this mechanism you can set function specific arguments that control the star formation history, i.e. this allows you to set values to non default even when not fitting these parameters. You almost certainly want to at least specify Data$arglist$z (the object redshift) since this almost always differ from the fiducial default of \option{z}=0.1.
}
  \item{Data$speclib}{
Object; required. Pass in the spectral library directly. Must be one of \code{\link{BC03lr}}, \code{\link{BC03hr}}, \code{\link{EMILES}}.
}
  \item{Data$Dale}{
Object; required. Pass in the Dale dust library directly (should be one of \option{Dale_NormTot}, \option{Dale_NormSFR} otherwise will run, but not work correctly).
}
  \item{Data$AGN}{
Pass in the AGN template to use. Must be one of \code{\link{AGN_UnOb}} or \code{\link{AGN_UnOb_Sparse}}.
}
  \item{Data$filtout}{
Object; required. Pass in the photometric filters directly. If \code{filters} is a vector of filter names, this can be achieved efficiently with a command like: filtout=foreach(i = filters)\%do\%{getfilt(i)} (see \code{\link{getfilt}}).
}
  \item{Data$fit}{
Character scalar; one of 'optim', 'LD' or 'check'. See Value for full information on what is output in each case. If missing it is set to 'optim'.
}
  \item{Data$like}{
Character scalar; specifies the type of likelihood. One of 'norm' (Normal), 'st' (Student-T) or 'chisq' (Chi-Square). If not one of these then all useful output is returned in a list. If missing it is set to 'st'.
}
  \item{Data$parm.names}{
Object; required. Character vector; correctly named and position matched parameters for fitting. Must match length of \option{parm} in \code{ProSpectSEDlike}, and elements should correctly describe the respective values.
}
  \item{Data$logged}{
Logical vector; match length of \option{parm} in \code{ProSpectSEDlike}, and elements should say which values should be raised to the power 10 when fitting. This is useful since scale parameters are better fit in log space, so Data$logged=TRUE for those elements.
}
  \item{Data$prior}{
Function; a prior function that takes in the \option{parm} vector and returns a log-likelihood to be added to LL to form LP. If not provided then it is ignored internally. The user can expect the vector \option{parm} to have elements specified by \option{Data$parm.names} for ease.
}
  \item{Data$intervals}{
List; has elements \option{lo} and \option{hi} which if provided must be the same length as \option{parm} in \code{ProSpectSEDlike}. These specify the lower and upper limits in the same linear/log space as \option{parm}, i.e. set \option{lo} to -1 if you do not want to logged unit to go below 0.1. 
}
  \item{Data$constraints}{
Function, a constraints function that takes in the\option{parm} vector and modifies it as required, updating with the new \option{parm} which is passed out. Using this you can control physical relationships between certain parameters if required. \option{Data$intervals} is passed after this, and always over-ride the constraints if there is a conflict. The user can expect the vector \option{parm} to have elements specified by \option{Data$parm.names} for ease.
}
  \item{Data$mon.names}{
Character vector; should either be length 1 (basically ignored), length 8 ("dustmass_birth", "dustmass_screen", "dustmass_AGN", "dustmass_total", "dustlum_birth", "dustlum_screen", "dustlum_AGN", "dustlum_total") or length 9 where the first string is 'LP'.
}
  \item{Data$verbose}{
Logical; should verbose messages be shown. If missing it is set to TRUE.
}
}
}
\value{
If \option{returnall}=FALSE then we get a vector containing the photometry in Jansky.

If \option{returnall}=TRUE, then we get a list containing:

\item{photom_out}{Vector; photometry in Jansky}
\item{FinalFlux}{Two column data.frame; Wavelength in Ang, Flux in Jansky}
\item{FinalLum}{Two column data.frame; Wavelength in Ang, Total SED luminosity in Lsol / Ang}
\item{StarsAtten}{Two column data.frame; Wavelength in Ang, Attenuated stellar SED luminosity in Lsol / Ang}
\item{StarsUnAtten}{Two column data.frame; Wavelength in Ang, Unattenuated stellar SED luminosity in Lsol / Ang}
\item{DustEmit}{Two column data.frame; Wavelength in Ang, Re-emitted stellar dust SED luminosity in Lsol / Ang}
\item{AGN}{Two column data.frame; Wavelength in Ang, Attenuated and re-emitted AGN SED luminosity in Lsol / Ang}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{SFH}}, \code{\link{magABcalc}}, \code{\link{dustmass}}, \code{\link{optim}}
}

\examples{
data('BC03lr')
data('Dale_NormTot')
data('AGN_UnOb_Sparse')
data('Dale_M2L_func')
filters=c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA',
'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1_WISE' , 'W2_WISE', 'W3_WISE', 'W4_WISE',
'P100_Herschel', 'P160_Herschel', 'S250_Herschel' , 'S350_Herschel', 'S500_Herschel')
filtout={}
for(i in filters){filtout=c(filtout,list(getfilt(i)))}

#Try playing around with the below to see what impact different AGN, dust etc have:

testSED=ProSpectSED(SFH=SFHp4, AGNlum=1e43, tau_birth=1, tau_screen=0.3, tau_AGN=2,
alpha_SF_birth=1, alpha_SF_screen=3, alpha_SF_AGN=0, speclib=BC03lr, Dale=Dale_NormTot,
AGN=AGN_UnOb_Sparse, filtout=filtout, Dale_M2L_func=Dale_M2L_func, returnall=TRUE)

plot(testSED$FinalLum, log='xy', xlim=c(1e2,1e7), ylim=c(1e2,1e7), xlab='Wave / Ang',
ylab='Lum / Lsol/A', type='l', lwd=5)
lines(testSED$StarsUnAtten, col='blue', lty=2)
lines(testSED$StarsAtten, col='green')
lines(testSED$DustEmit, col='darkgreen')
lines(testSED$AGN, col='brown')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SED }% use one of  RShowDoc("KEYWORDS")