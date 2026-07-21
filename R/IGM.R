## Inoue+14 IGM attenuation 
## https://arxiv.org/abs/1402.0677
## Lyman continuum and Lyman series absorption given distribution function of z and NHI
## Consider both Lyman alpha forest and damped lyman systems

## Load in Lyman series coefficients for Eqs, 21-22
## Data obtained from https://github.com/gbrammer/eazy-py/tree/master/eazy/data

tau_IGM_LSLAF <- function(wave, z, Inoue14_LAFcoef = NULL) {
  if (is.null(Inoue14_LAFcoef)) {
    data("Inoue14_LAFcoef", envir = environment())
  }
  
  LJ  <- Inoue14_LAFcoef$lambda
  AJ1 <- Inoue14_LAFcoef$AJ1
  AJ2 <- Inoue14_LAFcoef$AJ2
  AJ3 <- Inoue14_LAFcoef$AJ3
  n   <- length(LJ)
  
  tau_LAF <- numeric(length(wave))
  zfac <- 1 + z
  
  for (i in seq_len(n)) {
    Lj <- LJ[i]
    upper <- Lj * zfac
    
    valid <- wave > Lj & wave < upper
    if (!any(valid)) next          # skip lines with no contribution at all
    
    w <- wave[valid]
    ratio <- w / Lj
    
    ttau <- numeric(length(w))
    r1 <- ratio < 2.2
    r2 <- !r1 & ratio < 5.7
    r3 <- !r1 & !r2
    
    if (any(r1)) ttau[r1] <- AJ1[i] * ratio[r1]^1.2
    if (any(r2)) ttau[r2] <- AJ2[i] * ratio[r2]^3.7
    if (any(r3)) ttau[r3] <- AJ3[i] * ratio[r3]^5.5
    
    tau_LAF[valid] <- tau_LAF[valid] + ttau
  }
  
  tau_LAF
}

tau_IGM_LSDLA <- function(wave, z, Inoue14_DLAcoef = NULL) {
  if (is.null(Inoue14_DLAcoef)) {
    data("Inoue14_DLAcoef", envir = environment())
  }
  
  LJ  <- Inoue14_DLAcoef$lambda
  AJ1 <- Inoue14_DLAcoef$AJ1
  AJ2 <- Inoue14_DLAcoef$AJ2
  n   <- length(LJ)
  
  tau_DLA <- numeric(length(wave))
  zfac <- 1 + z
  
  for (i in seq_len(n)) {
    Lj <- LJ[i]
    upper <- Lj * zfac
    
    valid <- wave > Lj & wave < upper
    if (!any(valid)) next
    
    w <- wave[valid]
    ratio <- w / Lj
    
    ttau <- numeric(length(w))
    r1 <- ratio < 3.0
    r2 <- !r1
    
    if (any(r1)) ttau[r1] <- AJ1[i] * ratio[r1]^2.0
    if (any(r2)) ttau[r2] <- AJ2[i] * ratio[r2]^3.0
    
    tau_DLA[valid] <- tau_DLA[valid] + ttau
  }
  
  tau_DLA
}

tau_IGM_LCLAF <- function(wave, z) {
  
  ## Lyman continuum from Lyman alpha forest
  
  LL <- 911.75  ## hardcode Lyman limit
  zshift_LL <- LL * (1 + z)
  
  zLAF1 <- 1.2
  zLAF2 <- 4.7
  
  tau_LAF <- numeric(length(wave))
  
  if (z < zLAF1) {
    ## z < 1.2
    idx <- wave < zshift_LL & wave > LL
    if (any(idx)) {
      w <- wave[idx] / LL
      tau_LAF[idx] <- 0.325 * (w^1.2 - (1 + z)^(-0.9) * w^2.1)
    }
    
  } else if (z < zLAF2) {
    ## 1.2 <= z < 4.7
    idx1 <- wave < 2.2 * LL & wave > LL
    if (any(idx1)) {
      w <- wave[idx1] / LL
      tau_LAF[idx1] <- 2.55e-2 * (1 + z)^1.6 * w^2.1 + 0.325 * w^1.2 - 0.250 * w^2.1
    }
    idx2 <- wave >= 2.2 * LL & wave < zshift_LL
    if (any(idx2)) {
      w <- wave[idx2] / LL
      tau_LAF[idx2] <- 2.55e-2 * ((1 + z)^1.6 * w^2.1 - w^3.7)
    }
    
  } else {
    ## z >= 4.7
    idx1 <- wave < 2.2 * LL & wave > LL
    if (any(idx1)) {
      w <- wave[idx1] / LL
      tau_LAF[idx1] <- 5.22e-4 * (1 + z)^3.4 * w^2.1 + 0.325 * w^1.2 - 3.14e-2 * w^2.1
    }
    idx2 <- wave >= 2.2 * LL & wave < 5.7 * LL
    if (any(idx2)) {
      w <- wave[idx2] / LL
      tau_LAF[idx2] <- 5.22e-4 * (1 + z)^3.4 * w^2.1 + 0.218 * w^2.1 - 2.55e-2 * w^3.7
    }
    idx3 <- wave >= 5.7 * LL & wave < zshift_LL
    if (any(idx3)) {
      w <- wave[idx3] / LL
      tau_LAF[idx3] <- 5.22e-4 * ((1 + z)^3.4 * w^2.1 - w^5.5)
    }
  }
  
  tau_LAF
}

tau_IGM_LCDLA <- function(wave, z) {
  LL <- 911.75
  zshift_LL <- LL * (1 + z)
  zDLA <- 2.0
  
  tau_DLA <- numeric(length(wave))
  
  if (z < zDLA) {
    idx <- wave < zshift_LL & wave > LL
    if (any(idx)) {
      w <- wave[idx] / LL
      tau_DLA[idx] <- 0.211 * (1 + z)^2.0 - 7.66e-2 * (1 + z)^2.3 * w^-0.3 - 0.135 * w^2.0
    }
  } else {
    idx1 <- wave < 3.0 * LL & wave > LL
    if (any(idx1)) {
      w <- wave[idx1] / LL
      tau_DLA[idx1] <- 0.634 + 4.7e-2 * (1 + z)^3 - 1.78e-2 * (1 + z)^3.3 * w^-0.3 -
        0.135 * w^2.0 - 0.291 * w^-0.3
    }
    idx2 <- wave >= 3.0 * LL & wave < zshift_LL & wave > LL
    if (any(idx2)) {
      w <- wave[idx2] / LL
      tau_DLA[idx2] <- 4.7e-2 * (1 + z)^3 - 1.78e-2 * (1 + z)^3.3 * w^-0.3 - 2.92e-2 * w^3.0
    }
  }
  tau_DLA
}

tau_IGM_Tot = function(wave, z, Inoue14_LAFcoef = NULL, Inoue14_DLAcoef = NULL){
  tau = tau_IGM_LSLAF(wave, z, Inoue14_LAFcoef) +
    tau_IGM_LSDLA(wave, z, Inoue14_DLAcoef) + 
    tau_IGM_LCLAF(wave, z) + 
    tau_IGM_LCDLA(wave, z)
  return(tau)
}

Inoue14_IGM = function(wave, z, Inoue14_LAFcoef = NULL, Inoue14_DLAcoef = NULL){
  ## Below wave_observed<911.75 is unconstrained, manually set to 0
  igm = ifelse(
    wave<911.75,
    0,
    exp(-1 * tau_IGM_Tot(wave, z, Inoue14_LAFcoef, Inoue14_DLAcoef))
  )
  return(igm)
}