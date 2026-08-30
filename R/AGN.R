Fritz_interp = function(lum = 1e+44, ct = 40, rm = 60, an = 30, ta = 1, al = 4, be = -0.5,
                        Fritz = NULL){
  if(is.null(Fritz)){
    Fritz=NULL
    data('Fritz', envir = environment())
  }

  ctmix = interp_quick(ct, Fritz$ct)
  almix = interp_quick(al, Fritz$al)
  bemix = interp_quick(be, Fritz$be)
  tamix = interp_quick(ta, Fritz$ta)
  rmmix = interp_quick(rm, Fritz$rm)
  anmix = interp_quick(an, Fritz$an)

  slice = Fritz$Aspec[
    c(ctmix[1:2]),
    c(almix[1:2]),
    c(bemix[1:2]),
    c(tamix[1:2]),
    c(rmmix[1:2]),
    c(anmix[1:2]),
  ]

  weights = rep(1,64)
  weights = weights * rep(ctmix[3:4], each=1, times=32)
  weights = weights * rep(almix[3:4], each=2, times=16)
  weights = weights * rep(bemix[3:4], each=4, times=8)
  weights = weights * rep(tamix[3:4], each=8, times=4)
  weights = weights * rep(rmmix[3:4], each=16, times=2)
  weights = weights * rep(anmix[3:4], each=32, times=1)
  tempmat = matrix(as.numeric(slice),64,178)
  return(data.frame(wave=Fritz$Wave, lum=lum * colSums(tempmat * weights)))
}

SKIRTOR_interp = function(lum = 1e+44, ct = 40, rm = 60, an = 30, ta = 1, p = 1, q = 1,
                          SKIRTOR = NULL){

  if(is.null(SKIRTOR)){
    data('SKIRTOR', envir = environment())
  }


  tmix = interp_quick(ta, SKIRTOR$ta)
  pmix = interp_quick(p, SKIRTOR$p)
  qmix = interp_quick(q, SKIRTOR$q)
  oamix = interp_quick(ct, SKIRTOR$ct)
  rmmix = interp_quick(rm, SKIRTOR$rm)
  imix = interp_quick(an, SKIRTOR$an)

  slice = SKIRTOR$Aspec[c(tmix[1:2]), c(pmix[1:2]), c(qmix[1:2]),
                        c(oamix[1:2]), c(rmmix[1:2]), 1, c(imix[1:2]), ]

  weights = rep(1, 64)
  weights = weights * rep(tmix[3:4], each = 1, times = 32)
  weights = weights * rep(pmix[3:4], each = 2, times = 16)
  weights = weights * rep(qmix[3:4], each = 4, times = 8)
  weights = weights * rep(oamix[3:4], each = 8, times = 4)
  weights = weights * rep(rmmix[3:4], each = 16, times = 2)
  weights = weights * rep(imix[3:4], each = 32, times = 1)

  tempmat = matrix(as.numeric(slice), 64, 132)
  agn_spectrum = (colSums(tempmat * weights))/SKIRTOR$Wave


  agn_spectrum = (agn_spectrum/(3.839e33 * 1E11))*lum
  return(data.frame(wave=SKIRTOR$Wave, lum = agn_spectrum))

}

LRDMBB_interp = function(lum = 1e+44, teff = 5000, beta = 0, taV = 1, powV = -0.7, ct = 40, rm = 60, an = 30, ta = 1, p = 1, q = 1,
                        LRDBBSKIRTOR = NULL){
  
  if(is.null(LRDBBSKIRTOR)){
    data('LRDBBSKIRTOR', envir = environment())
  }
  
  wave = LRDBBSKIRTOR$Wave ## Ang
  
  BB = greybody_norm(wave = wave, Temp = teff, beta = beta, norm = 1) ## thermal dense gas distribution 
  
  balmer_break = 1 / (1 + exp(-1 * (wave - 3646)))
  BB_balmer_break = BB * balmer_break
  BB_balmer_break_norm = BB_balmer_break / sum(c(0, diff(wave)) * BB_balmer_break)
  
  BB_atten = CF_atten(
    wave = wave, 
    flux = BB_balmer_break_norm, 
    tau = taV,
    pow = powV
  ) ## Attenuate by dust along LOS
  
  dust_emit = SKIRTOR_interp(
    lum = BB_atten$total_atten, 
    ct = ct, 
    rm = rm, 
    an = an, 
    ta = ta, 
    p = p,
    q = q,
    SKIRTOR = LRDBBSKIRTOR
  ) ## Re-emit in the IR
  
  dust_emit_norm = dust_emit$lum * BB_atten$total_atten / sum(c(0, diff(wave)) * dust_emit$lum)
  
  out = addspec(
    wave1 = wave, 
    flux1 = BB_atten$flux,
    wave2 = wave,
    flux2 = dust_emit_norm
  )
  
  return(data.frame(wave = out$wave, lum = out$flux * lum))
}

LRDLIU_interp = function(lum = 1e+44, teff = 5000, logg = -2.0, taV = 1, powV = -0.7, ct = 40, rm = 60, an = 30, ta = 1, p = 1, q = 1, 
                         LRDLIU = NULL){

  if(is.null(LRDLIU)){
    data('LRDLIU', envir = environment())
  }

  lrd_wave_raw = LRDLIU$Wave ## Ang
  
  teffmix = interp_quick(teff, LRDLIU$Teff)
  loggmix = interp_quick(logg, LRDLIU$logg)
  
  slice = LRDLIU$Aspec[c(teffmix[1:2]), c(loggmix[1:2]), ]
  
  weights = rep(1, 4)
  weights = weights * rep(teffmix[3:4], each = 1, times = 2)
  weights = weights * rep(loggmix[3:4], each = 2, times = 1)
  
  tempmat = matrix(as.numeric(slice), 4, length(lrd_wave_raw))
  lrd_spectrum_raw = (colSums(tempmat * weights))
  
  lrd_spectrum_rebin = specReBin(
    wave = lrd_wave_raw,
    flux = lrd_spectrum_raw,
    bin = 1e-3
  )

  lrd_wave = lrd_spectrum_rebin$wave
  lrd_spectrum = lrd_spectrum_rebin$flux
  
  lrd_atten = CF_atten(
    wave = lrd_wave,
    flux = lrd_spectrum,
    tau = taV,
    pow = powV
  ) ## Attenuate by dust along LOS

  dust_emit = SKIRTOR_interp(
    lum = lrd_atten$total_atten,
    ct = ct,
    rm = rm,
    an = an,
    ta = ta,
    p = p,
    q = q,
    SKIRTOR = LRDLIU$SKIRTORDUST
  ) ## Re-emit in the IR
  dust_emit_norm = dust_emit$lum / sum(c(0, diff(dust_emit$wave)) * dust_emit$lum) * lrd_atten$total_atten ## normalise to total attenuated energy

  agn_spectrum = addspec(
    wave1 = lrd_wave,
    flux1 = lrd_atten$flux,
    wave2 = dust_emit$wave,
    flux2 = dust_emit$lum,
    extrap = 0
  )
  
  # magplot(
  #   NA, 
  #   xlim = c(100, 1e8),
  #   ylim = c(1e-10, 1),
  #   log = "xy"
  # )
  # lines(
  #   lrd_wave,
  #   lrd_spectrum, 
  #   col = "cornflowerblue", 
  #   lwd = 2, 
  #   lty = 2
  # )
  # lines(
  #   lrd_wave,
  #   lrd_atten$flux,
  #   col = "blue", 
  #   lty = 2, 
  #   lwd = 2
  # )
  # lines(
  #   dust_emit$wave,
  #   dust_emit$lum, 
  #   col = "red", 
  #   lwd = 2
  # )
  # lines(
  #   agn_spectrum$wave, 
  #   agn_spectrum$flux
  # )

  return(data.frame(wave = agn_spectrum$wave, lum = agn_spectrum$flux * lum))
}

SIROCCO_interp = function(lum = 1e+44, mdot = 1, Zagn = -1.0, taV = 1, powV = -0.7, ct = 40, rm = 60, an = 30, ta = 1, p = 1, q = 1, 
                         SIROCCO = NULL){
  
  if(is.null(SIROCCO)){
    data('SIROCCO', envir = environment())
  }
  
  lrd_wave_raw = SIROCCO$Wave ## Ang
  
  mdotmix = interp_quick(mdot, SIROCCO$mdot)
  Zmix = interp_quick(Zagn, SIROCCO$Z)
  
  slice = SIROCCO$Aspec[c(mdotmix[1:2]), c(Zmix[1:2]), ]
  
  weights = rep(1, 4)
  weights = weights * rep(mdotmix[3:4], each = 1, times = 2)
  weights = weights * rep(Zmix[3:4], each = 2, times = 1)
  
  tempmat = matrix(as.numeric(slice), 4, length(lrd_wave_raw))
  lrd_spectrum_raw = (colSums(tempmat * weights))
  
  lrd_spectrum_rebin = specReBin(
    wave = lrd_wave_raw,
    flux = lrd_spectrum_raw,
    bin = 1e-3
  )
  
  lrd_wave = lrd_spectrum_rebin$wave
  lrd_spectrum = lrd_spectrum_rebin$flux
  
  lrd_atten = CF_atten(
    wave = lrd_wave,
    flux = lrd_spectrum,
    tau = taV,
    pow = powV
  ) ## Attenuate by dust along LOS
  
  dust_emit = SKIRTOR_interp(
    lum = lrd_atten$total_atten,
    ct = ct,
    rm = rm,
    an = an,
    ta = ta,
    p = p,
    q = q,
    SKIRTOR = SIROCCO$SKIRTORDUST
  ) ## Re-emit in the IR
  dust_emit_norm = dust_emit$lum / sum(c(0, diff(dust_emit$wave)) * dust_emit$lum) * lrd_atten$total_atten ## normalise to total attenuated energy
  
  agn_spectrum = addspec(
    wave1 = lrd_wave,
    flux1 = lrd_atten$flux,
    wave2 = dust_emit$wave,
    flux2 = dust_emit$lum,
    extrap = 0
  )
  
  # magplot(
  #   NA,
  #   xlim = c(100, 1e8),
  #   ylim = c(1e-10, 1),
  #   log = "xy"
  # )
  # lines(
  #   lrd_wave,
  #   lrd_spectrum,
  #   col = "cornflowerblue",
  #   lwd = 2,
  #   lty = 2
  # )
  # lines(
  #   lrd_wave,
  #   lrd_atten$flux,
  #   col = "blue",
  #   lty = 2,
  #   lwd = 2
  # )
  # lines(
  #   dust_emit$wave,
  #   dust_emit$lum,
  #   col = "red",
  #   lwd = 2
  # )
  # lines(
  #   agn_spectrum$wave,
  #   agn_spectrum$flux
  # )

  return(data.frame(wave = agn_spectrum$wave, lum = agn_spectrum$flux * lum))
}