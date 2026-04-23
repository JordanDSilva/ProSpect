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

LRD_interp = function(lum = 1e+44, temp = 5000, taV = 1, powV = -0.7, ct = 40, rm = 60, an = 30, ta = 1, p = 1, q = 1,
                      LRD = NULL){
  
  wave = LRD$Wave ## Ang
  
  BB = blackbody_norm(wave = wave, Temp = temp, norm = 1) ## thermal dense gas distribution 
  
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
    SKIRTOR = LRD
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

LIU_interp = function(lum = 1e+44, temp = 5000, logg = -2.0, taV = 1, powV = -0.7, ct = 40, rm = 60, an = 30, ta = 1, p = 1, q = 1, LRDLIU = NULL){

  if(is.null(LRDLIU)){
    data('LRDLIU', envir = environment())
  }
  
  lrd_wave = LRDLIU$Wave ## Ang
  
  teffmix = interp_quick(temp, LRDLIU$Teff)
  loggmix = interp_quick(logg, LRDLIU$logg)
  
  slice = LRDLIU$Aspec[c(teffmix[1:2]), c(loggmix[1:2]), ]
  
  weights = rep(1, 4)
  weights = weights * rep(teffmix[3:4], each = 1, times = 2)
  weights = weights * rep(loggmix[3:4], each = 2, times = 1)
  
  tempmat = matrix(as.numeric(slice), 4, length(lrd_wave))
  lrd_spectrum = (colSums(tempmat * weights))
  
  # lrd_spectrum_slope = mean(tail(diff(log10(lrd_spectrum))/diff(log10(lrd_wave)), 10))
  lrd_spectrum_slope = -4
  
  ## put Rayleigh Jeans tail to extrapolate
  RJ_wave = 10^seq(3, 8, 0.0001)
  RJ_tail = 10^(log10(tail(lrd_spectrum, 1)) + lrd_spectrum_slope*(log10(RJ_wave) - log10(tail(lrd_wave,1))))
  # RJ_tail = blackbody_norm(RJ_wave, Temp = temp)
  
  agn_spectrum = c(
    lrd_spectrum,
    RJ_tail[RJ_wave > tail(lrd_wave,1)]
  )
  waveout = c(
    lrd_wave,
    RJ_wave[RJ_wave > tail(lrd_wave,1)]
  )
  # magplot(
  #   waveout, agn_spectrum,
  #   log = "xy",
  #   type = "l",
  #   # xlim = c(1e4, 1e5),
  #   # ylim = c(1e-7, 1e-5),
  #   lwd = 2
  # )
  # lines(
  #   lrd_wave, lrd_spectrum, col = "red"
  # )
  # lines(
  #   RJ_wave, RJ_tail, col = "purple"
  # )

  agn_atten = CF_atten(
    wave = waveout,
    flux = agn_spectrum,
    tau = taV,
    pow = powV
  ) ## Attenuate by dust along LOS

  dust_emit = SKIRTOR_interp(
    lum = agn_atten$total_atten,
    ct = ct,
    rm = rm,
    an = an,
    ta = ta,
    p = p,
    q = q,
    SKIRTOR = LRDLIU$SKIRTORDUST
  ) ## Re-emit in the IR
  dust_emit_norm = dust_emit$lum * agn_atten$total_atten / sum(c(0, diff(dust_emit$wave)) * dust_emit$lum)

  out = addspec(
    wave1 = waveout,
    flux1 = agn_atten$flux,
    wave2 = dust_emit$wave,
    flux2 = dust_emit_norm,
    extrap = 0
  )
  
  # magplot(
  #   waveout, agn_spectrum,
  #   log = "xy",
  #   type = "l",
  #   xlim = c(1e2, 1e6),
  #   ylim = c(1e-15, 1e-3),
  #   lwd = 2
  # )
  # lines(
  #   waveout, agn_spectrum
  # )
  # lines(
  #   out$wave, out$flux, col = "blue"
  # )

  return(data.frame(wave = out$wave, lum = out$flux * lum))
  # return(data.frame(wave = waveout, lum = agn_spectrum * lum))
}
