# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

.wave_rebin_cpp <- function(wave, wave_bin, wave_bin_lo, wave_bin_hi, logbin = FALSE) {
    .Call(`_ProSpect_wave_rebin_cpp`, wave, wave_bin, wave_bin_lo, wave_bin_hi, logbin)
}

.spec_rebin_cpp <- function(wave_in, flux_in, wave_out, logbin_in = FALSE, logbin_out = FALSE) {
    .Call(`_ProSpect_spec_rebin_cpp`, wave_in, flux_in, wave_out, logbin_in, logbin_out)
}

