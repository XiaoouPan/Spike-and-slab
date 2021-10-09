### Multiple testing with spike and slab priors on censored data

post_sas = function(Y, n, p, mu_h0, n.adapt = 5000, n.burn = 5000, n.iter = 10000) {
  p1 = 0.25 + 0.5 * as.numeric(rowMeans(response) > p0 + 0.15)
  dat = list(Y = Y,
             n = n,
             p = p,
             mu_h0 = mu_h0)
  thismodel = try(jags.model(file = "bugs/sas.txt", 
                             data = dat, 
                             inits = list(Z = Z,
                                          diff1 = 1,
                                          diff2 = 1,
                                          ss1 = ss10,
                                          ss2 = ss20,
                                          tau1 = 0.001,
                                          tau2 = 0.001),
                             n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
  try(update(thismodel, n.burn, progress.bar = "none"), silent = TRUE)
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c("mu1", "mu2", "rho"),
                              n.iter = n.iter, progress.bar = "none"), silent = TRUE)
  mu1_rec = matrix(res.bugs$mu1, N, n.iter)
  mu2_rec = matrix(res.bugs$mu2, N, n.iter)
  rho = matrix(res.bugs$rho, N, n.iter)
  return (list("mu1_rec" = mu1_rec, "mu2_rec" = mu2_rec, "rho" = rho))
}
