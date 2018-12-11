

############################################
### Start of Function 
############################################

fun_pr <- function(occu_capa, det_capa, model){
  L_occu <-(model$coefs$Value[1] +  # L pb
             (occu_capa@layers[[1]] * model$coefs$Value[2]) + 
             (occu_capa@layers[[2]] * model$coefs$Value[3]))# +
             # (occu_capa@layers[[3]] * model$coefs$Value[4]) +
             # (occu_capa@layers[[4]] * model$coefs$Value[5]) +
             # (occu_capa@layers[[5]] * model$coefs$Value[6]) +
             # (occu_capa@layers[[6]] * model$coefs$Value[7]) +
             # (occu_capa@layers[[7]] * model$coefs$Value[8]) +
             # (occu_capa@layers[[8]] * model$coefs$Value[9])) 
  
  L_det <-  (model$coefs$Value[10] + 
             (det_capa@layers[[1]] * model$coefs$Value[11]) +
             (det_capa@layers[[2]] * model$coefs$Value[12]))# +
             # (det_capa@layers[[3]] * model$coefs$Value[13]) +
             # (det_capa@layers[[4]] * model$coefs$Value[14]))
  
  L_int <- exp(L_occu * L_det)
  
    #sum(W.po %*% alpha) - sum(log(1 + exp(W.po %*% alpha))) - sum(mu*p)
  
  return(L_int)
}


############################################
### End of Function 
############################################

# use the function to predict
prediction <- fun_pr(occu_capa = s.occupancy, 
                     det_capa = s.detection, 
                     model = poANDso.fit)

### plogis function?
prediction2 <- plogis(as.matrix(prediction))

# Visualize output
# plot(prediction) # scaled values
psi.mat <- prediction2
psi.raster <- raster(psi.mat)
extent(psi.raster) <- extent(prediction)
plot(psi.raster, main="integrated Likelihood") # transformed using plogis?
plot(pb, add=T)
# hist(psi.raster, xlab = 'integrated Likelihood')



########################################
############################################
### <Cody´s solution 
############################################


linearPredictor.poANDpa.sdm = as.matrix(s.occupancy) %*% poANDso.fit$coefs[1:dim(as.matrix(s.occupancy))[2]]
predict.poANDpa.sdm= exp(linearPredictor.poANDpa.sdm)
raster.poANDpa.sdm = raster(sgrid)
cells = cellFromXY(raster.poANDpa.sdm, xy)
raster.poANDpa.sdm[cells] = predict.poANDpa.sdm

# Where, 
# X.back is a matrix of the environmental predictors
# and fit.poANDpa$par are the coefficient estimates from the model

#### what is sgrid?   xy?
