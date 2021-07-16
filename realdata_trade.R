# Food trade application for "Latent space models for multiplex
# networks with shared structure"

# import multiness package
# devtools::install_github('peterwmacd/multiness')
library(multiness)

#### import data and fit model ####
set.seed(1994)

data <- log(1+agri_trade)
node_labels <- dimnames(data)[[1]]
# manual edits to labels
node_labels[39] <- 'U.A.E.'
node_labels[40] <- 'U.K.'
node_labels[41] <- 'U.S.A.'
# dimensions 
n <- dim(data)[1]
m <- dim(data)[3]

#  fit model with CV tuning
fit <- multiness_fit(data,model="gaussian",self_loops=FALSE,
                     refit=TRUE,tuning='cv',
                     tuning_opts=list(layer_wise=TRUE,
                                      penalty_const_vec=c(2,2.5,3,3.5,4,4.5,5),
                                      penalty_const_alpha=1,
                                      p_cv=.2,
                                      N_cv=5,
                                      refit_cv=TRUE,
                                      verbose_cv=TRUE),
                     optim_opts=list(max_rank=50,verbose=TRUE,return_posns=TRUE))

# recover cv constant
print(fit$lambda / (sqrt(n)*apply(data,3,denoiseR::estim_sigma,method="MAD")))
# C = 2.5 (delta = 0.5)

#### summarize latent dimensions ####

# total dimensions
print(fit$d1)
print(fit$d2)

# signs
print(fit$d1)
print(attr(fit$V_hat,'signs')) # 39
print(sum(attr(fit$V_hat,'signs')>0))
# 25 assortative, 14 disassortative
print(fit$d2)
print(sapply(1:m,function(kk){sum(attr(fit$U_hat[[kk]],'signs')>0)}))
# wine layer:
print(attr(fit$U_hat[[12]],'signs'))
# 18 total, 9 assortative
# chocolate layer:
print(attr(fit$U_hat[[3]],'signs'))
# 7 total, 5 assortative

# proportions of singular values
# common
dV <- diag(crossprod(fit$V_hat))
print(cumsum(dV)/sum(dV))
# first 4 account for 47%

# wine
dwine <- diag(crossprod(fit$U_hat[[12]]))
print(cumsum(dwine)/sum(dwine))
# first 2 account for 36%

# chocolate
dchoc <- diag(crossprod(fit$U_hat[[3]]))
print(cumsum(dchoc)/sum(dchoc))
# first 2 account for 48%

#### predefined continents for color/pch ####
asia <- c(2,7,8,9,16,17,21,24,27,29,32,36,42,91,97,101,119,122,123,127,128,132,133,141,144)
north_america <- c(6,15,41,45,52,55,58,69,72,74,89,93,95,102,114,115,120,138)
south_america <- c(5,44,50,59,86,87,88,98,99,104,110,116)
africa <- c(48,49,51,53,60,64,68,71,83,84,96,100,103,105,108,109,112,121,126,129,131,134,135,137,139,140,142,143,145)
europe <- c(3,4,10,12,13,14,19,20,22,25,28,30,33,34,35,37,38,40,43,46,47,54,56,57,61,62,65,66,67,70,73,75,76,77,78,79,80,82,90,92,94,106,107,111,117,125)
mid_east <- c(1,11,18,23,26,31,39,63,81,113,118,124,130,136)

colour_vec <- rep(NA,n)
colour_vec[asia] <- "red"
colour_vec[north_america] <- "blue"
colour_vec[south_america] <- "green"
colour_vec[africa] <- "orange"
colour_vec[europe] <- "purple"
colour_vec[mid_east] <- "cyan"

pch_vec <- rep(NA,n)
pch_vec[asia] <- 2
pch_vec[north_america] <- 13
pch_vec[south_america] <- 5
pch_vec[africa] <- 3
pch_vec[europe] <- 0
pch_vec[mid_east] <- 1

#### plot setup ####

# set an indicator for either color or bw plots
bw <- TRUE

# rules for which labels to produce:
# check the 'box' immediately to the right
label_subset <- function(X,x_space=.5,y_space=.15){
  temp <- rep(NA,nrow(X))
  for(ii in 1:nrow(X)){
    box_x <- c(X[ii,1],X[ii,1]+x_space)
    box_y <- c(X[ii,2]-.03,X[ii,2]+y_space)
    hits <- sum(box_x[1] < X[,1] & box_x[2] > X[,1] & box_y[1] < X[,2] & box_y[2] > X[,2])
    if(hits > 0){
      temp[ii] <- FALSE
    }
    else{
      temp[ii] <- TRUE
    }
  }
  # never label unspecified
  temp[85] <- FALSE
  return(temp)
}

# set margins
par(oma=c(1,1,0,0))
if(bw){
  par(mar=c(3,3,2,2))
}
if(!bw){
  par(mar=c(2,2,2,2))
}
doublepane <- TRUE
if(doublepane){
  par(mfrow=c(1,2))
}
if(!doublepane){
  par(mfrow=c(1,1))
}

#### plot first common dimensions ####

# plot 1st and 2nd common dimensions
V12 <- fit$V_hat[,1:2]
plot(V12,xlim=c(0,5),ylim=c(-1.8,1.7),
     main=ifelse(bw,NA,"Common latent dimensions 1+2"),
     xlab="",ylab="",
     col=ifelse(rep(bw,n),1,colour_vec),
     pch=ifelse(rep(bw,n),pch_vec,1))
do_label <- label_subset(V12)
do_label[c(4,10,35,46,47,54,57,62,66)] <- FALSE
text(V12[do_label,1],V12[do_label,2],node_labels[do_label],pos=4,cex=.6)

# plot 3rd and 4th common dimensions
V34 <- fit$V_hat[,3:4]
plot(V34,xlim=c(-1.5,2.5),
     main=ifelse(bw,NA,"Common latent dimensions 3+4"),
     xlab=,ylab="",
     col=ifelse(rep(bw,n),1,colour_vec),
     pch=ifelse(rep(bw,n),pch_vec,1))
do_label <- label_subset(V34,x_space=.45,y_space=.075)
do_label[c(20,25,27,43,53,79,80,81,118,123,129,141)] <- FALSE
text(V34[do_label,1],V34[do_label,2],node_labels[do_label],pos=4,cex=.6)

if(bw){
mtext('assortative',side=1,line=-1,outer=T,at=.27)
mtext('assortative',side=1,line=-1,outer=T,at=.77)
mtext('assortative',side=2,line=-1,outer=T,at=.5)
mtext('assortative',side=2,line=2,outer=F,at=.05)
}

#### plot wine layer individual dimensions ####

U_wine <- fit$U_hat[[12]][,1:2]
U_wine[,1] <- -U_wine[,1]
plot(U_wine,xlim=c(-2.9,4.2),ylim=c(-1.5,3.8),
     main=ifelse(bw,NA,"Individual latent dimensions 1+2, wine"),
     xlab="",ylab="",
     col=ifelse(rep(bw,n),1,colour_vec),
     pch=ifelse(rep(bw,n),pch_vec,1))
if(!bw){axis(1,col='red',col.axis='red')}
do_label4 <- do_label2 <- rep(FALSE,nrow(U_wine))
do_label4[c(88,44,75,100,97,104,17,11,31,41)] <- TRUE
do_label2[c(40,13,8,36,2)] <- TRUE
text(U_wine[do_label4,1],U_wine[do_label4,2],node_labels[do_label4],pos=4,cex=.6)
text(U_wine[do_label2,1],U_wine[do_label2,2],node_labels[do_label2],pos=2,cex=.6)

# compare to ase embedding
U_wine_ase <- ase(data[,,12],2)
U_wine_ase[,1] <- -U_wine_ase[,1]
plot(U_wine_ase[,2],U_wine_ase[,1],ylim=c(-.5,5),xlim=c(-1.5,4.4),
     main=ifelse(NA,"ASE latent dimensions 2+1, wine"),
     xlab="",ylab="",
     col=ifelse(rep(bw,n),1,colour_vec),
     pch=ifelse(rep(bw,n),pch_vec,1))
if(!bw){axis(1,col='red',col.axis='red')}
do_label <- rep(FALSE,nrow(U_wine_ase))
do_label[c(13,20,88,2,41,40,25,4,97,44,62,39,104,47,33)] <- TRUE
text(U_wine_ase[do_label,2],U_wine_ase[do_label,1],node_labels[do_label],pos=4,cex=.6)

if(bw){
  mtext('disassortative',side=1,line=-1,outer=T,at=.27)
  mtext('disassortative',side=1,line=-1,outer=T,at=.77)
  mtext('assortative',side=2,line=-1,outer=T,at=.5)
  mtext('assortative',side=2,line=2,outer=F,at=2.1)
}

#### plot chocolate layer individual dimensions ####

U_choc <- fit$U_hat[[3]][,1:2]
U_choc <- -U_choc
plot(U_choc,ylim=c(-1.5,2.5),xlim=c(-1.3,3.5),
     main=ifelse(bw,NA,"Individual latent dimensions 1+2, chocolate"),
     xlab="",ylab="",
     col=ifelse(rep(bw,n),1,colour_vec),
     pch=ifelse(rep(bw,n),pch_vec,1))
if(!bw){axis(2,col='red',col.axis='red')}
do_label <- label_subset(U_choc,x_space=.5,y_space=.15)
do_label[c(20,33,58,27,47,54)] <- FALSE
text(U_choc[do_label,1],U_choc[do_label,2],node_labels[do_label],pos=4,cex=.6)
text(U_choc[37,1],U_choc[37,2],node_labels[37],pos=2,cex=.6)

#compare to ase embedding
U_choc_ase <- ase(data[,,3],2)
#U_choc_ase[,1] <- -U_choc_ase[,1]
plot(U_choc_ase,xlim=c(-.4,5),ylim=c(-2,2.4),
     main=ifelse(bw,NA,"ASE latent dimensions 1+2, chocolate"),
     xlab="",ylab="",
     col=ifelse(rep(bw,n),1,colour_vec),
     pch=ifelse(rep(bw,n),pch_vec,1))
do_label <- label_subset(U_choc_ase,x_space=.5,y_space=.15)
do_label[4] <- TRUE
do_label[c(17,21,26,118,9,29,67,65,76,82,57,106,92,47,38,20)] <- FALSE
text(U_choc_ase[do_label,1],U_choc_ase[do_label,2],node_labels[do_label],pos=4,cex=.6)

if(bw){
  mtext('assortative',side=1,line=-1,outer=T,at=.27)
  mtext('assortative',side=1,line=-1,outer=T,at=.77)
  mtext('disassortative',side=2,line=-1,outer=T,at=.5)
  mtext('assortative',side=2,line=2,outer=F,at=0)
}

# save each plot as 9x6 potrait pdf
# common[_bw].pdf
# wine[_bw].pdf
# choc[_bw].pdf

#### edge imputation task ####

set.seed(1996)

p_holdout <- .05*1:8
n_holdout <- length(p_holdout)
nreps <- 5

# two approaches: multiness and svd
errors_mness <- array(NA,c(m,n_holdout,nreps))
errors_svd <- array(NA,c(m,n_holdout,nreps))

for(ll in 1:m){
  for(ii in 1:n_holdout){
    for(jj in 1:nreps){
      # hold out a stratified sample of zero and non-zero edges
      # use MP_nz (non-zero) to evaluate holdout error,
      # fit with MP (both non-zero and zero edges held out)
      MP <- array(TRUE,c(n,n,m))
      MP_nz <- multiness:::missing_mat_sym(n,p_holdout[ii],subset=(data[,,ll]>0))
      MP_z <- multiness:::missing_mat_sym(n,p_holdout[ii],subset=(data[,,ll]==0))
      MP[,,ll] <- as.logical(MP_nz*MP_z)

      # fit the multiness_model
      fit <- multiness_fit(data,model='gaussian',
                           self_loops=FALSE,refit=TRUE,
                           tuning='adaptive',
                           tuning_opts=list(layer_wise=TRUE,
                                            penalty_const=2.5),
                           optim_opts=list(max_rank=100,
                                           missing_pattern=MP,
                                           verbose=FALSE))
      fit_mness <- list(F_hat=fit$F_hat,G_hat=list(fit$G_hat[[ll]]))

      # fit a marginal low-rank matrix completion without the other layers
      # (hard singular value thresholding after inflating entries,
      # using a slightly larger lambda and only keeping positive eigenvalues)
      fit_svd <- list()
      inflation <- 1/(1-p_holdout[ii])
      temp <- multiness:::sv_thresh_f(inflation*data[,,ll]*MP[,,ll],1.2*fit$lambda[ll],
                                      max_rank=100,soft=FALSE,pos=TRUE)
      fit_svd$F_hat <- matrix(0,n,n)
      fit_svd$G_hat <- list()
      fit_svd$G_hat[[1]] <- multiness:::einfo_to_mat(temp)

      error_mness <- multiness:::holdout_error(array(data[,,ll],c(n,n,1)),fit_mness,
                                    hollow=TRUE,misspattern = array(MP_nz,c(n,n,1)))

      error_svd <- multiness:::holdout_error(array(data[,,ll],c(n,n,1)),fit_svd,
                                      hollow=TRUE,misspattern = array(MP_nz,c(n,n,1)))

      # RMSE normalized and scaled
      print(c(ll,ii,jj))
      errors_mness[ll,ii,jj] <- sqrt(error_mness/sum(!MP_nz))
      errors_svd[ll,ii,jj] <- sqrt(error_svd/sum(!MP_nz))
    }
  }
}

# collapse repeated dimension
errors_mness_mean <- apply(errors_mness,c(1,2),mean)
rownames(errors_mness_mean) <- dimnames(data)[[3]]
colnames(errors_mness_mean) <- p_holdout

errors_svd_mean <- apply(errors_svd,c(1,2),mean)
rownames(errors_svd_mean) <- dimnames(data)[[3]]
colnames(errors_svd_mean) <- p_holdout

#### PLOTTING ####

# plot results

par(mfrow=c(2,3),oma=c(2,2,0,0),mar=c(2,2.5,2,2.5))

p_holdout <- .05*1:8
n_holdout <- length(p_holdout)
nlay_plot <- 6
layer_select <- c(12,3,6,10,11,2)
layer_names <- c('Wine','Chocolate','Distilled Alcohol',
                 'Tea','Sugar','Prepared Foods')
# prepared foods including things like soups and condiments

# black and white indicator for plots
bw <- FALSE

for(ii in 1:nlay_plot){
  plot(p_holdout,errors_svd_mean[layer_select[ii],],
       main=layer_names[ii],xlab="",ylab="",type='b',
       col=ifelse(bw,'black','red'),pch=1,lty=2,
       ylim=c(0,3.5))
  lines(p_holdout,errors_mness_mean[layer_select[ii],],type="b",
        col=ifelse(bw,'black','blue'),
        pch=17,lty=1)
  #lines(p_holdout,errors_partial_mean[layer_select[ii],],type="b",col="orange")
  #lines(p_holdout,errors_extra_mean[layer_select[ii],],type="b",col="green")
  if(ii==4 & !bw){
    legend(0.06,1.3,legend=c('SVD','MultiNeSS+'),col=c('red','blue'),
           lty=c(2,1),pch=c(1,17),bty='n',ncol=1,cex=.9,xjust=0)
  }
}
mtext(ifelse(bw,'Root mean squared error','RMSE'),2,outer=T,line=.2)
mtext('Proportion missing',1,outer=T,line=.2)

# export as 5x7 landscape pdf
# layer_holdout[_bw].pdf

