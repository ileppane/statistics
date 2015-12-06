# 28.10.2014
#
# Epps-Singleton (1986) nonparametric two-sample omnibus test for differences in empirical characteristic functions.
#
# x,y are the data vectors
# tau and J are parameters that have the default 'omnibus' values calculated by Epps and Singleton
#  you can also set them yourself if you don't like them
# scalecor is a logical option for doing a scale invariance correction (default=TRUE).
#
# Translated and modified from Matlab code of Christian Rojas (University of Massachusetts).
#
# This is a 'two-sided' test. If you have an idea how to make this 'one-sided' and allow
# for directional hypothesis testing, please let me know! 
#
# You can read this into your R analysis file by command: source("http://users.aalto.fi/ileppane/EStest.R")
#
# (ilkka.j.leppanen@aalto.fi)
#
es.test <- function(x,y,tau=0.4,J=2,scalecor=TRUE) {
 Tx <- length(x);
 Ty <- length(y);
 # Create the scale invariance correction for the t's
 sorted <- sort(c(x,y));
 if (scalecor) {
  sigmahat <- 0.5*(  (sorted[floor((Tx+Ty)*0.25)] + sorted[floor((Tx+Ty)*0.25)+1])/2 + (sorted[floor((Tx+Ty)*0.75)] + sorted[floor((Tx+Ty)*0.75)-1])/2  );
  that <- (1:J)*tau/sigmahat;
 } else that <- (1:J)*tau;
 # Create the 2J x Ti matrices (J=2 and Ti is the size of sample i), equation (3) in paper
 J <- length((1:J)*tau);
 gxm <- array(0, dim=c(2*J,Tx))
 gym <- array(0, dim=c(2*J,Ty))
 for (i in 1:Tx) {
  gxm[seq(1,2*J,by=2),i] <- cos(that*x[i]);
  gxm[seq(2,2*J,by=2),i] <- sin(that*x[i]);
 }
 for (i in 1:Ty) {
  gym[seq(1,2*J,by=2),i] <- cos(that*y[i]);
  gym[seq(2,2*J,by=2),i] <- sin(that*y[i]);
 }
 # Create the g_i vectors
 gx <- (1/Tx)*rowSums(gxm);
 gy <- (1/Ty)*rowSums(gym);
 # Covariances (these are the biased versions and they are used by the test)
 Sxhat <- ((Tx-1)/Tx)*cov(t(gxm));
 Syhat <- ((Ty-1)/Ty)*cov(t(gym));
 # Compute Omega hat
 omegahat <- ((Tx+Ty)/Tx)*Sxhat + ((Tx+Ty)/Ty)*Syhat;
 # Compute the ES statistic
 if (Tx<25 & Ty<25) {   # use the small sample correction
  ssc <- TRUE;
  cTxTy <- 1/(1+(Tx+Ty)^(-0.45)+10.1*(Tx^(-1.7)+Ty^(-1.7)));
  W2 <- cTxTy*(Tx+Ty)%*%t((gx-gy))%*%solve(omegahat)%*%(gx-gy);
  } else {
   ssc <- FALSE;
   W2 <- (Tx+Ty)*t((gx-gy))%*%solve(omegahat)%*%(gx-gy);
  }
 # Compute p-value using the statistic and the degrees of freedom
 # pvalue <- 1-pchisq(W2,qr(omegahat)$rank);
 pvalue <- 1-pchisq(W2,qr(omegahat)$rank);
 meth <- ifelse(ssc,"Epps-Singleton two-sample omnibus test using small sample correction","Epps-Singleton two-sample omnibus test");
 names(W2) <- "Test statistic";
 rval <- list(statistic=W2,
              p.value=pvalue,
              alternative=paste("ECDF of", deparse(substitute(x)), "does not equal ECDF of", deparse(substitute(y))),
              method=meth,
              data.name=paste(deparse(substitute(x)), "and", deparse(substitute(y))));
 class(rval) <- "htest";
 
 return(rval);
}



## 
# (This part does not belong to the function)
# Data from Goerg and Kaiser (2009, The Stata Journal)
# germany <- c( 0, 1, 2, 2,  2,  2, 3, 3, 3, 3,  4, 5, 5, 5, 5, 6, 10, 10, 10, 10);
# china <-   c(10, 4, 0, 5, 10, 10, 0, 5, 6, 7, 10, 3, 1, 7, 0, 8,  1,  5,  8, 10);
# es.test(germany,china)
##
