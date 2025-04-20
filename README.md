# LatentOutcomes
R package for Estimating the Causal Effects with Latent Outcomes

Author: Jiawei Fu; Donald P. Green

Contributor: Jiawei Fu; Donald P. Green

The package is used to estimate causak effects with latent outcomes.

To install and use the latest version of the package, try the following codes:
```r
install.packages("devtools")
devtools::install_github("Jiawei-Fu/LatentOutcomes")
library(LatentOutcomes)
```
Some of users will receive warnings about time zone, it will not affect the installation.

There are two functions and one data set in the package now.

Example:
```r
data(test_data)  # input data
estlatent(test_dat$Z,test_dat[,1:3],X=NULL,eta = 1,method="sem",IV_Y=T,tau=T)   # use SEM estimation
estlatent(test_dat$Z,test_dat[,1:3],X=NULL,eta = 1,method="gmm_equal",IV_Y=T,tau=T) use equally weighted index estimation
estlatent(test_dat$Z,test_dat[,1:3],X=NULL,eta = 1,method="gmm_opt",IV_Y=T,tau=T) # use optimally weighted index estimation
```

If you have any comments, suggestions, or findings, appreciate sending me the email: jiawei.fu@yale.edu

