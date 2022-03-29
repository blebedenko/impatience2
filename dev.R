library(devtools)
imports <- c("base", "stats", "graphics", "utils", "svMisc", "svDialogs", "withr", "dplyr", "doParallel", "magrittr", "foreach")
lapply(imports, devtools::use_package)


par1 <- exampleParams()
gl_L0T <- oneKnownGrad(which_known = "gamma",params = par1,AWX = AWX)
gl_L0T(1,1)

gradOneKnown(par12 = c(lambda_0=1,theta=1),params = par1, AWX=AWX)


params <- exampleParams()
names(GL0T(params))
AWX <- exampleDataAWX()
point <- GL0T(params)[c("lambda_0","theta")]
which_known <- "gamma"
grad(point = point,AWX = AWX,params = params)
point <- GL0T(params)[c("lambda_0","gamma")]
grad(point = point,AWX = AWX,params = params)
point <- GL0T(params)[c("gamma","lambda_0")]

grad(point = point,AWX = AWX,params = params)
point <- GL0T(params) + runif(3)
grad(point = point,AWX = AWX,params = params)

mle1(AWX = AWX,params = params,which_known = "gamma")
