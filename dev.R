library(devtools)
imports <- c("base", "stats", "graphics", "utils", "svMisc", "svDialogs", "withr", "dplyr", "doParallel", "magrittr", "foreach")
lapply(imports, use_package)



####### Example usage
rm(list = ls())
library(impatience2)
res <- resSimAWX(n_thousands = 3,params = exampleParams())
makeAWXDirectories()
negLogLik(1,1,1,res)
mleFull(res,exampleParams())
mleLiron(res)
mleOneKnown(res,exampleParams(),"gamma")
mleOneKnown(res,exampleParams(),"lambda_0")
mleOneKnown(res,exampleParams(),"theta")
f <- oneKnownLik(which_known = "gamma",exampleParams(),AWX = res)
f(c(3,20))
params <- exampleParams()
mleTwoKnown(res,params,"gamma")
mleTwoKnown(res,params,"theta")
mleTwoKnown(res,params,"lambda_0")
