library(devtools)
imports <- c("base", "stats", "graphics", "utils", "svMisc", "svDialogs", "withr", "dplyr", "doParallel", "magrittr", "foreach")
lapply(imports, use_package)


res <- resSimAWX(n_thousands = 3,params = exampleParams())
makeAWXDirectories()
