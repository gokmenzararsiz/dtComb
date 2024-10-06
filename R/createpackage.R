install.packages(c('devtools','rgl','glmnet','gam','fastAdaboost','doParallel','adabag','plotrix','Formula','plotmo','TeachingDemos','earth','mda','rprojroot','diffobj','rematch2','brio','callr','desc','pkgload','praise','processx','ps','waldo','testthat','minqa','nloptr','lme4','abind','coda','ada','mvtnorm','stabs','nnls','quadprog','import','libcoin','inum','partykit','mboost','bitops','caTools','RWeka','party','coin','strucchange','multcomp','sandwich','modeltools','matrixStats','TH.data','kerndwd','xgboost','RSpectra','rARPACK','HDclassif','kknn','klaR','questionr','labelled','styler','haven','R.cache','readr','R.utils','vroom','R.oo','bit64','combinat','rstudioapi','miniUI','forcats','R.methodsS3','clipr','bit','optimx','monmlp','RSNNS','ncvreg','msaenet','naivebayes','pamr','randomForest','pls','dotCall64','gridExtra','backports','statnet.common','maps','SparseM','MatrixModels','permute','carData','network','spam','viridis','broom','vegan','quantreg','sna','fields','pbkrtest','bipartite','car','plsRglm','stepPlr','ordinalNet','RRF','LiblineaR','DEoptimR','pcaPP','robustbase','rrcov','truncnorm','mclust','Rsolnp','robustDA','rotationForest','kohonen','entropy','corpcor','fdrtool','sda','sdwd','lars','elasticnet','sparseLDA','spls','deepnet','gbm','evtree','wsrf','readxl','stringr','ggplot2','dplyr','DT','shinyBS','plyr','pROC','epiR','OptimalCutpoints','ROCR','pls','mda','shinyalert','combinat','rapport','arm','Cubist','C50','kernlab','misc3d','plot3D','prim','supervisedPRIM','munsell'))

library(devtools)
pkgload::load_all()
styler::style_pkg()
usethis::use_mit_license()
devtools::document()
devtools::build_manual()
devtools::build(args = "--compact-vignettes=gs+qpdf")
