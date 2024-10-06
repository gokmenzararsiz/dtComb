install.packages(c('devtools'))

library(devtools)
pkgload::load_all()
styler::style_pkg()
usethis::use_mit_license()
devtools::document()
devtools::build_manual()
devtools::build(args = "--compact-vignettes=gs+qpdf")
