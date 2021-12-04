exampleData1 <- read.table("data-raw/ddimer.txt", header = TRUE)

devtools::use_data(exampleData1, overwrite = TRUE)
