#devtools::load_all()
#library("statTarget", lib.loc = "C:/Users/bguppy/Documents/statTarget")
datpath <- "inst/extdata"
samPeno <- paste(datpath,"MTBLS79_sampleList.csv", sep="/")
samFile <- paste(datpath,"MTBLS79.csv", sep="/")
file <- paste(datpath,"data_example_two_groups.csv", sep="/")
shiftCor(samPeno, samFile, Frule = 0.8, MLmethod = "QCRLSC", QCspan = 0.75, degree = 2, imputeM = "KNN", plot = TRUE)
