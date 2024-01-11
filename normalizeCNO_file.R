library("ggplot2")
library("CellNOptR")
library("dplyr")


# Input file as a MIDAS formatted .CSV file. Read more about MIDAS formatted file attached: MIDAS-datarail-supplement.pdf
# Read data

data = readMIDAS("data/model_fitting_conc.csv")
# data = readMIDAS("data/model_val_conc.csv")


# makeCNOlist creates a CNOlist of non-normalized data in MIDAS-formatted .csv file.

data = makeCNOlist(data, subfield=FALSE, verbose=TRUE)

# normaliseCNOlist does normalization of data between 0 and 1 using a multi-step non-linear method
# read more: https://doi.org/10.1038/msb.2009.870.75

normdata <- normaliseCNOlist(data, EC50Data=0.5, HillCoef=2, EC50Noise=0, detection=0,
                             saturation=Inf, changeTh=0, mode="time",
                             options=list(rescale_negative = T), verbose=FALSE) # default reaction parameters EC50, n.


# normalized data values
getSignals(normdata)

# save to MIDAS-formatter .csv file
writeMIDAS(normdata,"normalized_data.csv")


