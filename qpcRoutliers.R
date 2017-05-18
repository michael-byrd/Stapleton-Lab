#looking for qpcR outliers

#ultimately want to use KOD, but first we need to construct a modlist
#arguments we need for a modlist
# dataframe containing qpcR data
ucnfmod<- modlist(ucnf, cyc = 2, fluo = mufcol, model = l5, remove = "KOD")
plot(ucnfmod, which = "single")

ucnfKOD<- KOD(ucnfmod, method = "uni1", remove = TRUE, plot = TRUE)
plot(ucnfKOD)
is.outlier(ucnfKOD)
#this does not mark any of them as outliers, i assume this is because the majority are low, and the sample
#does not have a high proportion of successful tests
#let's try this technique on the other data set with 3/4 successful trials out of 5 instead
ucnhmod<- modlist(ucnh, cyc = 2, fluo = mufcol, model = l5)
plot(ucnhmod, which = 'single')
ucnhKOD<- KOD(ucnhmod, method = "uni1", plot = TRUE)
is.outlier(ucnhKOD)
#This actually gave us D8 as an outlier