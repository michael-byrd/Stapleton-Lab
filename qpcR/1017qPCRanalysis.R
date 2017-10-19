#####10/17/17 qPCR Analysis#####
library(qpcR)
f1017<- read.csv(file = "C:/Users/Thomas/Documents/PCR data/1017qPCR/fam1017.csv", header = TRUE)
f1017[,1] = NULL
m <- modlist(x = f1017, cyc = 1, model = l4, remove = "none")
groups = vector('numeric')
for(x in (1:30)){
  groups = append(groups, rep(x,4))
}
groups
freps <- replist(m, group = groups, model = l4, remove = "none")
plot(freps)
fef = lapply(freps,function(x){
  efficiency(object = x, plot = FALSE)
})
fusing <-list('A5','C5','D1','E5','F1','G5','G9','H1',
                 'I5','I9','J1','K5','K9','L1','P5','P9')

h1017<- read.csv(file = "C:/Users/Thomas/Documents/PCR data/1017qPCR/hex1017.csv", header = TRUE)
h1017[,1] = NULL
hm <- modlist(x = h1017, cyc = 1, model = l4)
groups = vector('numeric')
for(x in (1:30)){
  groups = append(groups, rep(x,4))
}
groups
hreps <- replist(hm, group = groups, model = l4)
plot(hreps)
hef = lapply(hreps,function(x){
  efficiency(object = x, plot = FALSE)
})

husing = list('A1','A9','A13','A17','B5','C1','C9','C13','C17','D5','E1','E9','E13','E17','F5',
              'G1','G9','G13','G17','H5','I1','I9','I13','I17','J5','K1','K9','L5','P1','P9')

shared = setdiff(husing,setdiff(husing,fusing))
#which runs in hex are shared 
hshared = hef[which(is.na(match(husing,shared)) == FALSE)]
fshared = fef[which(is.na(match(fusing,shared)) == FALSE)]
ratios = lapply(1:length(hshared), function(x){
  hshared[[x]]$cpD1-fshared[[x]]$cpD1
})
unlist(ratios)
feff = lapply(m,function(x){
  efficiency(x,plot = FALSE) 
})
heff= lapply(hm,function(x){
  efficiency(x,plot = FALSE)
})
ratiosf = lapply(1:length(feff), function(x){
   c(heff[[x]]$cpD1-feff[[x]]$cpD1,heff[[x]]$cpD2-feff[[x]]$cpD2)
})
#####fam cpd1 vector#####
fcpd1 = vector(); fcpd1 = lapply(1:length(feff),function(x){
  fcpd1 = c(fcpd1, feff[[x]]$cpD1)
})
fcpd1 = unlist(fcpd1)
#####fam cpd2 vector#####
fcpd2 = vector(); fcpd2 = lapply(1:length(feff),function(x){
  fcpd2 = c(fcpd2, feff[[x]]$cpD2)
})
fcpd2 = unlist(fcpd2)
#####hex cpd1 vector#####
hcpd1 = vector(); hcpd1 = lapply(1:length(heff),function(x){
  hcpd1 = c(hcpd1, heff[[x]]$cpD1)
})
hcpd1 = unlist(hcpd1)
#####hex cpd2 vector#####
hcpd2 = vector(); hcpd2 = lapply(1:length(heff),function(x){
  hcpd2 = c(hcpd2, heff[[x]]$cpD2)
})
hcpd2 = unlist(hcpd2)

eff.frame = as.data.frame(ratiosf, col.names = colnames(f1017)[-1], row.names = c("cpD1rat","cpD2rat"))
eff.frame = rbind(fcpd1,fcpd2,hcpd1,hcpd2,eff.frame)
row.names(eff.frame) = c("FAM cpD1", "FAM cpD2",
              "HEX cpD1", "HEX cpD2",
              "cpD1 ratio", "cpD2 ratio")

write.csv(eff.frame, file = "Github/Stapleton-Lab/qPCR/1017ratios.csv")
