#####11/9/17 qPCR Analysis#####
library(qpcR)
f1117<- read.csv(file = "C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/qpcR/1117fam2.csv", header = TRUE)
f1117[,1] = NULL
m <- modlist(x = f1117, cyc = 1, model = l4, remove = "none")
h1117<- read.csv(file = "C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/qpcR/1117hex2.csv", header = TRUE)
h1117[,1] = NULL
hm <- modlist(x = h1117, cyc = 1, model = l4)
feff = lapply(1:length(m),function(x){
  efficiency(m[[x]],plot = FALSE) 
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

eff.frame = as.data.frame(ratiosf, col.names = colnames(f1117)[-1], row.names = c("cpD1rat","cpD2rat"))
eff.frame = rbind(fcpd1,fcpd2,hcpd1,hcpd2,eff.frame)
row.names(eff.frame) = c("FAM cpD1", "FAM cpD2",
                         "HEX cpD1", "HEX cpD2",
                         "cpD1 ratio", "cpD2 ratio")

write.csv(eff.frame, file = "Github/Stapleton-Lab/qPCR/1117ratios.csv")

#####Getting specific simplex FAM results####
wanted = c(18:21,38:41,58:61,
           78:81,98:101,110:113)
wm <- modlist(x = f1117,fluo = wanted, cyc = 1, model = l4, remove = "none")
wcpd1 = lapply(wm,function(x){
  efficiency(x)$cpD1
})
wcpd1 = unlist(wcpd1)
wcpd2 = lapply(wm,function(x){
  efficiency(x)$cpD2
})
wantedcp = rbind(wcpd1,wcpd2)
rownames(wantedcp) = c("cpD1","cpD2")
colnames(wantedcp) = c("B1","B2","B3","B4",
                       "D1","D2","D3","D4",
                       "F1","F2","F3","F4",
                       "H1","H2","H3","H4",
                       "J1","J2","J3","J4",
                       "L1","L2","L3","L4")
write.csv(x = wantedcp,file = "Simplex\ Data1117.csv")
