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
)

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
              'G1','G9','G13','G17','H5','I1','I9','I13','I17','J5','K1','K9','L5','P1','P9'
  
)
shared = setdiff(husing,setdiff(husing,fusing))
which(husing) = shared[1]
ratios = lapply(1:length(hef), function(x){
  hef[[x]]$cpD1/fef[[x]]$cpD1
})
