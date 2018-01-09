install.packages("bnlearn")
library(bnlearn)
letters[1:3]
dag = model2network("[A][C][F][B|A][D|A:C][E|B:F]")
plot(dag)
fitted = bn.fit(dag, learning.test)
fitted
dag = model2network("[A][B][E][G][C|A:B][D|B][F|A:D:E:G]")
plot(dag)
fitted = bn.fit(dag, gaussian.test)
fitted
learn.net = empty.graph(names(learning.test))
modelstring(learn.net) = "[A][C][F][B|A][D|A:C][E|B:F]"
learn.net
plot(learn.net)
score(learn.net, learning.test)
#####blacklists and unsupervised learning#####
bl = matrix(c("A", "B", "B", "A"), ncol = 2, byrow = TRUE)
plot(hc(learning.test, blacklist = bl))

#####bootstrapping#####
unlist(bn.boot(learning.test, statistic = narcs,
               algorithm = "hc", R = 10))
x = list(1,3,5,7,9)
x= c(x,11)
x
class(x)
x[[1]] = matrix(1:8, ncol = 2, byrow = TRUE)

res = empty.graph(names(alarm))
modelstring(res) = paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",
                           "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]",
                           "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]",
                           "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]",
                           "[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]",
                           "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]",
                           "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")
plot(res)
