fam = read.csv(file = url("https://raw.githubusercontent.com/tbillman/Stapleton-Lab/master/vQTL%20Random%20and%20Family/data/tidied/Family.csv"), header = TRUE)

chr1.nums = unlist(lapply(4:length(fam[1,]),function(x){
 a = NULL
  if(as.character(fam[1,x]) == "1"){
    a = c(a,x)
 } else{a = a}
 return(a)
}))
chr1.names = colnames(fam)[chr1.nums]
index = 1:length(chr1.nums)
chr1.pos = unlist(lapply(1:length(chr1.nums),function(x){
  as.numeric(as.character(fam[2,chr1.nums[x]]))
}))
chr1df = cbind(index, chr1.names, chr1.pos)
write.table(chr1df, file = "Github/Stapleton-Lab/mQTL\ Random\ and\ Family/mapchr1.txt",quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
