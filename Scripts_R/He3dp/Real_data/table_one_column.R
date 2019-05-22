# Table 1 column 
require(xtable)
NA_I_latex <- read.csv("NA_I.csv")



xl <- NA_I_latex 
df1 = data.frame(xl[1:46,c(1,3,2,4,5)])


x.big <- xtable(df1[,c(1,2,5)], type = "latex",display=c("e","g","E","E"),digits=4,caption= "Case I")




short <- xtable(df1[,c(1,2,5)], type = "latex",display=c("e","g","E","g"),digits=4,caption= "Case I")


print(short, include.rownames=FALSE)

write.csv(df1[,c(1,2,5)],"NV_one_column.csv",row.names = F)

write.table(short,"NV_one_column_latex.tex",row.names = F)

