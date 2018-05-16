
require(xtable)
NA_I_latex <- read.csv("NA_I.csv")
NA_II_latex <- read.csv("NA_II.csv")


xl <- NA_I_latex 
df1 = data.frame(xl[1:30,c(1,3,2,4,5)])
df2 = data.frame(xl[31:60,c(1,3,2,4,5)])
dfFull = data.frame(df1,df2)
x.big <- xtable(dfFull, type = "latex",display=c("e","g","E","E","E","g",
                                                 "g","E","E","E","g"),digits=4,caption= "Case I")
print(x.big, include.rownames=FALSE)

write.csv(gg2data,"NV_case_I.csv",row.names = F)
