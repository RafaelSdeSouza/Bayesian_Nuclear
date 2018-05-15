
require(xtable)
gg2data <- read.csv("NV_case_I.csv")


df1 = data.frame(gg2data[1:30,])
df2 = data.frame(gg2data[31:60,])
dfFull = data.frame(df1,df2)

x.big <- xtable(gg2data, type = "latex",display=c("e","g","E","E","E","g"
                                                 ),digits=4,caption= "Case I")
print(x.big, include.rownames=FALSE)



x.big <- xtable(dfFull, type = "latex",display=c("e","g","E","E","E","g",
                                                 "g","E","E","E","g"),digits=4,caption= "Case I")
print(x.big, include.rownames=FALSE)

write.csv(gg2data,"NV_case_I.csv",row.names = F)
