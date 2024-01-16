#install.packages("forestploter")
#install.packages("grid")

library(grid)
library(forestploter)

setwd("G:\\Others\\Program\\MR")

mydata=read.table("Forest.txt",header = T,sep = "\t")


mydata$pval=ifelse(mydata$pval<0.001, "<0.001", sprintf("%.4f", mydata$pval))
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- sprintf("%.4f (%.4f - %.4f)",
                                   mydata$or, mydata$or_lci95, 
                                   mydata$or_uci95)

forest(mydata[,c(1:3,7:8)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.5,
       ci_column =4 ,
       ref_line = 1,
       xlim = c(0.75, 1.2))


tm1 <- forest_theme(core=list(fg_params=list(hjust = 0.9, x = 0.9),
                              bg_params=list(fill = c("#edf8e9", "#c7e9c0", "#a1d99b"))),
                    colhead=list(fg_params=list(hjust=0.5, x=0.5)))
tm2 <- forest_theme(core=list(fg_params=list(hjust=c(1, 0, 0, 0.5),
                                            x=c(0.9, 0.1, 0, 0.5)),
                             bg_params=list(fill = c("#f6eff7", "#d0d1e6", "#a6bddb", "#67a9cf"))),
                   colhead=list(fg_params=list(hjust=c(1, 0, 0, 0, 0.5),
                                               x=c(0.9, 0.1, 0, 0, 0.5))))
tm3=forest_theme(ci_Theight = 0.2)


forest(mydata[,c(1:2,6:7)],
            est = mydata$or,
            lower =mydata$or_lci95, 
            upper = mydata$or_uci95,
            sizes =0.3,
            ci_column =3 ,
            ref_line = 1,
            xlim = c(0.05, 1.3),
            theme = tm3)