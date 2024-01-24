#R 3.5.1
library(ggplot2)

finalallTEinCRE.reldist <-read.table("reldist.allCRE.combine.final.txt",header=T)
ggplot(finalallTEinCRE.reldist,aes(x=reldist,y=fraction,color=type))+geom_line()+xlim(c(0,0.49))+facet_grid(TE~cCRE)+theme_bw()

