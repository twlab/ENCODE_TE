#R 3.5.1
library(ggplot2)

PLS_ENCOTE<-read.table("ENCODTE.25.PLS.TE.5classes.txt",header=T)
pELS_ENCOTE<-read.table("ENCODTE.25.pELS.TE.5classes.txt",header=T)
dELS_ENCOTE<-read.table("ENCODTE.25.dELS.TE.5classes.txt",header=T)
CTCF_ENCOTE<-read.table("ENCODTE.25.CTCF.TE.5classes.txt",header=T)
DNase_ENCOTE<-read.table("ENCODTE.25.DNase.TE.5classes.txt",header=T)

ggplot(PLS_ENCOTE,aes(x=total,fill=class))+geom_histogram(position = "fill",bins=25)+theme_bw()+ggtitle("TE contribution to tissue-specific PLS sites")+xlab("number of cell types")
ggplot(pELS_ENCOTE,aes(x=total,fill=class))+geom_histogram(position = "fill",bins=25)+theme_bw()+ggtitle("TE contribution to tissue-specific pELS sites")+xlab("number of cell types")
ggplot(dELS_ENCOTE,aes(x=total,fill=class))+geom_histogram(position = "fill",bins=25)+theme_bw()+ggtitle("TE contribution to tissue-specific dELS sites")+xlab("number of cell types")
ggplot(CTCF_ENCOTE,aes(x=total,fill=class))+geom_histogram(position = "fill",bins=25)+theme_bw()+ggtitle("TE contribution to tissue-specific CTCF sites")+xlab("number of cell types")
ggplot(DNase_ENCOTE,aes(x=total,fill=class))+geom_histogram(position = "fill",bins=25)+theme_bw()+ggtitle("TE contribution to tissue-specific DNaseI sensitive regions")+xlab("number of cell types")

