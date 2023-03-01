rm(list=ls())
library(Cairo)
library(Hmisc)
library(plotrix)
library(ggplot2)

traitname <- "g2mvuse2"
traits.main1 <- "/restricted/projectnb/cvdmrx/cacscnt2/"
trait.results <-paste("/restricted/projectnb/cvdmrx/cacscnt2/gen2/rst")
## go to each rst directory cat *csv> "traitmodel"rst.csv ##
load("/restricted/projectnb/cvdmrx/cacscnt2/gen2/rst/g2mvuse2.RData")
 colnames(g2mvuse2) <- c("IlmnID","beta","se","p","totaln")

## load annotation data
load("/restricted/projectnb/fhs-methylation/cleaned/scripts/shuo/HumanMethylation450_15017482_simple.Rdata")
post.dir <- "/restricted/projectnb/cvdmrx/cacscnt2/gen2/rst/" 
png_file_name <- "g2mvuse2"
## merge results with annotation data
dd2   <- merge(g2mvuse2,annotM,by="IlmnID")

## calculate lambda
print("Lambda:")
lambda.print <- median(qchisq(dd2$p,df=1,lower.tail=F),na.rm=T)/0.455
print(lambda.print)

## write out annotated datafile
write.table(dd2, file=paste(post.dir,"g2mvuse2705.csv"),sep=",",row.names=F,quote=F)
b    <- dd2[order(dd2$p),]

trait <- "g2mvuse2"
MAIN  <- "g2mvuse2"

out4 <- b
CHR <- as.character(out4$CHR)
CHR[CHR=="X"] <- "23"
CHR[CHR=="Y"] <- "24"
out4$CHR  <- as.integer(CHR)

chr <- c(1:23)

# QQ plot PNG
QQplot_png <- function(p, png_file_name)
   {
   obs <- -log(p,10)
   N <- length(obs) ## number of p-values

   ## create the null distribution (-log10 of the uniform)
   null <- -log(1:N/N,10)
   MAX <- max(c(obs,null))

   ## the jth order statistic from an Uniform(0,1) sample
   ## has a beta(j,n-j+1) distribution
   ## (Casella & Berger, 2002, 2nd edition, pg 230, Duxbury)
   c95 <- qbeta(0.95,1:N,N-(1:N)+1)
   c05 <- qbeta(0.05,1:N,N-(1:N)+1)

   lc95 <- -log(c95,10)
   lc05 <- -log(c05,10)

   # Axis labels
   x.lab <- expression(Expected~~-log[10](italic(p-value)),cex.axis=.4)
   y.lab <- expression(Observed~~-log[10](italic(p-value)),cex.axis=.4)

   png(png_file_name, width=600, height=600, type="cairo")

   ## plot confidence lines
   plot(null, lc95, ylim=c(0,MAX), xlim=c(0,max(null)), type="l", axes=FALSE, xlab="", ylab="", col="gray")
   par(new=T)
   plot(null, lc05, ylim=c(0,MAX), xlim=c(0,max(null)), type="l", axes=FALSE, xlab="", ylab="", col="gray")
#   polygon(c(null,rev(null)), c(lc05, rev(lc95)), border="white")
   lines(c(0,max(null)),c(0, max(null)), lwd=2)
   ## add diagonal
   par(new=T)
   ## add qqplot
   qqplot(null,obs, ylim=c(0,MAX),xlim=c(0,max(null)), main=traitname[2], pch=16, xlab=x.lab, ylab=y.lab, cex=1)
   text(1,MAX-1,substitute(lambda==x, list(x=round(lambda.print,2))),cex=2)
  dev.off()
   }

QQplot_png(g2mvuse2$p, paste(post.dir,"/QQ_plot_","gen2ascac2",".png",sep=""))



chr <- c(1:23)

#############
b2 <- b[b[,"p"] < 1e-05,]
b2 <- b2[order(b2[,"p"]),]

write.table(b2, file =paste(post.dir,"/","g2mvuse2","_top_cpgs.csv", sep=""),sep=",",row.names=F,quote=F)

# Manhattan plot PNG
MH_plot_PNG <- function(A, MH_file_name = "MH_plot.png", mytitle = "", color1 = "gray50", color2 = "gray25", pch.type = 16)
   {
   A <- A[order(A$CHR, A$MAPINFO),]
   
   A$log10p <- -log10(A$p)
   n.snps <- nrow(A)
   ymax <- max(A$log10p)+1
   A$abs_position <- 1:n.snps

   png(MH_file_name, width=1500, height=900, type="cairo")

   par(mar=c(3,5,1,1))
   plot(A$abs_position, A$log10p, type="n", yaxt="n", xaxt="n", xlab="", ylab="", main=mytitle, xlim=c(0,n.snps), ylim=c(0,ymax))
   chr <-  1:23
   middle.points <- numeric(22)
   for (i in 1:23)
      {
      idx <- A$CHR==i
      points(A$abs_position[idx], A$log10p[idx], col=ifelse(i%%2==0, color1, color2), pch=pch.type, cex=0.8)   
      middle.points[i] <- median(A$abs_position[idx])
      }
   axis(side=2, at=seq(from=0, to=ymax, by=1), labels=seq(from=0, to=ymax, by=1), tick=T, cex.axis=1, las=1)
   axis(side=1, at=middle.points, labels=as.character(c(1:23)), tick=F, cex.axis=1.1, las=1, pos=0)
   mtext(text="Chromosomes", side=1, cex=1.25, line=1.15)
   mtext(text=expression(-log[10]* (p-value)), side=2, line=3, cex=1.5)
   abline(h=-log10(1*10^-10),lty=2,col="black")
   dev.off()
   }

MH_plot_PNG(out4, paste(post.dir,"/MH_plot_","g2mvuse2",".png", sep=""), mytitle="g2mvuse2")



## dev.off()
