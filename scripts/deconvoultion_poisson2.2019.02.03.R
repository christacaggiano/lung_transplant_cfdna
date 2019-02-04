setwd("~/Desktop/deconvolution methylation")
library(reshape2)
library(ggplot2)

md.full <- read.csv("lung_meth.txt", sep = "\t")
# Remove sites not contributing to leukocytes vs epithelial cells.
md<- md.full[! (is.na(md.full$leuckocytes) | is.na(md.full$lung.epithelial)), ]

# Group data into dataframes 
data.methylated <- md[,grepl("\\.meth",colnames(md))]
data.methylated[is.na(data.methylated)]<-0
data.unmethylated <- md[,grepl("\\.unmeth",colnames(md))]
data.unmethylated[is.na(data.unmethylated)]<-0

all.methylation.types <-  c("leuckocytes" , "lung.epithelial" , "health.lung.cfna" 
                            , "Pulmonary.endothelial.cells", "pulmonary.fibroblasts" )
target.cell.types <- c("lung.epithelial","leuckocytes" , "Pulmonary.endothelial.cells","pulmonary.fibroblasts"  )
probs<-md[,target.cell.types]
probs[is.na(probs)]<-0

# Fit a linear model
fit <- lm (as.matrix(data.methylated) ~ probs[,1]+probs[,2]+probs[,3]+probs[,4])
x<-round(prop.table(coef(fit),2),5)[(-1),]
out <-t(x)
rownames(x)<- target.cell.types
colnames(x)<-names(data.methylated)
print(x)

# Alternatively optimize to minimize sum of squares. Error is worse and this is slower.
z<- c()
for (i in 1:ncol(data.methylated)){
  opt.fn <- function (X){
    methylated.vals <- t( X %*%t(probs))
    unmethylated.vals <- t( (X) %*%t(1-probs))
    return(sum( (methylated.vals - data.methylated[,i])^2 +
                  (unmethylated.vals - data.unmethylated[,i])^2
                , na.rm = T) * (1+sum(X<0))) # Include error term for values <0
  }
  opt <- optim( rep(.25,ncol(probs)), opt.fn)
  z<- rbind(z, opt$par)
}
z <- round(prop.table(z,margin = 1),digits = 5)

rownames(z)<-substr(colnames(data.methylated), 1,6)
colnames(z)<-target.cell.types
print (z)

###  Estimate error versus coverage using linear model ###
# This takes a while... #

md<- md.full[! (is.na(md.full$leuckocytes) | is.na(md.full$lung.epithelial)), ]
target.cell.types <- c("lung.epithelial","leuckocytes" , "Pulmonary.endothelial.cells","pulmonary.fibroblasts"  )
probs<-md[,target.cell.types]
probs[is.na(probs)]<-0
n.samples = 4
samples <- c("Epi", "Leuk", "Leuk/Stromal", "Leuk+others")
per.snp <- 0.8
# s1 is all lung, #s2 is all Immune, 

simulated.dist <- matrix(c(1,0,0,0,
                           0,1,0,0,
                           0,.5,0,.5,
                           .15,.5,.15,.2),ncol = n.samples, byrow = T)
colnames(simulated.dist) <- target.cell.types
rownames(simulated.dist) <- samples 

predicted.values <- t(simulated.dist%*%t(probs))

out <- c()
for (per.snp in 5*runif(500)){
  test.data.methylated <- apply(predicted.values, 1:2, function(X) rpois(1, per.snp*X))
  fit <- lm (test.data.methylated ~ probs[,1]+probs[,2]+probs[,3]+probs[,4])
  x<-round(prop.table(coef(fit),2),5)[(-1),]
  err <- sum(abs(t(x)-simulated.dist))/length(x)
  out<-rbind(out, cbind(per.snp, err))
}
print (out)
plot(out)
gg.df <- as.data.frame(out)
gg.df$log.snp<-log(gg.df$per.snp)
ggplot(gg.df, aes( y = err, x=per.snp))+
  geom_point(color="darkgrey")+
  geom_hline(yintercept = 0.01, color="lightgrey")+
  geom_vline(xintercept = .55, color="lightgrey")+
  geom_smooth(formula = y~log(x))+
  theme_classic(base_size = 12)+
  xlim(0,5)+
  ylab("Mean Error")+xlab("Mean Reads per SNP")

write.csv(out, file="error vs snps-lm-unif-2.csv")

### Make a pretty heatmap for top 95% of genes separating cell types ##
library(pheatmap)
probs<-md[,target.cell.types]
probs[is.na(probs)]<-0

# calculate informative regions by sum of squares vs. average.
probs$diff <- apply(probs, MARGIN=1, FUN=function (X) sum((X-mean(X))^2))

cutoff <- quantile(probs$diff,c(0.95), na.rm = T)
probs.sample <- probs[probs$diff>cutoff,]
probs.sample <- probs.sample[,c(-5)]
colnames(probs.sample)<- c("Epithelial","Leuckocytes" ,"Endothelial",  "Stromal" )
pheatmap(probs.sample[,c(-5)],show_rownames = F, treeheight_row = 0, treeheight_col = 0)

### After multiplying by starting cfDNA concentration, can graph ##

conc.data <- read.csv("deconvoluted mcg.ml2.csv", header = T)
conc.data$X<-factor(conc.data$X, levels =c("2 days"   , "2 weeks"   ,"CLAD-free","CLAD"      ) )

gg.df <- melt(conc.data)
gg.df$value <- gg.df$value/120 # because 25 ul elution from a 3 ml prep

ggplot(gg.df, aes(y = value, x=X, group=variable,color=variable))+
  geom_line()+
  geom_point(aes(color=variable))+
  theme_classic(base_size = 12)+
  ylab("Cell-free DNA (pg/ul)")+
  xlab("")+
  labs(color = "Cell type")+
  facet_grid(rows=vars(gg.df$variable!="Leukocyte"), scales = "free")
