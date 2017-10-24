# Run VAMF with 10 percent subset data
## Load VAMF and useful
library(vamf)
library(useful)
## Read in subset matrix. Because duplicate row names exist, set row.names = NULL first (R will create the row names column using numbers)
S <- read.table("10percentCellSubset_allgenes.txt", header=TRUE, quote = "\"'", row.names = NULL)
## Check corners of the table
corner(S)
## Make the gene count mask for gene filtering
geneCountsSum <- rowSums(S[,2:ncol(S)]) # the 1st column still carries gene names, and need to be skipped for row sum calc
## Remove genes (rows) with no expression
S <- S[geneCountsSum != 0,]
## Keep genes only when it is detected in more than 10 cells
Sfilter <- S[rowSums(S[,2:1331] > 0) > 10, ]
## Check the dimension post filtering 
dim(Sfilter) # Here for this subset is (15665, 1331)
## Set gene names column as row names for the table, then remove the column with names from the table
rownames(Sfilter) <- Sfilter$row.names
Sfilter <- Sfilter[,2:ncol(Sfilter)]
## Check the type of table, and change its type to matrix prior processing
typeof(Sfilter) # Should be "list"
Smtx <- as.matrix(Sfilter)
typeof(Smtx) # Should be "integer"
## Calculate PCA, and export factors to file
pca <- prcomp(t(Smtx), center=TRUE, scale=TRUE)
pca_factor <- pca$x
write.table(pca_factor, file = "subset_pca_factor.txt", quote=FALSE, sep="\t")
## Calculate VAMF factors, and export to file
vamf <- vamf(Smtx,10,nrestarts=2,log2trans=TRUE)
vamf_factor <- vamf$factors
## Save vamf_factor to file
write.table(vamf_factor, file="subset_vamf_factor.txt", quote=FALSE, sep="\t")
## Combine cell metadata with the two factor file in command line (join)
## Read in joined file and plot
pcacol <- read.table("subset_pca_factor_w_meta.txt", header = FALSE, sep = '\t')
pcacol$cols <- as.numeric(as.factor(pcacol$V2))
legend.cols <- unique(pcacol$V2)
plot(pcacol[,3:4], col=pcacol$cols)
legend("topleft", legend=unique(pcacol$V2), pch=16, col=legend.cols)
## Same for VAMF
vamfcol <- read.table("subset_vamf_factor_w_meta.txt", header = FALSE, sep = '\t')
plot(vamfcol[,4:5], col = as.numeric(as.factor(vamfcol$V2)))
legend("bottomright", legend=unique(vamfcol$V2), pch=16, col=unique(vamfcol$V2))
## Calculate percentage explained in each component. Require the entire object created by prcomp()
eigs <- pca$sdev^2
perexplain <- rbind(
    SD = sqrt(eigs),
   Proportion = eigs/sum(eigs),
   Cumulative = cumsum(eigs)/sum(eigs))
perexplain
## Same for VAMF. As mentioned in the paper, used dimension learning instead
barplot(colNorms(vamf$factors),xlab="Dimension",ylab="L2 norm", main="Dimension Learning")
## Are PCs correlated with detection rate? Following VAMF example workflow. pca calculation part can be replaced by above.
### Using their defined functions rm_zero_rowcol(), pca()
rm_zero_rowcol<-function(Y){
     #remove all rows and columns containing all zeros
     Y<-Y[rowSums(Y>0)>0,] #remove rows with zeros all the way across
     Y<-Y[,colSums(Y>0)>0]
     Y
}
pca<-function(Y,L=2,center=TRUE,scale=TRUE){
  Y<-rm_zero_rowcol(Y)
   #if(scale) Y<-scale(Y)
     factors<-prcomp(as.matrix(t(Y)),center=center,scale=scale)$x
     factors<-factors[,1:L]
     colnames(factors)<-paste0("dim",1:L)
     as.data.frame(factors)
}
### Calculate first 2 principal components
factors <- pca(Smtx,2)
### Combine metadata and store in ref
ref <- pcacol[,1:2]
colnames(ref) <- c("id","batch") # The function defined requires the column name of ref as id and batch
### Calculate detection rate for each cell (column)
cens_rates_obj <- Matrix::colMeans(Smtx==0)
factors$detection_rate<-1-cens_rates_obj
### Plot! Use batch as group color label, and set transparency to 0.5 in alpha. Use gradient rainbow color palette.
ggplot(cbind(factors,ref),aes(x=detection_rate,y=dim1,color=batch))+geom_point(size=3, alpha = 0.5) + theme_classic()+ggtitle("PC1 vs. cell-specific detection")+labs(x="Detection Rate",y="Dimension 1") + theme(legend.position="top") + scale_colour_gradientn(colours=rainbow(6))
## Are VAMF dimensions correlated with detection rate?
### Plot vamf using similar work flow. Extract first 2 dimensions.
factors <- vamf_factor[,1:2]
factors$detection_rate <- 1-cens_rates_obj
ggplot(cbind(factors,ref),aes(x=detection_rate,y=dim1,color=batch))+geom_point(size=3, alpha = 0.5) + theme_classic()+ggtitle("Dim1 vs. cell-specific detection")+labs(x="Detection Rate",y="Dimension 1") + theme(legend.position="top") + scale_colour_gradientn(colours=rainbow(6))


```
[102217 version]
## Calculate PCA factors, and export to file
pca_factor <- prcomp(t(Smtx), center=TRUE, scale=TRUE)$x
write.table(pca_factor, file = "subset_pca_factor.txt", quote=FALSE, sep="\t")
## Calculate VAMF factors, and export to file
vamf_factor <- vamf(Smtx,10,nrestarts=2,log2trans=TRUE)$factors
## Save vamf_factor to file
write.table(vamf_factor, file="Group1_vamf_factor.txt", quote=FALSE, sep="\t")
## In shell, join cell metadata with the PCA matrix. Assign colors to them by using metadata column as factor
###pca_factors <- read.table("subset_pca_factor_w_meta.txt", header=FALSE, sep = "\t")
###plot(pca_factors[,2:3] , col=pca_factors$V1)
###legend("topleft", legend=levels(pca_factors$V1), pch=16, col=unique(pca_factors$V1))
## Similarly, join cell metadata with VAMF matrix
with(vamf_factors,plot(dim1,dim2, col=))
```
