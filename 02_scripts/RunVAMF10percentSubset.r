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
## Calculate PCA factors, and export to file
pca_factor <- prcomp(t(Smtx), center=TRUE, scale=TRUE)$x
write.table(pca_factor, file = "subset_pca_factor.txt", quote=FALSE, sep="\t")
## Calculate VAMF factors, and export to file
vamf_factor <- vamf(Smtx,10,nrestarts=2,log2trans=TRUE)$factors
## Save vamf_factor to file
write.table(vamf_factor, file="Group1_vamf_factor.txt", quote=FALSE, sep="\t")
## Read in cell metadata, and join this with the PCA matrix? Assign colors to them by using metadata column as factor
## referring to plot(sample(1:10,20,TRUE), col=factor(sample(letters[1:3],20,TRUE)))
## Plot separately
plot(pca_factors[,1:2])
with(vamf_factors,plot(dim1,dim2))
