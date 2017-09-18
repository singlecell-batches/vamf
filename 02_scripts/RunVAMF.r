# requirement: R >= 3.4
# Install Rstan first
Sys.setenv(MAKEFLAGS = "-j4") 
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies=TRUE)
# Restart R and test with:
fx <- inline::cxxfunction( signature(x = "integer", y = "numeric" ) , '
     return ScalarReal( INTEGER(x)[0] * REAL(y)[0] ) ;
' )
fx( 2L, 5 ) # should be 10
# Install VAMF
install.packages("devtools")
devtools::install_github("willtownes/vamf")
# Load vamf
library(vamf)
# Install useful
install.packages(useful)
library(useful)

# Read in Group1 matrix
G1 <- read.table("Group1.matrix",header=FALSE,sep="\t")

# Check the matrix
corner(G1)
# Set Emsembl gene ID column as row names of the matrix
rownames(G1) <- G1$V1
# Remove Ensembl ID and gene name columns
G1 <- G1[,3:ncol(G1)]
# Check the right hand corner
corner(right(G1))
# Remove rows (genes) with no expression
G1 <- G1[rowSums(G1) != 0, ] 

# More stringent filtering (at least one transcript in more than 10 cells):
 G1test <- G1[rowSums(G1 > 0) > 10, ]

# Convert array to matrix
typeof(G1test) # it should be "list"
G1testmtx <- as.matrix(G1test)
typeof(G1testmtx) # now it will be "integer"

# Do Principal Component Analysis
pca_factor <- prcomp(t(G1testmtx), center=TRUE, scale=TRUE)$x

# Save pca_factor to file
write.table(pca_factor, file = "Group1_pca_factor.txt", quote=FALSE) # if quote = TRUE, the name labels will have quote marks.

# Run VAMF
vamf_factor <- vamf(G1testmtx,10,nrestarts=2,log2trans=TRUE)$factors

# Save vamf_factor to file
write.table(vamf_factor, file="Group1_vamf_factor.txt", quote=FALSE)

# Plot separately.
