## 091717

Start running VAMF locally. Installed VAMF in R 3.4.1 enviroment. Got an err message when running:
```
vamf_factor <- vamf(G1testmtx,10,nrestarts=2,log2trans=TRUE)$factors
```
> Error in new_CppObject_xp(fields$.module, fields$.pointer, ...) : 
>  no valid constructor available for the argument list
> trying deprecated constructor; please alert package maintainer

My test matrix is filtered Group1 matrix (only genes with more than 1 transcript in over 10 cells are kept), with 14918 rows X 1000 columns. Submitted this as issue to the developer.
