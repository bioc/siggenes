# Copyright (c) 2002 Holger Schwender

# If there are missing values, this function will replace these with the rowwise mean.

# X: a matrix

na.replace<-function(X){
    for(i in 1:nrow(X))
        X[i,  ] <- replace(X[i,  ], which(is.na(X[i,  ])), mean(X[i,  ], na.rm = TRUE))
    return(X)
}
