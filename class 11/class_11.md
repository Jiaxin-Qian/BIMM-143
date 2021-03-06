class 11: Structutal Bioinformatics pt.1
================

\#\#PDB Statistics

Here we inspect the types of structure in the main database for 3D
biomolecular data

> Q1. determine the percentage of structures solved by X-Ray and
> Electron Microscopy. Also can you determine what proportion of
> structures are protein?

``` r
#read the file 
stats <- read.csv("Data Export Summary.csv", row.names = 1)
#calculate
#calculate the percentage for each method 
ans <- stats$Total/sum(stats$Total) * 100
names(ans) <- rownames(stats)
round(ans, 2)
```

    ##               X-Ray                 NMR Electron Microscopy 
    ##               88.95                8.04                2.72 
    ##               Other        Multi Method 
    ##                0.19                0.10

``` r
protein <- round(sum(stats$Proteins)/sum(stats$Total) * 100, 2)
protein
```

    ## [1] 92.69

> Q2: Type HIV in the PDB website search box on the home page and
> determine how many HIV-1 protease structures are in the current PDB?
