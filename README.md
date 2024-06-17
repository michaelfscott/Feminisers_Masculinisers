# One-factor sex determination evolves without linkage between feminising and masculinising mutations

Scripts to replicate figures and numerical analysis in DOI:10.1098/rspb.2024.0693.

## Getting Started

### Dependencies

* These scripts were created using R version 4.1.0, which will need to be [installed](https://cran.r-project.org/doc/manuals/r-patched/R-admin.html)
* The R packages used for manipulating data and creating figures are listed in scripts/functions.R. These packages can be installed using

```
install.packages(c("readr", "stringr", "dplyr", "tidyr", "ggplot2", "viridis", "ggpubr", "ggtext"))
```

### Quick Start Example

We provide an [example script](scripts/example.R) that will numerically iterate invasion by a feminiser and then a masculiniser, allowing exploration of different parameters and producing a plot like this: 

![Example of a transition from cosexuality to dioecy](./figures/example.png)

## Author

Michael F Scott

## Reference

DOI:10.1098/rspb.2024.0693
