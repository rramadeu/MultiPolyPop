# diaQTLSimulate
Wrapper R package to coordinate simulations for classic breeder mating designs in outcrossing plants. It also simulates genotypic data from GBS. This package relies on genotypic simulations from `PedigreeSim` and uses some `updog` package functions for the simulation of the genotypes and GBS data. To use `diaQTLSimulate` you need to have Java installed in your machine (https://java.com/en/download/help/download_options.html). If you use this package you also should cite `PedigreeSim` and `updog` articles, please cite them properly.

# Disclaimer
This software was originally called `PedigreeSimR` and it has been renamed for `diaQTLSimulate`. Please note that `diaQTLSimulate` and this repository are not affiliated with `PedigreeSim` software or its authors. If you are looking for the PedigreeSim software please go to https://www.wur.nl/en/show/Software-PedigreeSim.htm.

# Citation

I created `diaQTLSimulate` as an auxiliary package to do the simulations for these two manuscripts:

Amadeu, R. R., Munoz, P., Zheng, C., & Endelman, J. B. (2021). QTL Mapping in Outbred Tetraploid (and Diploid) Diallel Populations. Genetics, 219(3). [https://doi.org/10.1093/genetics/iyab124](https://doi.org/10.1093/genetics/iyab124)

Zheng, C., Amadeu, R. R., Munoz, P., & Endelman, J. B. (2021). Haplotype reconstruction in connected tetraploid F1 populations. Genetics, 219(2). [https://doi.org/10.1093/genetics/iyab106](https://doi.org/10.1093/genetics/iyab106)

As there is no proper manuscript for this package, you can cite either of those two articles as reference source. If you use this package in your research, you also should cite `PedigreeSim` and `updog` articles as follow.

# References:

[PedigreeSim](https://www.wur.nl/en/show/Software-PedigreeSim.htm):

Voorrips, R. E., & Maliepaard, C. A. (2012). The simulation of meiosis in diploid and tetraploid organisms using various genetic models. BMC bioinformatics, 13(1), 248.

[updog](https://CRAN.R-project.org/package=updog):

Gerard D., Ferrão L., Garcia A., Stephens M. (2018). Genotyping Polyploids from Messy Sequencing Data. Genetics, 210(3), 789–807. ISSN 0016-6731

Gerard D., & Ferrão L. (2020). Priors for Genotyping Polyploids. Bioinformatics, 36(6), 1795–1800. ISSN 1367-4803



# Basic usage in R
```R
## Installing this diaQTLSimulation package
devtools::install_github("rramadeu/diaQTLSimulate")
library(diaQTLSimulate)

## Creating a fake map (here you add a real map with haplotypes)
map = 1:100
haplotypes = fake_haplo(n=50,m=100,seed=1234)

## Diallel pedigree with 3 parents, no selfs, and total pop size = 90
# Creates 3 fullsib pops from 3 parents (A, B, C), AxB (n=30), AxC (n=30), BxC (n=30)
pedigree = diallel_pedigree(parents=3,popsize=100)

## Diallel pedigree with 3 parents, no selfs, and total pop size = 100
# Creates 3 fullsib pops from 3 parents (A, B, C), AxB (n=334), AxC (n=333), BxC (n=333)
# If the number of populations is not a factor of the popsize
# it creates pops with difference of 1 unit to fit the popsize specified
pedigree = diallel_pedigree(parents=3,popsize=100)

# however, if nextinteger=TRUE, it look for the next integer based on the number of pops
# creates 3 fullsib pops from 3 parents (A, B, C), AxB (n=334), AxC (n=334), BxC (n=334)
pedigree = diallel_pedigree(parents=3,popsize=100,nextinteger = TRUE)

# you can also set the total size of each cross:
pedigree = diallel_pedigree(parents=3,subpopsize=333)

## Diallel pedigree with 6 parents, no selfs, and total pop size = 100
pedigree = diallel_pedigree(parents=6,popsize=1000)

## Diallel pedigree with 6 parents, 2 selfs, and total pop size = 100
pedigree = diallel_pedigree(parents=6,selfs=2,popsize=100)

## Diallel pedigree with 6 parents and 2 selfs, and each pop (cross) with size = 100
pedigree = diallel_pedigree(parents=6,selfs=2,subpopsize=100)

## Single-Round-Robin with 6 parents and total pop size = 100
pedigree = round_pedigree(parents=6,popsize=100)

## Single-Round-Robin with 6 parents and total pop size = 100
pedigree = round_pedigree(parents=6,popsize=100,nextinteger=TRUE)

## Single-Round-Robin with 5 parents and each pop (cross) with size = 100
pedigree = round_pedigree(parents=6,subpopsize=100)

## Simulating genotypes
pedigreesimR(map,haplotypes,pedigree,ploidy=4)

## Simulating genotypes and GBS data (avg depth=60, seqerror = 0.001, allelic bias = 0.7, overdispersion = 0.005)
pedigreesimR(map,haplotypes,pedigree,ploidy=4,GBS=TRUE,GBSavgdepth = 60,GBSseq = 0.001,GBSbias = 0.7,GBSod = 0.005)

## Simulating genotypes with GBS data and doing SNP calling with updog with 2 cores considering the general model (it takes a while)
pedigreesimR(map,haplotypes,pedigree,ploidy=4,GBS=TRUE,GBSavgdepth = 60,GBSseq = 0.001,GBSbias = 0.7,GBSod = 0.005,GBSsnpcall=TRUE,GBSnc = 2)
```


# Simulation of QTL effect and genome scan internally with diaQTL
```R
setwd("~/Documents/diaQTLSimulation_Example") #choose a folder for the simulations
library(diaQTLSimulate)

## Simulating a scenario of 3 parents, autotetraploid, 200 individuals, 0.3 QTL h2
## Setting parameters
parents=3
ploidy=4
popsize=200
map = seq(0,100,.1)
QTLh2 = 0.3

## Simulation of QTL effects

## Simulating a scenario of 3 parents, autotetraploid, 200 individuals, 0.3 QTL h2
## Setting parameters
parents=3 #number of parents in the diallel
ploidy=4 #2 or 4
popsize=200 #total population size
map = seq(0,100,.1) #one chromosome at time
QTLh2 = 0.3 #if interested in FPR investigation, set QTLh2 to 0

## Simulation
haplotypes = fake_haplo(n=50,m=length(map),seed=1234)
pedigree = diallel_pedigree(parents=parents,popsize=popsize)
pedigreesimR(map,haplotypes,pedigree,ploidy=ploidy,workingfolder = "diaQTLSimulation_files")
QTLsim(parents=parents,
       ploidy=ploidy,
       workingfolder="PedigreeSimR_files",
       QTLmarker=NULL,
       QTLh2=QTLh2,
       run_diaQTL=FALSE) #FALSE if just interested in creating diaQTL inputs
## "diaQTLSimulation_files" folder is now populated with the files necessary for PolyOrigin sofware (`polyorigin` prefix) and for diaQTL software (`QTL` prefix). 

## To perform genome scan with diaQTL:
setwd("~/Documents/diaQTLSimulation_Example/diaQTLSimulation_files")
library(diaQTL)
data1 <- read_data(genofile="QTLgeno.csv",
                   pedfile="QTLped.csv",
                   phenofile="QTLpheno.csv",
                   ploidy=4,
                   dominance=1)                   
set_params(data1,trait="pheno1")
ans <- scan1(data=data1,trait="pheno1",n.core=1,dominance = 1)
scan1_summary(ans)

## Or: read_data, set_params, scan1, scan1_summary functions are encapsulated within QTLsim function when run_diaQTL=TRUE
setwd("~/Documents/diaQTLSimulation_Example/")
QTLsim(parents=parents,
       ploidy=ploidy,
       workingfolder="diaQTLSimulation_files",
       QTLmarker=NULL,
       QTLh2=QTLh2,
       run_diaQTL=TRUE) 
 ```
