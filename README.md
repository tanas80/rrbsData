
rrbsData
========

[![GitHub release](https://img.shields.io/github/release/tanas80/rrbsData.svg)](https://github.com/tanas80/rrbsData/releases)


# Introduction 

*rrbsData* is an [R](http://en.wikipedia.org/wiki/R_%28programming_language%29) project
for DNA methylation analysis and feature selection with S4 class structure. The package is designed to deal with 
[bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) output files.

## Current Features

 * Loads methylation calling data based on data.frame sample description
 * Easy R-style filtering object with indexes of specimens or CpG-dinucleotides
 * Two layered matrix data design (B-value and coverage).
 * Some reports
## Installation
In R console:
```r
install.packages("BiocManager")
library(devtools)
install_github("tanas80/rrbsData", build_vignettes=FALSE, 
               dependencies=TRUE, ref="init", repos=BiocManager::repositories())
```
-------
# How to Use

See sample.R file.

Type in R console:
```r
source("https://github.com/tanas80/rrbsData/raw/master/samples/sample.R")
```
or
```console
$ wget https://github.com/tanas80/rrbsData/raw/master/samples/sample.R
```
in Unix console to get source code of sample.R

-------
# Citing rrbsData
If you used rrbsData please cite:

 * Alexander S Tanas, Vladimir O Sigin, Alexey I Kalinkin, Nikolai V Litviakov,
Elena M Slonimskaya, Marina K Ibragimova, Mona A Frolova, Ekaterina O
Ignatova, Olga A Simonova, Ekaterina B Kuznetsova, Tatiana V Kekeeva,
Sergey S Larin, Elena V Poddubskaya, Ivan D Trotsenko, Viktoria V
Rudenko, Kristina O Karandasheva, Kseniya D Petrova, Irina V Deryusheva,
Polina V Kazantseva, Artem V Doroshenko, Natalia A. Tarabanovskaya, Galina
G Chesnokova, Marina I Sekacheva, Marina V Nemtsova, Vera L Izhevskaya,
Sergey I Kutsev, Dmitry V Zaletaev, and Vladimir V Strelnikov 
*[Genome-wide methylotyping resolves breast cancer epigenetic heterogeneity and suggests novel therapeutic perspectives](http://lab.epigenetic.ru/)* // _Epigenomics_. - _(2019)_. - 0(0). - p.000-000.
