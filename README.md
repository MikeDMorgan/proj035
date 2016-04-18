# Project035 analysis pipeline
This pipeline chains together tasks required to analyse RNA-seq and small RNA-seq datasets used in project035.  It is not a general purpose pipeline, and as such some variables are hard-coded.  It contains all of the necessary custom Python and R scripts written for the project.  It relies heavily on the CGAT code collection and ruffus, the former is available here on GitHub:
https://github.com/CGATOxford/cgat

Third party tool (including R packages) requirements:
* DESeq2
* gplots
* RColorBrewer
* ggplot2
* featureCounts v1.4.6
