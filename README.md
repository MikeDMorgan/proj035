# Project035 analysis pipeline
This pipeline chains together tasks required to analyse RNA-seq and small RNA-seq datasets used in project035.  It is not a general purpose pipeline, and as such some variables are hard-coded.  It contains the necessary custom Python and R scripts written for the project.  It relies heavily on the CGAT code collection and ruffus, the former is available here on GitHub:
https://github.com/CGATOxford/cgat

Third party tool (including R packages) requirements:
* DESeq2
* gplots
* RColorBrewer
* ggplot2
* featureCounts v1.4.6

Prior read QC and mapping was carried out using the CGATPipelines, pipeline_readqc.py and pipeline_mapping.py found in the CGATPipelines repo:
https://github.com/CGATOxford/CGATPipelines

###Input
The pipeline expects compressed alignment files (.bam), with files named according to the convention:
  `condition-tissue-replicate.file_suffix`

It also requires a gtf file of annotations of interest.
