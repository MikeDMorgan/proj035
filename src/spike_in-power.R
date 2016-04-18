##########
# spike-in power analysis
##########
# DESeq2 analysis on spike-in experiment
library(DESeq2)
spike.table <- read.table("/ifs/projects/proj035/pipeline_project035_mRNA/spike_in.dir/Rngerm-genesetAll-spike3.tsv", 
                          h=T, row.names=1)
spike.design <- read.table("/ifs/projects/proj035/pipeline_project035_mRNA/design.tsv", h=T, row.names=1)
spike.design <- data.frame(spike.design[,-c(1,3)])
colnames(spike.design) <- "condition"
rownames(spike.design) <- colnames(spike.table)

spike.deseq <- DESeqDataSetFromMatrix(countData=spike.table,
                                      colData=spike.design,
                                      design= ~ condition)
spike.de <- DESeq(spike.deseq)
spike.res <- results(spike.de)
spike.df <- data.frame(spike.res)
# calculate the average success rate for each effect size observed
# bins=0, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0

bins <- seq(0, 2.0, 0.1)
power = numeric(length(bins))

for(i in 1:length(bins-1)){
  cbin = bins[i]
  nbin = bins[i+1]
  efs <- spike.df[abs(spike.df$log2FoldChange) >= cbin & abs(spike.df$log2FoldChange) < nbin,]
  sig.efs <- efs[efs$padj <= 0.05,]
  ntot = dim(efs)[1]
  sigtot = dim(sig.efs)[1]
  pwwr = sigtot/ntot
  power[i] <- pwwr
}

power[is.na(power)] <- 0.99999
pwr.df <- data.frame(bins, power)
colnames(pwr.df) <- c("log2FoldChange", "power")
pwr.df$group = "spike_in"
pwr_p <- ggplot(pwr.df, aes(x=log2FoldChange, y=power)) + 
  geom_point() +  geom_line(aes(group=group)) + theme_bw() + labs(x="log2FoldChange", y="Power") +
  geom_vline(mapping=aes(xintercept=0.6), linetype=4, colour="maroon") + 
  geom_hline(mapping=aes(yintercept=0.8), linetype=5, colour="grey")
pwr_p

# # perform same analysis for miRNAs
# mir.table <- read.table("/ifs/projects/proj035/pipeline_project035_smRNA/spike_in.dir/Rn_germ-rn5_mirbase-spike_in.tsv.gz", 
#                         h=T, row.names=1, sep="\t")
# mir.spike.design <- read.table("/ifs/projects/proj035/pipeline_project035_smRNA/design.tsv", h=T, row.names=1)
# mir.spike.design <- data.frame(mir.spike.design[,-c(1,3)])
# colnames(mir.spike.design) <- "condition"
# rownames(mir.spike.design) <- colnames(mir.spike.table)
# 
# mir.spike.deseq <- DESeqDataSetFromMatrix(countData=mir.table,
#                                           colData=mir.spike.design,
#                                           design= ~ condition)
# mir.spike.de <- DESeq(mir.spike.deseq)
# mir.spike.res <- results(mir.spike.de)
# mir.spike.df <- data.frame(mirspike.res)
# calculate the average success rate for each effect size observed
# bins=0, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0
mir.spike$padj[is.na(mir.spike$padj)] <- 1.0
mir.bins <- seq(0, 0.5, 0.05)
mir.power = numeric(length(mir.bins))

for(i in 1:length(mir.bins-1)){
  cbin = mir.bins[i]
  nbin = mir.bins[i+1]
  efs <- mir.spike[abs(mir.spike$log2FoldChange) >= cbin & abs(mir.spike$log2FoldChange) < nbin,]
  mir.sig.efs <- efs[efs$padj <= 0.05,]
  ntot = dim(efs)[1]
  sigtot = dim(mir.sig.efs)[1]
  pwwr = sigtot/ntot
  mir.power[i] <- pwwr
}

mir.power[is.na(mir.power)] <- 0.99
mir.pwr.df <- data.frame(mir.bins, mir.power)
colnames(mir.pwr.df) <- c("log2FoldChange", "power")
mir.pwr.df$group = "spike_in"
mir.pwr_p <- ggplot(mir.pwr.df, aes(x=log2FoldChange, y=power)) + 
  geom_point() +  geom_line(aes(group=group)) + theme_bw() + labs(x="log2FoldChange", y="Power") +
  geom_vline(mapping=aes(xintercept=0.26), linetype=4, colour="maroon") + 
  geom_hline(mapping=aes(yintercept=0.8), linetype=5, colour="grey")
mir.pwr_p