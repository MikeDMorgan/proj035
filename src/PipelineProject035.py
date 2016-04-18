import CGAT.Experiment as E
import pandas as pd
import pandas.rpy.common as com
from rpy2.robjects import pandas2ri as py2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import r as R
import rpy2.robjects as ro
import CGATPipelines.Pipeline as P
import itertools
import numpy as np


@P.cluster_runnable
def deseqAnalysis(counts_table,
                  design,
                  reference,
                  outfile):
    '''
    Perform differential expression analysis using DESeq2
    '''

    design_df = pd.read_table(design, sep="\t",
                              header=0, index_col=0)
    counts_df = pd.read_table(counts_table, sep="\t",
                              header=0, index_col=0, compression="gzip")

    E.info("setting up counts table")
    py2ri.activate()
    r_design = py2ri.py2ri_pandasdataframe(design_df)
    r_counts = py2ri.py2ri_pandasdataframe(counts_df)

    R.assign("design", r_design)
    R.assign("counts", r_counts)

    R('''sink(file="/dev/null")''')
    E.info("loading required R packages")
    R('''suppressPackageStartupMessages(library(DESeq2))''')
    R('''suppressPackageStartupMessages(library(gplots))''')
    R('''suppressPackageStartupMessages(library(RColorBrewer))''')
    R('''suppressPackageStartupMessages(library(ggplot2))''')

    R('''notZero <- counts[rowMeans(counts) > 1,]''')
    R('''dds <- DESeqDataSetFromMatrix(countData=notZero,'''
      '''colData=design, design=~group)''')
    E.info("performing differential expression testing")
    R('''de <- DESeq(dds, fitType="parametric")''')
    R('''res <- results(de)''')

    E.info("generating MA plot")
    # generate MAplots
    R('''png("images.dir/%s-MAplot.png", height=480, width=480)''' % reference)
    R('''plotMA(res, alpha=0.05)''')
    R('''dev.off()''')

    E.info("performing variance stabilising transformation")
    R('''vst <- data.frame(getVarianceStabilizedData(de))''')

    E.info("clustering samples and plotting heatmap")
    R('''cors <- cor(vst)''')
    R('''hmcol <- colorRampPalette(brewer.pal(9, "PuOr"))(100)''')
    R('''png("images.dir/%s-sample_clustering-heatmap.png", height=480, '''
      '''width=480)''' % reference)
    R('''heatmap.2(as.matrix(cors), col=hmcol, trace="none", '''
      '''breaks=seq(0, 1, 0.01), margins=c(10,10), cexRow=0.8,'''
      '''cexCol=0.8)''')
    R('''dev.off()''')

    E.info("performing principal components analysis")
    R('''pca <- prcomp(data.frame(t(vst)), scale=T, centre=T)''')
    R('''pcs <- data.frame(pca$x)''')
    R('''pcs$condition <- as.factor(design$group)''')
    R('''p_pca <- ggplot(pcs, aes(x=PC1, y=PC2, colour=condition)) + '''
      '''geom_point(size=6)''')
    R('''png("images.dir/%s-PCA_pc1-pc2.png", height=480, '''
      '''width=480)''' % reference)
    R('''print(p_pca)''')
    R('''dev.off()''')

    E.info("writing table of results")
    R('''res.df <- data.frame(res)''')
    ('''sink(file=NULL)''')
    out_df = com.load_data("res.df")
    out_df.to_csv(outfile, sep="\t", index_label="gene_id")


def enumerateCounts(counts_file, design_file, 
                    bin_size, max_bin):
    '''
    enumerate a counts file, then shuffle the rows and bin into
    log fold change range

    '''
    if counts_file.endswith("gz"):
        compression = "gzip"
    else:
        compression = None

    counts_df = pd.read_table(counts_file, header=0, index_col=0, sep="\t",
                              compression=compression)
    design_df = pd.read_table(design_file, header=0, index_col=0, sep="\t")

    n_genes = len(counts_df)
    genes = counts_df.index
    gene_combs = itertools.combinations(genes, 2)
    log_fold_changes = pd.DataFrame(columns=genes, index=genes,
                                    dtype=np.float64)

    # control diets are first 4 samples, hfd are last 4
    controls = counts_df.iloc[:,:4]
    hfd = counts_df.iloc[:,4:]

    bins = range(int(-max_bin*100), int(max_bin*100), int(bin_size*100))
    bins = [x/100.0 for x in bins]

    idx = 0
    for gene1, gene2 in gene_combs:
        cntrl = controls.loc[gene1]
        test = hfd.loc[gene2]
        try:
            res = np.mean(test)/np.mean(cntrl)
        except ZeroDivisionError:
            res = 2.0
        log_fold_changes.loc[gene1, gene2] = res

    log_fold_changes = log_fold_changes.fillna(0.0)
    log_fold_changes = log_fold_changes.apply(np.nan_to_num, axis=0)
    
    mean_filter = lambda x: np.mean(abs(x)) > 1.0
    cols = log_fold_changes.apply(mean_filter, axis=1)
    rows = log_fold_changes.apply(mean_filter, axis=0)

    col_df = log_fold_changes.loc[cols,:]
    row_df = col_df.loc[:,rows]
    
    return row_df
