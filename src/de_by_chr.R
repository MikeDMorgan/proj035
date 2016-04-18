# loop over all chromosomes, counting the number of unadjusted p < 0.05

de_by_chr <- list()
for(chr in unique(merge.de$CHR)){
  pvals <- (merge.de[merge.de$CHR == chr,])$sig
  de_by_chr[[chr]] <- (table(pvals)[2])/length(pvals)
}

de_chr_df <- t(data.frame(de_by_chr))
chr.reps <- gsub("X$", "XX", rownames(de_chr_df))
rownames(de_chr_df) <- gsub("^X", "", chr.reps)
de_chr_df <- data.frame(t(data.frame(de_by_chr)), rownames(de_chr_df))


colnames(de_chr_df) <- c("prop_DE", "CHR")


mir_by_chr <- list()
for(chr in unique(merge.mir$CHR)){
  pvals <- (merge.mir[merge.mir$CHR == chr,])$sig
  mir_by_chr[[chr]] <- table(pvals)[2]/length(pvals)
}

mir_chr_df <- t(data.frame(mir_by_chr))
mir.reps <- gsub("X$", "XX", rownames(mir_chr_df))
rownames(mir_chr_df) <- gsub("^X", "", mir.reps)
mir_chr_df <- data.frame(t(data.frame(mir_by_chr)), rownames(mir_chr_df))
