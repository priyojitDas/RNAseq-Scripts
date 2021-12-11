library(dplyr)
library(DESeq2)

extract.deseq <- function(dds,seed,alpha,condition,contrast1,contrast2) {
	
	# Example parameters and variables
	# dds: eg1.wt.ko.de
	# seed: "eg1.wt.ko.de"
	# alpha: 0.05
	# condition: "celltype"
	# contrast1: "KO"
	# contrast2: "WT" 
	
	x <- results(dds, alpha=alpha, contrast=c(condition, contrast1, contrast2));

	summary(x);
	
	x <- data.frame(x) %>%
	dplyr::mutate(genename=factor(rownames(x))) %>%
	dplyr::select(
		genename,
		lfc=log2FoldChange,
		padj=padj
	) %>%
	dplyr::mutate(bin=
		as.factor(
			ifelse(is.na(lfc) | is.na(padj), "no",
				ifelse(padj < alpha & lfc > 0, "up",
					ifelse(padj < alpha & lfc < 0, "dn", "no")
				)
			)
		)
	) %>%
	dplyr::mutate(bin=factor(bin, levels=c("up","dn","no")))
	
	names(x) <- c("genename", paste0(seed,".lfc"), paste0(seed,".padj"), paste0(seed,".3bin"))

	return(x);
}

flist <- read.delim('age13+_O.txt',header=FALSE)$V1

for (i in 1:length(flist)){
        print(i)
	tdf <- read.delim(as.character(flist[i]),sep='\t',header=FALSE)
        if (i == 1){
                df <- tdf
        }
        else{
                df <- cbind(df,tdf$V2)
        }
}
df <- df[1:26364,]
genenames <- df$V1
df <- df[,-1]
rownames(df) <- genenames

colData <- data.frame(
	cell.type=c(rep("Normal",15), rep("HGPS",7)),
	row.names=flist
)
# 7 13
# 6 7
# 14 13
# 15 13
# 14 7
# 15 7
rna.norm.hgps.de <- DESeqDataSetFromMatrix(
        countData=df,      
        colData=colData, 
        design= ~ cell.type
)

rna.norm.hgps.de <- DESeq(rna.norm.hgps.de)

rna.norm.hgps.de05 <- extract.deseq(rna.norm.hgps.de, "rna.norm.hgps.de05", 0.01, "cell.type", "HGPS", "Normal")

write.csv(rna.norm.hgps.de05,"ProgeriaDEG_age_13+_O_de01.csv")
