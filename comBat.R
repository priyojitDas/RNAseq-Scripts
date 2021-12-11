library(sva)

filelist <- list.files(".",pattern="*.tsv")
fbase_list <- c()
batch_list <- c()
s_list <- c()
mf_list <- c()

for (i in 1:length(filelist)){
	print(i)
	fbase <- substr(filelist[i],start=1,stop=nchar(filelist[i])-4)
	fmeta <- unlist(strsplit(fbase,"_"))
	fbase_list <- c(fbase_list,fbase)
	batch_list <- c(batch_list,fmeta[2])
	if (fmeta[3] == "P"){
		status <- 1
	}
	else{
		status <- 0
	}
	s_list <- c(s_list,status)
	if (fmeta[4] == "M"){
                mfstatus <- 0
        }
        else{
                mfstatus <- 1
        }
	mf_list <- c(mf_list,mfstatus)
	tdf <- read.delim(filelist[i],sep='\t',header=FALSE)
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
#adjusted_counts <- ComBat_seq(df, batch=batch_list)
#adjusted_counts <- ComBat_seq(df, batch=batch_list, group=s_list)
adjusted_counts <- ComBat_seq(df, batch=batch_list, group=NULL, covar_mod=cbind(s_list,mf_list))

rownames(adjusted_counts) <- genenames
colnames(adjusted_counts) <- fbase_list

write.csv(adjusted_counts,"ComBatProgeria_o_BSMF.csv")

