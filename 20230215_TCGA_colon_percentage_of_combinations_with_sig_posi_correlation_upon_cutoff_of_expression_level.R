# This script is to draw line graph about percentage of combinations with significant positive correlation upon cutoff miRNA and transcript expression  
# 2023/02/15 made

# import table about mean of miRNA and transcript expression
# this table is located at ""
setwd("C:/Rdata/20230206_TCGA_colon_confirmation_of_expression_levels_of_transcript_and_miRNA")
mean.df <-read.table("TCGA_colon_confirmation_of_expression_level_of_transcript_and_miRNA.txt",sep="\t",header = T,stringsAsFactors = F)

# delete unnecessary rows and columns, and do log2
mean.df <-mean.df[,c(1:5,9,13)]
mean.df <-mean.df[!is.na(mean.df[,3]),]
mean.df[,6] <-log2(mean.df[,6])
mean.df[,7] <-log2(mean.df[,7])

# set cutoff value about miRNA or transcript expression
m.cutoff <-seq(-10,10,2)
t.cutoff <-seq(-10,10,2)

# make table for cutoff of miRNA expression
m.sm <-as.data.frame(matrix(nrow = length(m.cutoff),ncol = 4))
colnames(m.sm) <-c("cutoff","number_of_sig_posi","total","percentage_of_sig_posi")

# calculate percentage upon cutoff of miRNA expression
for (i in 1:length(m.cutoff)) {
  m.cutoff.df <-mean.df[mean.df[,6]>=m.cutoff[i],]
  m.sm[i,1] <-m.cutoff[i]
  m.sm[i,2] <-nrow(m.cutoff.df[m.cutoff.df[,3]>0&m.cutoff.df[,4]<0.05,])
  m.sm[i,3] <-nrow(m.cutoff.df)
  m.sm[i,4] <-m.sm[i,2]/m.sm[i,3]*100
}

# make table for cutoff of transcript expression
t.sm <-as.data.frame(matrix(nrow = length(m.cutoff),ncol = 4))
colnames(t.sm) <-c("cutoff","number_of_sig_posi","total","percentage_of_sig_posi")

# calculate percentage upon cutoff of transcript expression
for (i in 1:length(t.cutoff)) {
  t.cutoff.df <-mean.df[mean.df[,7]>=t.cutoff[i],]
  t.sm[i,1] <-t.cutoff[i]
  t.sm[i,2] <-nrow(t.cutoff.df[t.cutoff.df[,3]>0&t.cutoff.df[,4]<0.05,])
  t.sm[i,3] <-nrow(t.cutoff.df)
  t.sm[i,4] <-t.sm[i,2]/t.sm[i,3]*100
  }

# make new directory
setwd("C:/Rdata")
dir.create("20230215_TCGA_colon_percentage_of_combinations_with_sig_posi_correlation_upon_cutoff_of_expression_level")
setwd("C:/Rdata/20230215_TCGA_colon_percentage_of_combinations_with_sig_posi_correlation_upon_cutoff_of_expression_level")

# draw each line graph
pdf("line_graph_of_percentage_of_combinations_with_sig_posi_correlation_upon_cutoff_of_miRNA_exprssion.pdf")
plot(m.sm[,1],m.sm[,4],type="b",ylim = c(0,100),xlab="cutoff for mean of miRNA expression",ylab="percentage of combinations with sig posi correlation",
    pch=19)
text(m.sm[,1], m.sm[,4], paste0(m.sm[,2],"/",m.sm[,3]), adj=c(0.5,-0.5),cex = 0.6)
dev.off()

pdf("line_graph_of_percentage_of_combinations_with_sig_posi_correlation_upon_cutoff_of_transcript_exprssion.pdf")
plot(t.sm[,1],t.sm[,4],type="b",ylim = c(0,100),xlab="cutoff for mean of transcript expression",ylab="percentage of combinations with sig posi correlation",
     pch=19)
text(t.sm[,1], t.sm[,4], paste0(t.sm[,2],"/",t.sm[,3]), adj=c(0.5,-0.5),cex = 0.6)
dev.off()

###

# import correspondence table between TCGA colon transcriptome bam and TCGA colon miRNA quantification file
# this table is located at "https://github.com/Ryosuke-Hirota/20221124_make_table_of_TCGA_transcriptome_bam_correspond_to_TCGA_miRNA_expression_data"
setwd("C:/Rdata/20221124_make_table_of_TCGA_transcriptome_bam_correspond_to_TCGA_miRNA_expression_data")
cor.table <-read.table("correspondence_table_between_TCGA_colon_transcriptome_bam_and_miRNA_qunatification_file.txt",sep="\t",header = T,stringsAsFactors = F)

# import list of transcripts that intersect with miRNAs in gencode v36
# this list is located at "https://github.com/Ryosuke-Hirota/20230110_TCGA_colon_transcriptome_bam_correlation_between_transcript_and_miRNA"
setwd("C:/Rdata")
primir.list <-read.table("TCGA_hg38_transcript_intersect_with_miRNA.txt",sep="\t",header = F,stringsAsFactors = F)
primir.list[,2] <-primir.list[,2]-1
primir.list <-primir.list[primir.list[,2]!=primir.list[,8]&primir.list[,3]!=primir.list[,9],]

# remove duplicated gene names and make list regarding miRNA name and gene name
mir_transcipt <-primir.list[,c(4,11)]
mir_transcipt <-subset(mir_transcipt,!duplicated(mir_transcipt))
mir_transcipt <-mir_transcipt[order(mir_transcipt[,1]),]
primir.list <-primir.list[,c(4,11,10)]
transcripts <-unique(primir.list[,3])

# list TCGA colon transcript quantification files
# these files are located at "\\fsw-q02\okamura-lab\20221006_TCGA_colon_salmon_quant_transcriptome"
setwd("C:/Rdata/20221006_TCGA_colon_salmon_transcriptome_quantification")
transcript.quant <-list.files(path = "C:/Rdata/20221006_TCGA_colon_salmon_transcriptome_quantification",pattern = ".txt")
t.file.id <-gsub("_quant.txt","",transcript.quant)

# make table summarized TCGA colon transcript quantification
for (i in 1:length(transcript.quant)){
  transcript.file <-read.table(transcript.quant[i],sep = "\t",header = T,stringsAsFactors = F)
  t <-match(transcripts,transcript.file[,1])
  transcript.file <-transcript.file[t,]
  transcript.file <-transcript.file[,c(1,4)]
  colnames(transcript.file)[2] <-transcript.quant[i]
  if(i==1){
    transcript.quant.table <-transcript.file
  }else{
    transcript.quant.table <-merge(transcript.quant.table,transcript.file,by="Name")
  }}

# import table summarized miRNA quantification files
# this table is located at ""
setwd("C:/Rdata/20230105_TCGA_colon_miRNA_quantification")
miRNA.quant.table <-read.table("table_of_TCGA_colon_miRNA_quantifications.txt",sep="\t",header = T,stringsAsFactors = F,check.names = F)

# make empty summary
sm <-as.data.frame(matrix(nrow = nrow(mir_transcipt),ncol = 3))
colnames(sm) <-c("miRNA","transcript","number_of_sample_without_expression")

# investigate quantile and mean of each expression level
for (i in 1:nrow(mir_transcipt)){
  # extract expression level of a certain miRNA 
  miRNA.df <-miRNA.quant.table[miRNA.quant.table[,1]==mir_transcipt[i,1],]
  miRNA.df <-as.data.frame(t(miRNA.df),stringsAsFactors = F)
  m.cor <-match(rownames(miRNA.df),cor.table[,5])
  miRNA.df[,2] <-rownames(miRNA.df)
  miRNA.df[,3] <-cor.table[m.cor,2]
  colnames(miRNA.df) <-c(mir_transcipt[i,1],"miRNA_file_name","transcript_file_name")
  rownames(miRNA.df) <-NULL
  miRNA.df <-miRNA.df[-1,]
  
  # extract expression level of a certain transcript
  transcript <-primir.list[primir.list[,1]==mir_transcipt[i,1]&primir.list[,2]==mir_transcipt[i,2],3]
  t <-match(transcript,transcript.quant.table[,1])
  transcript.df <-transcript.quant.table[t,]
  transcript.df <-transcript.df[,-1]
  transcript.df[1,] <-apply(transcript.df, 2, sum)
  transcript.df <-transcript.df[-2,] 
  transcript.df <-as.data.frame(t(transcript.df),stringsAsFactors = F)
  t.cor <-match(t.file.id,cor.table[,1])
  transcript.df[,2] <-cor.table[t.cor,2]
  colnames(transcript.df) <-c(mir_transcipt[i,2],"transcript_file_name")
  rownames(transcript.df) <-NULL
  
  # merge expression levels of miRNA and transcript
  mt.df <-merge(miRNA.df,transcript.df,by="transcript_file_name")
  mt.df <-subset(mt.df,!is.na(mt.df[,1]))
  mt.df <-mt.df[,c(4,2,1,3)]
  mt.df[,2] <-as.numeric(mt.df[,2])
  
  # extract row without miRNA or trasnscript expression
  zero.df <-mt.df[mt.df[,1]==0|mt.df[,2]==0,]
  
  # write summary
  sm[i,1] <-mir_transcipt[i,1]
  sm[i,2] <-mir_transcipt[i,2]
  sm[i,3] <-nrow(zero.df)
}

# import table about quantile and mean of miRNA/transcript expression
# this table is located at ""
setwd("C:/Rdata/20230206_TCGA_colon_confirmation_of_expression_levels_of_transcript_and_miRNA")
mean.df <-read.table("TCGA_colon_confirmation_of_expression_level_of_transcript_and_miRNA.txt",sep="\t",header = T,stringsAsFactors = F)

# edit table and do log2
mean.df <-mean.df[,c(1:5,9,13)]
mean.df <-mean.df[!is.na(mean.df[,3]),]
mean.df[,6] <-log2(mean.df[,6])
mean.df[,7] <-log2(mean.df[,7])

# merge table and summary
zero.mean.df <-merge(mean.df,sm,by=c("miRNA","transcript"))

# output merged table
setwd("C:/Rdata/20230215_TCGA_colon_percentage_of_combinations_with_sig_posi_correlation_upon_cutoff_of_expression_level")
write.table(zero.mean.df,"TCGA_colon_table_about_number_of_sample_without_expression.txt",sep="\t",row.names = F,quote = F)

# set cutoff for number of sample without expression
cutoff <-seq(0,280,20)

# make table for
sm <-as.data.frame(matrix(nrow = length(cutoff),ncol = 4))

# calculate percentage
for (i in 1:length(cutoff)) {
  cutoff.df <-zero.mean.df[zero.mean.df[,8]<=cutoff[i],]
  sm[i,1] <-cutoff[i]
  sm[i,2] <-nrow(cutoff.df[cutoff.df[,3]>0&cutoff.df[,4]<0.05,])
  sm[i,3] <-nrow(cutoff.df)
  sm[i,4] <-sm[i,2]/sm[i,3]*100
}

# draw line graph
pdf("line_graph_about_percentage_of_combinations_with_sig_posi_correlation_upon_cutoff_of_number_of_sample.pdf")
plot(sm[,1],sm[,4],xlab="cutoff for number of sample without expression",ylab="percentage of combinations with sig posi correlation",
     pch=19,type="b")
text(sm[,1],sm[,4], paste0(sm[,2],"/",sm[,3]), adj=c(0.5,-0.5),cex = 0.6)
dev.off()
