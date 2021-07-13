#Step 2: Check how many of the HSM-LD block associations we replicate with our data, and comparing this to the matched background
#Analyst: Rory Wilson
#Date: 02.07.2021

rm(list=ls())

library(dplyr)
library(ggplot2)

################## Background #########################

#Bell et al. had associations between
#LD blocks (more precisely, the SNPs within them) and HSMs contained within the LD blocks

#Testing replication in KORA data: 
#SNPs: matched precisely betweeen Bell et al. and KORA 
#Methylation: CpG sites from KORA mapped to HSMs from Bell et al.

#If a HSM-LD block contains multiple SNPs/multiple CpGs, we average the p-values as per Bell et al.

#Note that an LD block can be associated with multiple HSMs
#but HSMs can only be associated with one LD block 
#(except if it falls on the border of an LD block - occurs once)
#So a unique HSM counts as a unique association
##############################

#Loading the associations to test and the KORA results 
load("medip_associations_to_test.RData")
tested_res <- read.csv("medip_replication_results.txt",sep="\t",header=TRUE,stringsAsFactors = F)

dim(tested_res)    #[1] 845   7
dim(assoc_to_test) #[1] 845  12

#beginning the process of merging the two data-frames
assoc_to_test <- mutate(assoc_to_test,pair=paste0(precise_snp,"_",IlmnID))
tested_res    <- mutate(tested_res,pair=paste0(snp,"_",cpg))

######
#quality control: one HSM (CpG) overlaps two LD blocks (with one SNP)
#this CpG-SNP pair is in both dataframes

duplicated(assoc_to_test$pair) %>% sum #[1] 1

dup_pair <- assoc_to_test$pair[duplicated(assoc_to_test$pair)]
filter(assoc_to_test,pair==dup_pair)
#hsm_id   chr  hsm_low   hsm_up   LD_low    LD_up snps_hsm_assoc hsm_length LD_block_length     IlmnID  cpg_pos
#1  X2637 chr16 54114500 54115250 54114528 54114824 rs16953002_1st        750             296 cg06434738 54115019
#2  X2638 chr16 54114500 54115250 54114823 54208972 rs16953002_2nd        750           94149 cg06434738 54115019
#snp_in_LD_block  snp_pos                  pair
#1      rs16953002 54114824 cg06434738_rs16953002
#2      rs16953002 54114824 cg06434738_rs16953002

tested_res[duplicated(tested_res$pair),]
#cpg        snp       beta          se     pvalue          r2 nobs
#222 cg06434738 rs16953002 0.00400644 0.002062682 0.05226678 0.002300748 1638
#pair
#222 cg06434738_rs16953002
##### 

#putting them in the same order and merging them
tested_res    <- arrange(tested_res,pair)
assoc_to_test <- arrange(assoc_to_test,pair)
identical(tested_res$pair,assoc_to_test$pair) %>% stopifnot
tested_res$pair <- NULL

res <- cbind(assoc_to_test,tested_res) %>% arrange(hsm_id)

#handling the one duplicated result by adding a "v2" to the second
length(unique(res$pair)) #[1] 844
res$pair[duplicated(res$pair)] <- paste0(res$pair[duplicated(res$pair)],"_v2")

#qc on the dataframe
identical(res$IlmnID,res$cpg)      %>% stopifnot
identical(res$precise_snp,res$snp) %>% stopifnot

res <- select(res,-precise_snp,-IlmnID)

####
rm(tested_res,assoc_to_test,dup_pair)
####

################## HSM-SNP pairs we replicate #########################

#the number of hsms
unique(res$hsm_id) %>% length

#If there is more than one SNP/one CpG per HSM-LD block pair:
#1) we use the lowest p-value from our associations 
#as average p-value approach is unfair to background:
#background does not take into account any correlation structure, so unassociated results swamp out associations
#2) however, we also look at the results using the average p, as per Bell et al.

#getting the two types per HS;
hsm_pvalues_realpairs <- group_by(res,hsm_id) %>% summarize(avg_pval=mean(pvalue),low_pval=min(pvalue))


############# HSM-SNP pairs replicated using matched background #########################

#the matched background is made up of matched 100 CpG-SNP pairs

bckgrd <- read.csv(file="medip_matched_background_results.txt",sep="\t",header=T,stringsAsFactors = F)

#quality control: seeing that there are 100 matched background pairs per tested (original) pair
length(unique(res$pair)) #[1] 845
dim(res)                 #[1] 845  18
sapply(res$pair,function(x) sum(bckgrd$orig_pair==x)) %>% range #[1] 100 100

#the random matchings
random_matches <- lapply(c(1:100),function(i) filter(bckgrd,num==i))

#quality control 
sapply(random_matches,nrow) %>% range #[1] 845 845
sapply(random_matches,function(x) setequal(x$orig_pair,res$pair)) %>% all %>% stopifnot

#match of the snp-cpgs and hsm_id
random_matches <- lapply(random_matches, function(rm) {
  rm <- merge(rm,select(res,hsm_id,pair),by.x="orig_pair",by.y="pair")
  stopifnot(nrow(rm)==845)
  return(rm) })

#getting the average and lowest p-values for HSM within each random match
hsm_results_per_rm <- lapply(random_matches, function(rm) {
  group_by(rm,hsm_id) %>% summarize(avg_pval=mean(pvalue),low_pval=min(pvalue)) })

####
rm(res,random_matches,bckgrd)
####

#a function to take a cut-off, a pvalue type (low or avg), and the real associations and matched associations and return
compare_to_bg <- function(ctoff,pval_type="low_pval",ma=hsm_results_per_rm,ra=hsm_pvalues_realpairs) {
  
  real_sig_regions <- sum(ra[,pval_type] < ctoff)
  
  rm_sig_regions   <- sapply(ma,function(x) sum(x[,pval_type] < ctoff))
  
  pval             <- sum(rm_sig_regions >= real_sig_regions)/100
  
  return(list(real_sig_regions=real_sig_regions,rm_sig_regions=rm_sig_regions,pval=pval))
}

#performing the results for various significance cut-offs
cutoffs <- c(0.05,0.01,0.05/nrow(hsm_pvalues_realpairs),10^-14)

#using the low pvalue method and average pvalue method
res_per_cutoff_low <- lapply(cutoffs,compare_to_bg )
res_per_cutoff_avg <- lapply(cutoffs,compare_to_bg,pval_type="avg_pval")
names(res_per_cutoff_low) <- names(res_per_cutoff_avg) <- c("p05","p01","bonf_hsm","bonf_gw")

#####
#some summary statistics for the various p-level cutoffs, low p-value method
sapply(res_per_cutoff_low,function(x) x$real_sig_regions)       
sapply(res_per_cutoff_low,function(x) mean(x$rm_sig_regions))   
sapply(res_per_cutoff_low,function(x) median(x$rm_sig_regions)) 
sapply(res_per_cutoff_low,function(x) x$pval)

#enrichment (based on median)
sapply(res_per_cutoff_low,function(x) x$real_sig_regions/median(x$rm_sig_regions))
#[1] 1.218543 1.319672 1.377778 1.812500

#####
#some summary statistics for the various p-level cutoffs, avg p-value method
sapply(res_per_cutoff_avg,function(x) x$real_sig_regions)       
sapply(res_per_cutoff_avg,function(x) mean(x$rm_sig_regions))   
sapply(res_per_cutoff_avg,function(x) median(x$rm_sig_regions)) 
sapply(res_per_cutoff_avg,function(x) x$pval)

#enrichment (based on median)
sapply(res_per_cutoff_avg,function(x) x$real_sig_regions/median(x$rm_sig_regions))

####
rm(compare_to_bg,res_per_cutoff_avg)
####

############# Plotting results #########################

#the various p-value thresholds
levs <- sapply(cutoffs,function(x) paste0("P-cutoff: ",x))
levs[3] <- "P-cutoff: 0.00015" #manually changing the thrid one 0.05/328

#making the data frame of the results
plot_df <- data.frame("P.cutoff" = factor(as.vector(sapply(levs,function(x) rep(x,100))),levels=levs),
                     nums=as.vector(sapply(res_per_cutoff_low,function(x) x$rm_sig_regions)))

#the line of actual significant hits
df_line <- data.frame("P.cutoff" = factor(levs,levels=levs),line_pos = sapply(res_per_cutoff_low,function(x) x$real_sig_regions))

#the plot:
plt <- ggplot(plot_df, aes(x=nums)) + geom_histogram() +   theme_bw() + 
  ylab("Frequency of random samples") +
  xlab("Number of significant HSM-LD block pairs") + 
  geom_vline(data = df_line, aes(xintercept=line_pos,color="Actual sig. hits")) + facet_grid(P.cutoff ~ .) +
  scale_color_manual(values = c("Actual sig. hits" = "red")) + theme(legend.title = element_blank(),legend.position="top")

plt

#plot 3 looks like there is some with greater numbers of significant hits: this is an artefact of binning

pdf("Supplementary_Figure_MeDIP.pdf")
plt
dev.off()



######################################################
######################################################














