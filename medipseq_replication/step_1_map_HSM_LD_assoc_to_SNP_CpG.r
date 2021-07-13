#Step 1: Inputting the results from Bell et al. (PMID: 29295990) and map the HSM-LD blocks to CpGs and SNPs from our analysis
#Analyst: Rory Wilson
#Date: 02.07.2021

rm(list=ls())

library(dplyr); library(parallel)

################## Loading necessary data #########################

#loading the 450K annotation file
an450 <- data.table::fread(file="illumina_450K_annot_pruned.txt",sep="\t") %>% as.data.frame(stringsAsFactors=F) 

#cpgs actually used in our analysis (ie, after qc)
load("kora_cpgs.RData")
an450 <- filter(an450,IlmnID %in% cpgs)

#read in the medip results: a file received from the author
#Renaming and new variables:
#V1 is the chromosome, 
#V2, V3 are the start and end of the HSM peaks, 
#V4 chromosome of LD blocks
#V5, V6 start and end of LD blocks
#6th are the SNPs in the LD blocks
#hsm peaks and ld blocks are found on the same chromosomes
#hsm ID will just be based on rownames

medip_res <- read.csv("mergeFinal_hsm_LDblock_GWASSNPs.txt",
					sep="\t",header=FALSE,stringsAsFactors = FALSE) %>% 
					rename(chr=V1,hsm_low=V2,hsm_up=V3,LD_low=V5,LD_up=V6,snps_in_LD_block=V7)        %>% 
					select(-V4)                                                          %>% 
					mutate(hsm_length=hsm_up - hsm_low,LD_block_length = LD_up-LD_low,hsm_id = make.names(rownames(.)))

#Basic descriptive info:
dim(medip_res)                #[1] 7184    9
summary(medip_res$hsm_length)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  500.0   500.0   750.0   816.3  1000.0 32500.0
summary(medip_res$LD_block_length)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  60  118910  240445  379656  450922 9184438

#as stated by the authors, some HSM peaks overlap a couple LD blocks: here number of unique HSM peaks
paste0(medip_res$chr,"_",medip_res$hsm_low,"_",medip_res$hsm_up) %>% unique %>% length    #[1] 7173

#QC: Checking that all HSM peaks fall within their LD blocks: #three possibilities
p_1 <- with(medip_res,(hsm_up  >= LD_low) & (hsm_up <=LD_up))
p_2 <- with(medip_res,(hsm_low >= LD_low) & (hsm_low <=LD_up))
p_3 <- with(medip_res,(hsm_low <= LD_low) & (hsm_up >=LD_up)) 
all(p_1 | p_2 | p_3) %>% stopifnot

####
rm(p_1,p_2,p_3)
####


################## Checking the overlap of CpG sites we tested and HSM peaks #########################

#dividing each into chromosomes
medip_by_chr <- lapply(unique(medip_res$chr),function(x) filter(medip_res,chr==x))
an450_by_chr <- lapply(unique(medip_res$chr),function(x) filter(an450,chr==x))
names(medip_by_chr) <- names(an450_by_chr) <- unique(medip_res$chr)


#function to check for an overlap of HSM peak and a vector of CpG site positions: return T/F for a given HSM peak
check_medip_overlap <- function(pos_vec,hsm_low,hsm_up) {
	any((pos_vec >= hsm_low) & (pos_vec <= hsm_up))
}

#function to get the actual overlap: ie getting the cpgs sites themselves that fall within the HSM peak
get_overlap <- function(annot_df,hsm_low,hsm_up) {
	annot <- annot_df[(annot_df$MAPINFO >= hsm_low) & (annot_df$MAPINFO <= hsm_up),,drop=FALSE]
	return(annot)
}


#Main function: Getting the CpG sites that are contained within the HSM peaks, chromosome by chromosome
res_comp <- mcmapply(function(annot_chr,medip_chr) {

	medip_range <- medip_chr[,c("hsm_low","hsm_up")]
	is.numeric(medip_range[,"hsm_low"])  %>% stopifnot
	is.numeric(medip_range[,"hsm_up"])   %>% stopifnot

	overlaps <- apply(medip_range,1,function(x) {
		get_overlap(annot_chr,hsm_low=x["hsm_low"],hsm_up=x["hsm_up"])
	})
	names(overlaps) <- medip_chr$hsm_id

	return(overlaps)

},an450_by_chr,medip_by_chr,mc.cores=22)


#getting the HSM ids as part of the results, and binding together
for (i in 1:length(res_comp)) {
	for (j in 1:length(res_comp[[i]])) {
		if (nrow(res_comp[[i]][[j]])>0) res_comp[[i]][[j]]$hsm_id <- names(res_comp[[i]])[j]
} }
all_overlap <- lapply(res_comp,function(x) do.call(rbind,x)) %>% do.call(rbind,.)  

#quality control based on rownames of all_overlap
strsplit(rownames(all_overlap),"\\.") %>% sapply(function(x) x[[2]]) %>% identical(all_overlap$hsm_id) %>% stopifnot

#merging the results with the original MeDIP results
all_res <- merge(medip_res,all_overlap,by="hsm_id")
dim(all_res) #[1] 487  12


################## quality control of CpG results to this point #########################

all(all_res$chr.x==all_res$chr.y)       %>% stopifnot
all(all_res$hsm_low <= all_res$MAPINFO) %>% stopifnot
all(all_res$hsm_up >= all_res$MAPINFO)  %>% stopifnot

#one CpG is within two LD blocks:
cpg_dup <- filter(all_res,duplicated(IlmnID))$IlmnID
filter(all_res,IlmnID %in% cpg_dup)
#  hsm_id chr.x  hsm_low   hsm_up   LD_low    LD_up            snps_in_LD_block hsm_length
#1  X2637 chr16 54114500 54115250 54114528 54114824 rs16953002_1st        750
#2  X2638 chr16 54114500 54115250 54114823 54208972 rs16953002_2nd        750
#  LD_block_length     IlmnID  MAPINFO chr.y
#1             296 cg06434738 54115019 chr16
#2           94149 cg06434738 54115019 chr16

all_res <- mutate(all_res,chr=chr.x) %>% select(-chr.x,-chr.y)


####
rm(i,j,get_overlap,check_medip_overlap,res_comp,all_overlap,cpg_dup,an450,medip_by_chr,cpgs,an450_by_chr)
####

################## only the SNPs from our analysis #########################

#snps used in our analyses
load("kora_snps.RData")

#splitting each SNP list (separated by ";", sometimes with the tags "_1st" or "_2nd"), per row
split_snps <- lapply(all_res$snps_in_LD_block,function(x) {
  strsplit(split=";",x=x) %>% unlist %>% gsub("_1st|_2nd","",.) }
)

#qc:
(length(split_snps)==nrow(all_res)) %>% stopifnot

#expanding each row of all_res to one row per SNP
expand_res <- lapply(1:length(split_snps),function(i) {
  df_total <- c()
  df_row   <- all_res[i,,drop=F]
  for (j in 1:length(split_snps[[i]])) df_total <- rbind(df_total,df_row)
  df_total$precise_snp <- split_snps[[i]]
  return(df_total)
}) %>% do.call(rbind,.)

#filtering to only those SNPs we tested specifically
assoc_to_test <- filter(expand_res,precise_snp %in% snps)
dim(assoc_to_test) #[1] 845  12


################## saving #########################

save(assoc_to_test,file="medip_associations_to_test.RData")



