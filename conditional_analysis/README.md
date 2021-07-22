# Conditional analysis and Pruning

Code was executed in R version 3.2.0, gtool/0.7.5, plink/1.90

Input files: 
- sigtab.RData: table of cosmopolitan pairs results (SNP-CpG pairs with P<1E-14); 
header: pair; snp; snp.chr;snp.pos; A1; A2; eaf; cpg; cpg.chr; cpg.pos; beta.disco(very); se.disco; p.disco; beta.repl(ication); se.repl; p.repl; beta.comb(ined); se.comb; p.comb; pop.disco; pop.repli

- beta.RData: consists of "cpgs.beta" (vector of all cpgs on array) and "samples.beta" (vector of all samples included)
- betaQN_[cohort].RData: matrix (n by m) of quantile-normalised beta (methylation) levels for all cpgs on array (n) for all samples in cohort (m)
- [cohort]snps.RData: matrix of N snps by 5 columns in cohort (chr, rsid, BP, A1, A2, index)
- [cohort]gen.bin: gen (genotype) file of snps listed in [cohort]snps.RData in cohort 


```{r}
require(meta)
### 1. Genome-wide association
# For each CpG, identify the complete list of SNP-CpG pairs reaching P<10-14
# load table of significant SNP-CpG pairs (P<10-14)
load('sigtab.RData')
sigtab=cosmo
cpg=as.character(unique(sigtab$cpg)[tid])
rsigtab=sigtab[sigtab$cpg==cpg,]

#get data for conditional analyses; possible to split up and run in parallel  
cohorts=c('IA317', 'IA610', 'OmniX', 'OmniExome', 'nfbc66', 'nfbc86', 'KF4') 
columns.gen=c(987, 2748, 1788, 4062, 2196, 1542,4968)
names(columns.gen)=cohorts
for(cohort in cohorts){
	print(cohort)
	#get methylation data
	if(cohort %in% c('IA317', 'IA610', 'OmniX', 'OmniExome')){
		load('beta.RData')
		to.read = file('beta.bin', "rb")
		row.number=which(cpgs.beta==cpg)
		columns.meth=length(samples.beta)
		bits.skip=8*(row.number)*columns.meth   
		seek(to.read, origin='start', where=bits.skip)
		beta = readBin(to.read, double(), n=columns.meth)
		close(to.read)
		names(beta)=samples.beta 
		assign(paste(cohort, '.beta', sep=''), beta)
	}else{
		load(paste('betaQN_', cohort, '.RData', sep=''))	
		to.read = file(paste('betaQN_', cohort, '.bin', sep=''), "rb")
		row.number=which(cpgs.beta==cpg)
		columns.meth=length(samples.beta)
		bits.skip=8*(row.number)*columns.meth    
		seek(to.read, origin='start', where=bits.skip)
		beta = readBin(to.read, double(), n=columns.meth)
		close(to.read)
		assign(paste(cohort, '.beta', sep=''), beta)
	}

	#genotypes; converts to dosage
	to.read = file(paste('/project/lolipop_b/METHQTL2/STUFF/PRUNE/COND2/DATA/', cohort, '_gen.bin', sep=''), "rb") 
	load(paste(cohort, '_snps.RData', sep='')) 
	snp.anno=snps[snps$rsid %in% rsigtab$snp,]
	rownames(snp.anno)=as.character(snp.anno$rsid)
	snp.anno=snp.anno[as.character(rsigtab$snp),]
	rownames(snp.anno)=as.character(rsigtab$snp)
	
        for(s in 1:nrow(snp.anno)){
	        dos=rep(NA,columns.gen[cohort]/3)
		if(!is.na(snp.anno$index[s])){
			row.number=as.numeric(snp.anno$index[s])
			bits.skip=8*((row.number-1)*columns.gen[cohort])
			seek(to.read, origin='start', where=bits.skip)
			geno.raw = readBin(to.read, double(), n=columns.gen[cohort])
	                i=1
	                h=1
	                while(i<length(geno.raw)){
	                        dos[h]=geno.raw[i+1]+(2* geno.raw[i+2])
	                        i=i+3
	                        h=h+1
	                }
		}
                if(exists('dosage')){
                        dosage=rbind(dosage, dos)
                }else{
                        dosage=dos
                }
        }
	alleles=snp.anno[,c('A1', 'A2')]
	alleles[,1]=as.character(alleles[,1])
	alleles[,2]=as.character(alleles[,2])
	rownames(alleles)=rownames(snp.anno)
	samples =read.table(paste(cohort, '_pheno.txt', sep=''), header=T, sep=' ')
	dosage=matrix(dosage, nrow=nrow(snp.anno))
	colnames(dosage) = samples$ID_1
	rownames(dosage)=rownames(snp.anno)
	if(!(exists('ref.alleles'))){
	        ref.alleles=alleles
	}else{
	        alleles=alleles[rownames(ref.alleles),]
	        for(s in 1:nrow(alleles)){
			if(!is.na(alleles[s,1])){
		                if((ref.alleles[s,1]!= alleles[s,1]) | (ref.alleles[s,2]!= alleles[s,2])){
		                        if((ref.alleles[s,1]== alleles[s,2]) & (ref.alleles[s,2]== alleles[s,1])){
		                                dosage[s,]=2-dosage[s,]
		                        }else{
		                                print(cbind(snp.anno[s,], ref.alleles[s,], alleles[s,]))
		                        }
				}
	                }
	        }
	}
	if(cohort=='OmniExome'){
		colnames(dosage)[1350:1354] =c('12236','12248','4687','6256','7309')
		colnames(dosage)[856:857] =c('4121', '7747')
	}
	assign(paste(cohort, '.dosage', sep=''), dosage)
	rm(dosage)
	rm(alleles)
	close(to.read)
}
```

## Preliminary LD pruning to exclude SNPs that are almost identical (and therefore resistant to conditional analysis)
```{r}
rsigtab$maf.snp=ifelse(rsigtab$eaf<0.5, rsigtab$eaf, 1-rsigtab$eaf)
if(nrow(rsigtab)>1){
	rsq.thres=0.99
	proxies=as.character(rsigtab[with(rsigtab, order(abs(z), maf.snp,snp, decreasing=T)),'snp'])  #sort by z,maf,rsid
	cohorts=c('KF4', 'OmniExome', 'IA610', 'OmniX', 'IA317', 'nfbc86', 'nfbc66')
	for(cohort in cohorts){
		sentinels=NA
		dosage=eval(parse(text=paste(cohort, '.dosage', sep='')))
		rsq=cor(t(dosage[proxies,]), use="pairwise.complete.obs")^2
		while(length(proxies)>0){
		        p=proxies[1]
		        rsq.p=rsq[p,]
		        rsq.p=rsq.p[setdiff(names(rsq.p),p)]
		        idents=names(rsq.p[rsq.p>rsq.thres])
		        proxies=setdiff(proxies, c(p, idents))
		        sentinels=c(sentinels, p)
		}
		proxies=sentinels[2:length(sentinels)]
		if(length(proxies)==1){break}
	}
	sentinels=sentinels[2:length(sentinels)]
	proxies=sentinels
}else{
	proxies=as.character(rsigtab$snp)
}
print(paste('Preliminary Pruning: ', nrow(rsigtab), '->', length(proxies), sep=''))
```

## Conditional association testing 
```{r}
nproxies=length(proxies)
cohorts=c('IA317', 'IA610', 'OmniX', 'OmniExome', 'nfbc66', 'nfbc86', 'KF4') 
thres=1E-14
sentinels=NA
#counter=1
while(length(proxies)>0){
	sentinels = na.omit(sentinels)	
	for(cohort in cohorts){
		if(cohort %in% c('IA317', 'IA610', 'OmniX', 'OmniExome')){
			samples =read.table(paste(cohort, '_pheno.txt', sep=''), header=T, sep=' ')
			if('CHD' %in% colnames(samples)){samples$CHD=as.factor(samples$CHD)}
			colnames(samples)[3:ncol(samples)]=paste(colnames(samples)[3:ncol(samples)], '_genet', sep='')
			load('phe.RData')
			phe=merge(phe, samples, by.x='sample', by.y='ID_1')
			rownames(phe)=phe$sample
		}else{
			phe =read.table(paste(cohort, '_pheno.txt', sep=''), header=T, sep=' ')			
			rownames(phe)=phe$ID_1
		}
		dosage=eval(parse(text=paste(cohort, '.dosage', sep='')))
		phe=phe[colnames(dosage),]
		beta=eval(parse(text=paste(cohort, '.beta', sep='')))
		beta.sel=beta[as.character(phe$ID)]
		if(cohort=='IA317'){
			meth.covars=paste(paste('phe$', c('Age','CD8T','CD4T','NK','Bcell','Mono','Gran'),collapse=' + ', sep=''), paste('phe$PC', seq(1,10,1), '_cp', collapse=' + ', sep=''), sep=' + ')
		}else if(cohort %in% c('nfbc66', 'nfbc86')){
			meth.covars=paste(paste('phe$', c('Sex','CD8T','CD4T','NK','Bcell','Mono','Gran'),collapse=' + ', sep=''), paste('phe$PC', seq(1,10,1), '_cp', collapse=' + ', sep=''), sep=' + ')
		}else if(cohort=='KF4'){
			meth.covars=paste(paste('phe$', c('Sex','Age','CD8T','CD4T','NK','Bcell','Mono'),collapse=' + ', sep=''), paste('phe$PC', seq(1,10,1), '_cp', collapse=' + ', sep=''), sep=' + ')
		}else{
	        	meth.covars=paste(paste('phe$', c('Sex','Age','CD8T','CD4T','NK','Bcell','Mono','Gran'),collapse=' + ', sep=''), paste('phe$PC', seq(1,10,1), '_cp', collapse=' + ', sep=''), sep=' + ')
		}
		sentinel.covars=paste('as.vector(dosage[as.character(sentinels[', seq(1, length(sentinels)), ']),])',sep='')
		if(length(sentinels)>0){
			model <- as.formula(paste("beta.sel ~ as.vector(dosage[as.character(n),])", paste(paste(sentinel.covars, collapse='+'), meth.covars, sep='+'), sep='+'))
		}else{
			model <- as.formula(paste("beta.sel ~ as.vector(dosage[as.character(n),])", meth.covars, sep='+'))
		}
		results=rep(NA,5)
		for(n in proxies){
			tryCatch({fit=lm(model)}, error = function(error) {return(NA)})
                        if(!exists("fit")){
                        	coefs=rep(NA,4)
				nsamp=NA
			}else{
	        		coefs=summary(fit)$coefficients[2,]
	       			nsamp=nobs(fit)
				rm(fit)
			}
	        	results=rbind(results, c(coefs, nsamp))
		}
		results=results[2:nrow(results),]
		results=matrix(results,nrow=length(proxies))
		rownames(results)=paste(proxies, cpg, sep='_')
		assign(paste(cohort, '.results', sep=''), results)
		if(length(sentinels)==0){print(results)}
	}

	#meta-analysis
	effects=eval(parse(text=paste(cohorts[1],'.results', sep='')))[,1]
	results=ifelse(effects==0,min(abs(effects)),effects)
	for(c in 2:length(cohorts)){
		effects= eval(parse(text=paste(cohorts[c],'.results', sep='')))[,1]
		effects=ifelse(effects==0,min(abs(effects)),effects)
		results=cbind(results, effects)
	}
	for(c in 1:length(cohorts)){results=cbind(results, eval(parse(text=paste(cohorts[c],'.results', sep='')))[,2])}
	rownames(results)=rownames(eval(parse(text=paste(cohorts[1],'.results', sep=''))))
	colnames(results)=c(paste(cohorts, '.effect', sep=''), paste(cohorts, '.se', sep=''))
	meta=matrix(nrow=nrow(results), ncol=4)
	for(i in 1:nrow(results)){
	        meta[i,1:4]=unlist(summary(metagen(results[i,1:length(cohorts)], results[i,(length(cohorts)+1):ncol(results)], comb.random=F, prediction=F))$fixed)[c(1,2,5,6)]
	}
	meta[,4]= 2*pnorm(-abs(meta[,3]))
	rownames(meta)=rownames(IA317.results)
	colnames(meta)=c('effect', 'se', 'z', 'p')
	meta=cbind(snp=proxies, meta)
	meta=data.frame(meta)
	meta$p=as.numeric(as.character(meta$p))
	meta$z=as.numeric(as.character(meta$z))
	proxies=as.character(meta[meta$p<thres,'snp'])
	#print(meta[with(meta, order(abs(z), snp, decreasing=T)),])
	#print(sort(as.vector(meta[,5])))
	if(length(proxies)>0){
		sentinel=as.character(meta[with(meta, order(abs(z), snp, decreasing=T)),'snp'])[1]
		#print(sentinel)
		sentinels=c(sentinels, sentinel)
		proxies=setdiff(proxies, sentinels)
	}
	#assign(paste('meta_', counter, sep=''),meta)
	#counter=counter+1
}

print(paste('Pruned SNPs: ', nproxies, '->', length(sentinels), sep=''))
save(nfbc66.dosage, nfbc86.dosage, KF4.dosage, meth.covars, IA317.dosage, IA610.dosage, OmniX.dosage,OmniExome.dosage, beta, rsigtab, sentinels, ref.alleles, cpg, file=paste('res',format(tid, scientific=F),'.RData', sep=''))
save(sentinels, cpg, file=paste('sentinels',format(tid, scientific=F),'.RData', sep=''))

#if split up into parallel runs, combine back results (final.RData)

# make weighted R2 files
load('final.RData')
cond.snps=as.character(unique(rsigtab$snp))
write.table(cond.snps, file='cond.snps', quote=F, row.names=F, col.names=F)
```


## Input files
- [cohort].gen.gz: gen file for each cohort

Remove list of indels if applicable/available (using indels.txt) 
Converts via PLINK from .gen to .ped and .map format 
Applies a threshold of 0.98 for calling genotypes for calculating missing data
Generates LD statistics report based on list of snps from previous section (cond.snps) 
```
module load gtool/0.7.5
#remove indels if applicable
gtool -S --g /GENOS/nfbc66.gen.gz --exclusion indels.txt --s nfbc66_pheno.txt --og /TMP/nfbc66.gen --log nfbc66_1.log
gtool -S --g /GENOS/nfbc86.gen.gz --exclusion indels.txt --s nfbc86_pheno.txt --og /TMP/nfbc86.gen --log nfbc86_1.log
gtool -S --g /GENOS/KF4.gen.gz --exclusion indels.txt --s KF4_pheno.txt --og /TMP/KF4.gen --log KF4_1.log
gtool -S --g /GENOS/IA317.gen.gz --exclusion indels.txt --s IA317_pheno.txt --og /TMP/IA317.gen --log IA317_1.log
gtool -S --g /GENOS/IA610.gen.gz --exclusion indels.txt --s IA610_pheno.txt --og /TMP/IA610.gen --log IA610_1.log
gtool -S --g /GENOS/OmniX.gen.gz --exclusion indels.txt --s OmniX_pheno.txt --og /TMP/OmniX.gen --log OmniX_1.log
gtool -S --g /GENOS/OmniExome.gen.gz --exclusion indels.txt --s OmniExome_pheno.txt --og /TMP/OmniExome.gen --log OmniExome_1.log


gtool -G --g /TMP/nfbc66.gen.gz --s nfbc66_pheno.txt --ped /TMP/nfbc66.ped --map /TMP/nfbc66.map --snp --threshold 0.98 --log nfbc66_2.log
gtool -G --g /TMP/nfbc86.gen.gz --s nfbc86_pheno.txt --ped /TMP/nfbc86.ped --map /TMP/nfbc86.map --snp --threshold 0.98 --log nfbc86_2.log
gtool -G --g /TMP/KF4.gen.gz --s KF4_pheno.txt --ped /TMP/KF4.ped --map /TMP/KF4.map --snp --threshold 0.98 --log KF4_2.log
gtool -G --g /TMP/IA317.gen.gz --s IA317_pheno.txt --ped /TMP/IA317.ped --map /TMP/IA317.map --snp --threshold 0.98 --log IA317_2.log
gtool -G --g /TMP/IA610.gen.gz --s IA610_pheno.txt --ped /TMP/IA610.ped --map /TMP/IA610.map --snp --threshold 0.98 --log IA610_2.log
gtool -G --g /TMP/OmniX.gen.gz --s OmniX_pheno.txt --ped /TMP/OmniX.ped --map /TMP/OmniX.map --snp --threshold 0.98 --log OmniX_2.log
gtool -G --g /TMP/OmniExome.gen.gz --s OmniExome_pheno.txt --ped /TMP/OmniExome.ped --map /TMP/OmniExome.map --snp --threshold 0.98 --log OmniExome_2.log


module load plink/1.90
plink --file /TMP/nfbc66 --extract cond.snps --r2 --ld-snp-list cond.snps --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --noweb
plink --file /TMP/nfbc86 --extract cond.snps --r2 --ld-snp-list cond.snps --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --noweb
plink --file /TMP/KF4 --extract cond.snps --r2 --ld-snp-list cond.snps --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --noweb
plink --file /TMP/IA317 --extract cond.snps --r2 --ld-snp-list cond.snps --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --noweb
plink --file /TMP/IA610 --extract cond.snps --r2 --ld-snp-list cond.snps --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --noweb
plink --file /TMP/OmniX --extract cond.snps --r2 --ld-snp-list cond.snps --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --noweb
plink --file /TMP/OmniExome --extract cond.snps --r2 --ld-snp-list cond.snps --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --noweb
```

read into R
[cohort].ld: LD reports for each cohort (columns:CHR_A, BP_A, SNP_A, CHR_B, BP_B, SNP_B, R2)

```{r}

#shown for OmniX; repeat for all cohorts
ld=read.table('OmniX.ld', header=T, as.is=T)
ld=ld[ld$SNP_A!=ld$SNP_B,]
snpa=ifelse(ld$SNP_A>ld$SNP_B, ld$SNP_A, ld$SNP_B)
snpb=ifelse(ld$SNP_A>ld$SNP_B, ld$SNP_B, ld$SNP_A)
ld.mod=cbind(snpa, snpb, ld$R2)
OmniX.ld=unique(ld.mod)
save(OmniX.ld, file='OmniX_ld.RData')

load('nfbc66_ld.RData')
load('nfbc86_ld.RData')
load('KF4_ld.RData')
load('IA317_ld.RData')
load('IA610_ld.RData')
load('OmniX_ld.RData')
load('OmniExome_ld.RData')
pairs=unique(c(paste(nfbc66.ld[,1], nfbc66.ld[,2], sep='_'), paste(nfbc86.ld[,1], nfbc86.ld[,2], sep='_'), paste(KF4.ld[,1], KF4.ld[,2], sep='_'), paste(ia317.ld[,1], ia317.ld[,2], sep='_'),paste(ia610.ld[,1], ia610.ld[,2], sep='_'), paste(OmniX.ld[,1], OmniX.ld[,2], sep='_'), paste(OmniExome.ld[,1], OmniExome.ld[,2], sep='_')))
save(pairs, file='pairs.RData')

# combine
cohorts=c('IA317', 'IA610', 'OmniX', 'OmniExome', 'nfbc66', 'nfbc86', 'KF4') 
columns.gen=c(987, 2748, 1788, 4062, 2196, 1542,4968)
samples=columns.gen/3
names(samples)=cohorts
load('IA317_ld.RData')
ia317.ld[,3]=as.numeric(as.character(ia317.ld[,3])) * samples['IA317']
rownames(ia317.ld)=paste(ia317.ld[,1], ia317.ld[,2])
load('KF4_ld.RData')
KF4=as.numeric(as.character(KF4.ld[,3])) * samples['KF4']
names(KF4)=paste(KF4.ld[,1], KF4.ld[,2])
m1=merge(ia317.ld, KF4, by.x='row.names', by.y='row.names')
load('IA610_ld.RData')
ia610=as.numeric(as.character(ia610.ld[,3])) * samples['IA610']
names(ia610)=paste(ia610.ld[,1], ia610.ld[,2])
m2=merge(m1, ia610, by.x=1, by.y='row.names')
load('OmniX_ld.RData')
OmniX=as.numeric(as.character(OmniX.ld[,3])) * samples['OmniX']
names(OmniX)=paste(OmniX.ld[,1], OmniX.ld[,2])
m3=merge(m2, OmniX, by.x=1, by.y='row.names')
colnames(m3)[4:7]=c('ia317', 'kf4', 'IA610', 'OmniX')
save(m3, file='/project/lolipop_b/METHQTL2/STUFF/PRUNE/COND2/MERGE/m3_tmp.RData')
load('OmniExome_ld.RData')
oE=as.numeric(as.character(oE.ld[,3])) * samples['OmniExome']
oE=data.frame(cbind(ID=paste(oE.ld[,1], oE.ld[,2]), oE))
colnames(m3)[1]='ID'
require(plyr)
m4=join(m3, oE, by='ID', type='left')
load('nfbc66_ld.RData')
nfbc66=as.numeric(as.character(nfbc66.ld[,3])) * samples['nfbc66']
nfbc66=data.frame(cbind(ID=paste(nfbc66.ld[,1], nfbc66.ld[,2]), nfbc66))
m5=join(m4, nfbc66, by='ID', type='left')
load('nfbc86_ld.RData')
nfbc86=as.numeric(as.character(nfbc86.ld[,3])) * samples['nfbc86']
nfbc86=data.frame(cbind(ID=paste(nfbc86.ld[,1], nfbc86.ld[,2]), nfbc86))
m6=join(m5, nfbc86, by='ID', type='left')
m6[,4]=as.numeric(as.character(m6[,4]))/samples['IA317']
m6[,5]=as.numeric(as.character(m6[,5]))/samples['KF4']
m6[,6]=as.numeric(as.character(m6[,6]))/samples['IA610']
m6[,7]=as.numeric(as.character(m6[,7]))/samples['OmniX']
m6[,8]=as.numeric(as.character(m6[,8]))/samples['OmniExome']
m6[,9]=as.numeric(as.character(m6[,9]))/samples['nfbc66']
m6[,10]=as.numeric(as.character(m6[,10]))/samples['nfbc86']
weighted=rowMeans(m6[,4:10], na.rm=T)
ld=cbind(m6[,2:3],weighted)
save(ld, file='ld_combined.RData')
```

## R2 merging and pruning of SNPs

```{r}
# Merge/Prune SNP; can be parallelized by chr (tid=chr no.)  
load('ld_combined.RData')
betacor=data.frame(ld)
betacor$snpa=as.character(betacor$snpa)
betacor$snpb=as.character(betacor$snpb)
colnames(betacor)[3]='r2'
betacor$r2=as.numeric(as.character(betacor$r2))
r2.thres=0.02
cis.window=1000000
load('final.RData')
sigtab=rsigtab[rsigtab$snp.chr==tid,]
betacor=betacor[betacor$snpa %in% sigtab$snp & betacor$snpb %in% sigtab$snp,]
sigtab$z=sigtab$beta.comb/sigtab$se.comb
sigtab=sigtab[with(sigtab, order(abs(z), decreasing=T)),]
sigtab$snp.pos=as.numeric(as.character(sigtab$snp.pos))
sigtab=unique(sigtab[,c('snp','snp.chr', 'snp.pos')])
colnames(sigtab)=c('snp', 'chr', 'pos')
sigtab$snp=as.character(sigtab$snp)
sigtab$pos=as.numeric(as.character(sigtab$pos))
sigtab$chr=as.numeric(as.character(sigtab$chr))
sigtab$locus=rep(NA, nrow(sigtab))
for(i in 1:nrow(sigtab)){
        if(is.na(sigtab$locus[i])){
                markers=as.character(sigtab[sigtab$chr==sigtab$chr[i] & abs(sigtab$pos-sigtab$pos[i])<cis.window & is.na(sigtab$locus),'snp'])
                if(length(markers)>0){
                        r2.a=betacor[betacor$snpa==sigtab$snp[i] & betacor$snpb %in% markers,c(2,3)]
                        colnames(r2.a)=c('snp', 'r2')
                        r2.b=betacor[betacor$snpb==sigtab$snp[i] & betacor$snpa %in% markers,c(1,3)]
                        colnames(r2.b)=c('snp', 'r2')
                        r2=rbind(r2.a,r2.b)
                        r2.markers=as.character(r2[r2[,2]>r2.thres,1])
                        r2.markers=c(r2.markers,sigtab$snp[i])
                }else{
                        r2.markers=sigtab$snp[i]
                }
                sigtab[sigtab$snp %in% r2.markers & is.na(sigtab$locus),'locus']=i
        }
        print(i)
}

save(sigtab, file=paste(tid, '.RData', sep=''))
```

## Merge and pruning of CpGs

MARIE: What are these files
anno.RData: consists of "annometh" and "annosnp" (rsid/cpg, chr, pos)

```{r}
## can be parallelized by chr (tid=chr no.)
cis.window=1000000

#CpGs
load('final.RData')
cond.cpgs=as.character(unique(rsigtab[rsigtab$cpg.chr==tid,'cpg']))
load('anno.RData')
annometh.cond=annometh[annometh$rsid %in% cond.cpgs,]
rownames(annometh.cond)=as.character(annometh.cond$rsid)
annometh.cond=annometh.cond[as.character(cond.cpgs),]
cohorts=c('lolipop', 'nfbc66', 'nfbc86', 'KF4')
for(cohort in cohorts){
	betacor.all=rep(NA, 3)
	load(paste('betaQN_', cohort, '.RData', sep=''))
	nsamples=ncol(beta)
	for(i in 1:length(cond.cpgs)){
	        cpgs=intersect(cond.cpgs, annometh[ abs(annometh$pos-annometh.cond$pos[i])<cis.window,'rsid'])
	        betacor=matrix(nrow=length(cpgs), ncol=3)
	        for(c in 1:length(cpgs)){
	                betacor[c,1]=as.character(cond.cpgs[i])
	                betacor[c,2]=as.character(cpgs[c])
	                betacor[c,3]=summary(lm(beta[as.character(cond.cpgs[i]),]~beta[as.character(cpgs[c]),]))$r.squared
	        }
	        betacor.all=rbind(betacor.all, betacor)
	}
	colnames(betacor.all)=c('cpg1', 'cpg2', 'r2')
	betacor.all=data.frame(betacor.all)
	betacor.all$cpg1=as.character(betacor.all$cpg1)
	betacor.all$cpg2=as.character(betacor.all$cpg2)
	betacor.all[,3]=as.numeric(as.character(betacor.all[,3]))
	assign(paste(cohort, '.betacor', sep=''), betacor.all)
	assign(paste(cohort, '.nsamples', sep=''), nsamples)
}
save(lolipop.betacor,nfbc66.betacor, nfbc86.betacor, KF4.betacor, lolipop.nsamples,nfbc66.nsamples,nfbc86.nsamples, KF4.nsamples, file=paste('tmp_chr', tid, '.RData', sep=''))

tmp=((lolipop.betacor[,3]*lolipop.nsamples)+(nfbc66.betacor[,3]*nfbc66.nsamples)+(nfbc86.betacor[,3]*nfbc86.nsamples)+(KF4.betacor[,3]*KF4.nsamples))/(lolipop.nsamples+nfbc66.nsamples+nfbc86.nsamples+KF4.nsamples)
betacor.all=cbind(lolipop.betacor[,1:2], r2=tmp)
save(betacor.all, file=paste('chr', tid, '.RData', sep=''))

c=1
load(paste('chr', c, '.RData', sep=''))
betacor=betacor.all[2:nrow(betacor.all),]
for(c in 2:22){
        load(paste('chr', c, '.RData', sep=''))
        betacor=rbind(betacor, betacor.all[2:nrow(betacor.all),])
}
colnames(betacor)=c('cpg1', 'cpg2', 'r2')
betacor=data.frame(betacor)
betacor$cpg1=as.character(betacor$cpg1)
betacor$cpg2=as.character(betacor$cpg2)
betacor=betacor[betacor$cpg1!=betacor$cpg2,]
cpg1=ifelse(betacor$cpg1>betacor$cpg2, betacor$cpg1, betacor$cpg2)
cpg2=ifelse(betacor$cpg1>betacor$cpg2, betacor$cpg2, betacor$cpg1)
betacor.mod=cbind(cpg1, cpg2, r2=round(as.numeric(as.character(betacor$r2)), digits=5))
betacor=unique(betacor.mod)
save(betacor, file='cpg_r2.RData')

load('cpg_r2.RData')
betacor=data.frame(betacor)
betacor$cpg1=as.character(betacor$cpg1)
betacor$cpg2=as.character(betacor$cpg2)
betacor$r2=as.numeric(as.character(betacor$r2))
r2.thres=0.2
cis.window=1000000
load('final.RData')
sigtab=rsigtab
sigtab$z=as.numeric(as.character(sigtab$z))
sigtab$cpg.pos=as.numeric(as.character(sigtab$cpg.pos))
sigtab=sigtab[with(sigtab, order(abs(z), decreasing=T)),]
sigtab=unique(sigtab[,c('cpg','cpg.chr', 'cpg.pos')])
sigtab$cpg=as.character(sigtab$cpg)
sigtab$cpg.pos=as.numeric(as.character(sigtab$cpg.pos))
sigtab$cpg.chr=as.numeric(as.character(sigtab$cpg.chr))
sigtab$cpg.locus=rep(NA, nrow(sigtab))
for(i in 1:nrow(sigtab)){
        if(is.na(sigtab$cpg.locus[i])){
                markers=as.character(sigtab[sigtab$cpg.chr==sigtab$cpg.chr[i] & abs(sigtab$cpg.pos-sigtab$cpg.pos[i])<cis.window & is.na(sigtab$cpg.locus),'cpg'])
                if(length(markers)>0){
                        r2.a=betacor[betacor$cpg1==sigtab$cpg[i] & betacor$cpg2 %in% markers,c(2,3)]
                        colnames(r2.a)=c('cpg', 'r2')
                        r2.b=betacor[betacor$cpg2==sigtab$cpg[i] & betacor$cpg1 %in% markers,c(1,3)]
                        colnames(r2.b)=c('cpg', 'r2')
                        r2=rbind(r2.a,r2.b)
                        r2.markers=as.character(r2[r2[,2]>r2.thres,1])
                        r2.markers=c(r2.markers,sigtab$cpg[i])
                }else{
                        r2.markers=sigtab$cpg[i]
                }
        }
}

save(sigtab, file='cpg_loci_genomicCond_r2.RData')
```

