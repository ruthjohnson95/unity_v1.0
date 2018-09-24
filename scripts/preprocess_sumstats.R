#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(logging)
library(optparse)
library(R.utils)


# setup global logging 
basicConfig()

filter_MAF <- function(gwas, freq_head, snp_head, THRESH){

	if (freq_head == "NA"){
		loginfo("User did not provide minor allele frequency column header...inferring header")
		
		# search for column name for minor allele frequency (MAP, FREQ, maf, freq)
		
		if("MAF" %in% colnames(gwas))
		{
		  loginfo("Found 'MAF' as minor allele frequency column")
		  loginfo("Renaming column as 'FREQ'")
		  colnames(gwas)[colnames(gwas) == 'MAF'] <- 'FREQ'
		}
		else if("maf" %in% colnames(gwas))
		{
			loginfo("Found 'maf' as minor allele frequency column")
			loginfo("Renaming column as 'FREQ'")
			colnames(gwas)[colnames(gwas) == 'maf'] <- 'FREQ'
		}
		else if("FREQ" %in% colnames(gwas))
		{
			loginfo("Found 'FREQ' as minor allele frequency column")
		}
		else if("freq" %in% colnames(gwas))
		{
			loginfo("Found 'freq' as minor allele frequency column")
			loginfo("Renaming column as 'FREQ'")
			colnames(gwas)[colnames(gwas) == 'freq'] <- 'FREQ'
		}
		else{
			# gwas does not contain minor allele frequency
			loginfo("Did not find minor allele freq column...going to use maf from 1000G and match via rsid")

			# list of SNP freq from 1000G
			freq_table="misc/all.chr.frq"
			maf_df <- fread(freq_table, header=T, showProgress=FALSE)
			maf_df<-data.frame(maf_df$SNP, maf_df$MAF)
			colnames(maf_df)<-c("SNP", "FREQ")


			# changing SNP header to 'SNP' if not already 
			if('SNP' %in% colnames(gwas))
			{
				loginfo("Using 'SNP' as rsid column header")
			}
			else{
				if(is.null(snp_head)){

					# try alternate headers 

					if("RSID" %in% colnames(gwas)){
						loginfo("Found 'RSID' as rsid column")
						loginfo("Renaming rsid column as 'SNP'")
						colnames(gwas)[colnames(gwas) == 'RSID'] <- 'SNP'
					}
					else if("rsid" %in% colnames(gwas)){
						loginfo("Found 'rsid' as rsid column")
						loginfo("Renaming rsid column as 'SNP'")
						colnames(gwas)[colnames(gwas) == 'rsid'] <- 'SNP'
					}
					else if('snp' %in% colnames(gwas)){
						loginfo("Found 'snp' as rsid column")
						loginfo("Renaming rsid column as 'SNP'")
						colnames(gwas)[colnames(gwas) == 'snp'] <- 'SNP'
					}
					else{ # cannot infer SNP column 
						stop("Cannot infer rsid column!", call.=FALSE)
					}
				}
				else{
					# use user defined SNP header 
					loginfo("User provided %s as rsid column header", snp_head)

					if(snp_head != "SNP"){
						loginfo("Renaming rsid column as 'SNP'")
						colnames(gwas)[colnames(gwas) == snp_head] <- 'SNP'
					}
				}
			}

			gwas <- inner_join(gwas, maf_df, by="SNP")

		}

	}
	else{ # user proivded header for minor allele frequency column 
		loginfo("User provided %s as minor allele freq column name", freq_head)

		if(freq_head != 'FREQ'){
			loginfo("Renaming column as 'FREQ'")
			colnames(gwas)[colnames(gwas) == freq_head] <- 'FREQ'
		}
	}

	# Filter by MAF <= THRESH 

	# replace NA minor allele frequencies
	gwas$FREQ[is.na(gwas$FREQ)]	<- 0

	maf_inds <- which(gwas$FREQ	>= THRESH)

	gwas <- gwas[maf_inds, ]

	return(gwas)
}


zscores_to_betas <- function(gwas, z_head, N_head){

	# if zhead is null, try to infer header 
	if(is.null(z_head)){
		# try to infer
		if("Z" %in% colnames(gwas)){
			loginfo("Found 'Z' as zscore column")
		} 
		else if("Zscore" %in% colnames(gwas)){

			loginfo("Found 'Zscore' as zscore column")
			loginfo("Renaming zscore column as 'Z")
			colnames(gwas)[colnames(gwas) == 'Zscore'] <- 'Z'
		}
		else{
			stop("Could not infer column header for zscores!", call.=FALSE)
		}

	}
	else{
		if(z_head != 'Z'){
			# rename header 
			loginfo("Renaming zscore column as 'Z'")
			colnames(gwas)[colnames(gwas) == z_head] <- 'Z'
		}
	}

	# standardized effect sizes and standard errors  
	# SE = 1/sqrt(N)
	# Z = beta/SE(beta)
	# beta = Z*SE(beta)

	SE<-1/sqrt(gwas$N)
	gwas$BETA_STD <- gwas$Z*SE 
	gwas$SE_STD

	return(gwas)
}


overlap_sumstats <- function(gwas_1, gwas_2){

	# overlap SNPs by rsid 
	overlap_snps_1 <- which(gwas_1$SNP %in% gwas_2$SNP)
	overlap_snps_2 <- which(gwas_2$SNP %in% gwas_1$SNP)

	# filter by overlapping snps
	gwas_1_overlap <- gwas_1[overlap_snps_1, ]
	gwas_2_overlap <- gwas_2[overlap_snps_2, ]

	return(list(gwas_1_overlap, gwas_2_overlap))
}


main<-function(){
	# setup command line arguments 
	option_list = list(
		make_option(c("--gwas_file_1"), type="character", default=NULL, help="gwas file name for Trait 1", metavar="character"),
		make_option(c("--gwas_file_2"), type="character", default=NULL, help="gwas file name for Trait 2", metavar="character"),
		make_option(c("--snp_head"), type="character", default=NULL, help="header for rsid column"),
		make_option(c("--freq_head"), type="character", default=NULL, help="header for minor allele freq column", metavar="character"),
		make_option(c("--z_head"), type="character", default=NULL, help="header for zscore column", metavar="character"),
		make_option(c("--N_head"), type="character", default=NULL, help="header for sample size column", metavar="character"),
		make_option(c("--trait_1"), type="character", default=NULL, help="Name for Trait 1", metavar="character"),
		make_option(c("--trait_2"), type="character", default=NULL, help="Name for Trait 2", metavar="character"),
		make_option(c("--maf_thresh"), type="double", default=0.05, help="minor allele freq filetering threshold", metavar="double"),
		make_option(c("--outdir"), type="character", help="output directory for processed sumstats")
	); 

	opt_parser = OptionParser(option_list=option_list);
	opt = parse_args(opt_parser);

	loginfo("----- Pre-processing summary statistics -----")

	if (is.null(opt$gwas_file_1)){
	  stop("GWAS file for Trait 1 has not been provided!", call.=FALSE)
	}
	else if(is.null(opt$gwas_file_2)){
	  stop("GWAS file for Trait 2 has not been provided!", call.=FALSE)
	}
	else{
		gwas_file_1<-opt$gwas_file_1 
		loginfo('Reading in sumstats file for Trait 1: %s', gwas_file_1)
		gwas_file_2<-opt$gwas_file_2
		loginfo('Reading in sumstats file for Trait 2: %s', gwas_file_2)
	}

	# outdir assumes given file structure 
	out_dir="../results/processed_gwas"

	# check for .gz format 
	if(grepl('.gz', gwas_file_1, fixed=TRUE)){
		gwas_1 <- fread(input = sprintf("zcat %s", gwas_file_1), header = TRUE, showProgress = FALSE,fill=T)
	}
	else{
		gwas_1 <- fread(gwas_file_1, header=TRUE, showProgress=FALSE, fill=T)
	}

	if(grepl('.gz', gwas_file_2, fixed=TRUE)){
		gwas_2 <- fread(input = sprintf("zcat %s", gwas_file_2), header = TRUE, showProgress = FALSE,fill=T)
	}
	else{
		gwas_2 <- fread(gwas_file_2, header=T, showProgress=F, fill=T)
	}
	
	freq_head <- opt$freq_head
	snp_head <- opt$snp_head 
	z_head <- opt$z_head 
	N_head <- opt$N_head 
	trait_1 <- opt$trait_1 
	trait_2 <- opt$trait_2
	maf_thresh <- opt$maf_thresh
	outdir <- opt$outdir

	# add MAF 
	loginfo('STEP 0: Filtering by minor allele frequencies')
	loginfo('Using maf threshold: %.2f', maf_thresh)
	gwas_1 <- filter_MAF(gwas_1, freq_head, snp_head, maf_thresh)
	gwas_2 <- filter_MAF(gwas_2, freq_head, snp_head, maf_thresh)

	# covert zscores to standardized betas 
	loginfo('STEP 1: Converting zscores to standardized betas')
	gwas_1 <- zscores_to_betas(gwas_1, z_head, N_head)
	gwas_2 <- zscores_to_betas(gwas_2, z_head, N_head) 

	# overlap sumstats 
	loginfo('STEP 2: Overlapping sumstats by rsid')
	gwas_overlap_list <- overlap_sumstats(gwas_1, gwas_2) 
	gwas_1_overlap <- gwas_overlap_list[1]
	gwas_2_overlap <- gwas_overlap_list[2]

	# print to files 
	trait_pair_name_1 <- paste(trait_1, trait_2, sep='_')
	trait_pair_filename_1 <- paste(trait_pair_name_1, 'processed', sep='.')
	trait_pair_name_2 <- paste(trait_2, trait_1, sep='_')
	trait_pair_filename_2 <- paste(trait_pair_name_2, 'processed', sep='.')

	outname_1 <- paste(outdir, trait_pair_filename_1, sep='/')
	outname_2 <- paste(outdir, trait_pair_filename_2, sep='/')

	write.table(gwas_1_overlap, file=outname_1, quote=F, sep=' ', row.names=F)
	write.table(gwas_2_overlap, file=outname_2, quote=F, sep=' ', row.names=F)

}


if(getOption("run.main", default=TRUE)) {
   main()
}





