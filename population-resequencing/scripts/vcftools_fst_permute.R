#!/usr/bin/env Rscript

source('scripts/_include_options.sh')

arguments = commandArgs(trailing=TRUE)

this.step = arguments[1]
n.cores = if (length(arguments) > 1) as.integer(arguments[2]) else parallel::detectCores(logical=FALSE)

library(data.table)

# Read in windowed fst results
window.fst = read.delim(file.path('results/fst',paste0(dataset,'.',genome,'.bootstrap.all.step',this.step,'.windowed.weir.fst')))

# Read in site fst results
site.fst = read.delim(file.path('results/fst',paste0(dataset,'.',genome,'.bootstrap.all.step',this.step,'_persite.weir.fst')))

# Calculate all the unique numbers of variants (the null distribution for all windows with the same N_VARIANTS is identical so no need to waste effort computing them multiple times)
n.variants = sort(unique(window.fst$N_VARIANTS))

library(parallel)

# Read in background genome-wide distribution of site Fst
all.fst = na.omit(site.fst$WEIR_AND_COCKERHAM_FST)

# For each unique number of variants within window, calculate the mean Fst from N sampled variants. Do this 1e6 times
null.distribution = mclapply(n.variants,function(i) {
	replicate(1e6,mean(sample(all.fst,i)))
},mc.cores=n.cores)

# Name them by number of variants for easy retrieval
names(null.distribution) = n.variants

# Iterate over windows and use the null distribution for the given number of variants to calculate empirical P values
# Need to use mean fst instead of weighted fst because the numerators and denominators are not available
window.fst$pval = unlist(mclapply(1:nrow(window.fst),function(i) {
	with(window.fst[i,],mean(null.distribution[[as.character(N_VARIANTS)]] >= MEAN_FST))
},mc.cores=n.cores))

# Calculate q value
window.fst$qval = p.adjust(window.fst$pval,'fdr')

dir.create('checkpoints',showWarnings=FALSE)

saveRDS(window.fst,file=file.path('checkpoints',paste0(dataset,'.',genome,'.bootstrap.all.step',this.step,'.fst_windows.rds')))

# Merge consecutive bases

sig.windows = subset(window.fst,qval < 0.001)

# Split into chromosomes so we can calculate consecutive windows based on coordinates
sig.windows.split = split(sig.windows,sig.windows$CHROM)

sig.window.merged = do.call(rbind,mclapply(sig.windows.split,function(x) {
	# Calculate gap between start coordinates and subtract 1000. "Redundant" (overlapping and consecutive windows) will become <= 0
	x$gap = c(Inf,diff(x$BIN_START)-1000)
	# When gap is > 0, windows are true starts of a new block
	which.start = which(x$gap > 0)
	# Calculate indices from the start position to the position of the next block
	these.indices = cbind(which.start,diff(c(which.start,nrow(x)+1)))
	# Extract the ranges using the above coordinates
	these.ranges = apply(these.indices,1,function(y) y[1]:(y[1]+y[2]-1))
	# For each of the ranges, subset those windows and merge them. Calculate summary stats
	do.call(rbind,lapply(these.ranges,function(i) {
		these.windows = x[i,]
		out.window = with(these.windows,data.frame(
			CHROM=unique(CHROM),
			BIN_START=min(BIN_START),
			BIN_END=max(BIN_END),
			LENGTH=max(BIN_END)-min(BIN_START)+1,
			N_WINDOWS=length(CHROM),
			N_VARIANTS=sum(N_VARIANTS),
			WEIGHTED_FST_MEAN=mean(WEIGHTED_FST),
			MEAN_FST=sum(MEAN_FST*N_VARIANTS)/sum(N_VARIANTS),
			MEAN_PVAL=mean(pval),
			MIN_QVAL=min(qval),
			MAX_QVAL=max(qval)
		))
		out.window
	}))
},mc.cores=n.cores))

rownames(sig.window.merged) = NULL
sum(sig.window.merged$LENGTH)

saveRDS(sig.window.merged,file=file.path('checkpoints',paste0(dataset,'.',genome,'.bootstrap.all.step',this.step,'.fst_sig_windows.rds')))
