#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

id = arguments[1]

mt = read.delim(paste0('mt-calls/',id,'.summary.txt'))

# Subset to the cytb+dloop region of Zinner et al
cytb = subset(mt,Pos >=14177 & Pos < 15910)

key = c(
	'A'    = 'A',
	'C'    = 'C',
	'G'    = 'G',
	'T'    = 'T',
	'AG'   = 'r',
	'CT'   = 'y',
	'CG'   = 's',
	'AT'   = 'w',
	'GT'   = 'k',
	'AC'   = 'm',
	'CGT'  = 'b',
	'AGT'  = 'd',
	'ACT'  = 'h',
	'ACG'  = 'v',
	'ACGT' = 'n'
)

if (identical(cytb$Pos,14177:15909)) {
	out = apply(
		matrix(as.integer(gsub('[ACGT]:','',as.matrix(cytb[,c('FilteredNumA','FilteredNumC','FilteredNumG','FilteredNumT')]))),ncol=4),1,function(x) {
			this.call = paste(c('A','C','G','T')[which(x == x[which.max(x)])],collapse='')
			key[this.call]
		}
	)
	# cat(paste0('>',id,'\n'),paste(apply(matrix(as.integer(gsub('[ACGT]:','',as.matrix(cytb[,c('FilteredNumA','FilteredNumC','FilteredNumG','FilteredNumT')]))),ncol=4),1,function(x) if ((x[which.max(x)] - max(x[setdiff(1:4,which.max(x))])) <= 0) 'n' else c('A','C','G','T')[which.max(x)]),collapse=''),'\n',sep='')
} else {
	out = character(length(14177:15909))
	names(out) = 14177:15909
	out[as.character(setdiff(14177:15909,cytb$Pos))] = 'n'
	out[as.character(intersect(14177:15909,cytb$Pos))] = apply(
		matrix(as.integer(gsub('[ACGT]:','',as.matrix(cytb[,c('FilteredNumA','FilteredNumC','FilteredNumG','FilteredNumT')]))),ncol=4),1,function(x) {
			this.call = paste(c('A','C','G','T')[which(x == x[which.max(x)])],collapse='')
			key[this.call]
		}
	)
}
cat(paste0('>',id,'\n'),paste(out,collapse=''),'\n',sep='')
