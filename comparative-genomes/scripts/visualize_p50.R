#!/usr/bin/env Rscript

p50 = read.delim('data/p50.txt')

library(ggplot2)
library(RColorBrewer)

# Set taxon order
p50$taxon = factor(p50$taxon,levels=c('gelada','baboon','human'))
p50$set = factor(p50$set,levels=c('+ cofactors','- cofactors'))

# summary(lm(p50~taxon+set,data=p50))

# Split raw data by taxon and set (+/- cofactors)
p50.split = split(p50,list(p50$taxon,p50$set))

p50.summary = do.call(rbind,lapply(p50.split,function(x) {

	x$logp50 = log10(x$p50)
	out = unique(x[c('taxon','set')])
	out$p50.mean = 10^unlist(predict(lm(logp50~ph50,data=x),data.frame(ph50=7.4)))
	# Standard error of fitted value
	out$p50.sef = 10^unlist(predict(lm(logp50~ph50,data=x),data.frame(ph50=7.4),se.fit=TRUE)$se.fit)
	# Standard error of estimate
	out$p50.see = 10^summary(lm(logp50~ph50,data=x))$sigma
	# Confidence intervals
	out$conf.ub = 10^predict(lm(logp50~ph50,data=x),data.frame(ph50=7.4),se.fit=TRUE,interval='confidence')$fit[,'upr']
	out$conf.lb = 10^predict(lm(logp50~ph50,data=x),data.frame(ph50=7.4),se.fit=TRUE,interval='confidence')$fit[,'lwr']
	# For plotting, mean +/- the standard errors of the predicted values
	out$ub = with(out,p50.mean + p50.see)
	out$lb = with(out,p50.mean - p50.see)
	# sample size (number of measurements
	out$n = nrow(x)
	out
}))

# Function to test that the logP50 of group 1 is less than the logP50 of group 2
# Statistical basis explained in 10.2307/2348167
compare.fitted.values = function(data1,data2,x=7.4) {
	# data1 and data2 are a longform data frame with repeated measurements of P50 (p50) at different PH (column ph50) values
	# data1 = p50.split[['gelada.- cofactors']]
	# data2 = p50.split[['baboon.- cofactors']]
	# x = 7.4

	data1$logp50 = log10(data1$p50)
	lm1 = lm(logp50~ph50,data=data1)
	data2$logp50 = log10(data2$p50)
	lm2 = lm(logp50~ph50,data=data2)

#	# Making sure I got the inputs to the equations right: eq 1 (with real logp50 values)
#	# Start with the case of [x,y] = [0,y intercept]
#	eq.left = with(data1,logp50-coef(lm1)['(Intercept)'])
#	eq.right = with(data1,coef(lm1)['ph50'] * (ph50 - 0))
#	# Now do the case of [x,y] = [7.4,predicted values at pH 7.4] # or replace 7.4 with the variable x (argument of this function)
#	eq.left = with(data1,logp50-predict(lm1,data.frame(ph50=x)))
#	eq.right = with(data1,coef(lm1)['ph50'] * (ph50 - x))
#
#	# Proof that I got the equations right: eq 1 (with predicted logp50 values)
#	# Start with the case of [x,y] = [0,y intercept]
#	eq.left = with(data1,predict(lm1,data.frame(ph50=c(ph50)))-coef(lm1)['(Intercept)'])
#	eq.right = with(data1,coef(lm1)['ph50'] * (ph50 - 0))
#	# Now do the case of [x,y] = [7.4,predicted values at pH 7.4] # or replace 7.4 with the variable x (argument of this function)
#	eq.left = with(data1,predict(lm1,data.frame(ph50=c(ph50)))-predict(lm1,data.frame(ph50=x)))
#	eq.right = with(data1,coef(lm1)['ph50'] * (ph50 - x))
#
#	# eq.left and eq.right match up nearly exactly (differences are due to rounding) (success!!)
#
#	# Now prove equation 2 (looks good)
#	d.test = (coef(lm1)['(Intercept)']-coef(lm2)['(Intercept)']) + coef(lm1)['ph50']*(x-0) - coef(lm2)['ph50']*(x-0)
#	d.eval = predict(lm1,data.frame(ph50=x)) - predict(lm2,data.frame(ph50=x))

#	# Equation 2 produces the exact same result as the difference in predicted values through predict()
#	identical(d.test,d.eval) # TRUE

	# n1 and n2 are the sample sizes (number of measurements in the regressions) of data1 and data2
	n1 = nrow(data1)
	n2 = nrow(data2)

	# For the t-test, keep all values in log scale

	# s1 and s2 are the standard errors of the estimates (SEE) from the regressions for data1 and data2
	s1 = summary(lm1)$sigma
	s2 = summary(lm2)$sigma

	# Estimate the combined variance (S2): see equation 4 from 10.2307/2348167
	S2 = ((n1-2)*s1^2 + (n2-2)*s2^2) / (n1 + n2 - 4)

	# Estimate the standard error of the difference (SED): see equation 5 from 10.2307/2348167
	SED = sqrt(S2) * sqrt(1/n1 + 1/n2 + (x - 0)^2/sum((data1$ph50 - 0)^2) + (x - 0)^2/sum((data2$ph50 - 0)^2) )

	# Calculate the difference in predicted values (D)
	D = predict(lm1,data.frame(ph50=x)) - predict(lm2,data.frame(ph50=x))

	# Calculate the t statistic (T) and degrees of freedom (DF): see equation 6 from 10.2307/2348167
	T = D/SED
	DF = n1 + n2 - 4

	# Because we expect that log P50 for data1 is less than log P50 for data2, reverse the sign of T and run a one-sided test
	pval = pt(-T,DF,lower.tail=FALSE)
	as.numeric(pval)
}

compare.fitted.values(p50.split[['gelada.+ cofactors']],p50.split[['baboon.+ cofactors']],x=7.4)
compare.fitted.values(p50.split[['gelada.- cofactors']],p50.split[['baboon.- cofactors']],x=7.4)

compare.fitted.values(p50.split[['gelada.+ cofactors']],p50.split[['human.+ cofactors']],x=7.4)
compare.fitted.values(p50.split[['gelada.- cofactors']],p50.split[['human.- cofactors']],x=7.4)


p = ggplot(p50.summary,aes(set,p50.mean,fill=taxon)) +
	geom_bar(stat='identity',position='dodge') +
	geom_errorbar(aes(ymin=lb,ymax=ub),position=position_dodge2(width = 0.25, padding = 0.75)) +
	theme_classic(base_size=24) +
	scale_fill_manual(values=c('#377eb8','#4daf4a','#984ea3')) +
	scale_y_continuous(trans='log10',breaks=seq(0,15,5)) +
	ylab(expression(log[10](italic(P['50'])))) +
	theme(axis.title.x=element_blank(),legend.title = element_blank())
ggsave(p,file='figures/gelada_p50_log.pdf',useDingbats=FALSE)

p = ggplot(p50.summary,aes(set,p50.mean,fill=taxon)) +
	geom_bar(stat='identity',position='dodge') +
	geom_errorbar(aes(ymin=lb,ymax=ub),position=position_dodge2(width = 0.25, padding = 0.75)) +
	theme_classic(base_size=24) +
	scale_fill_manual(values=c('#377eb8','#4daf4a','#984ea3')) +
	scale_y_continuous(breaks=seq(0,15,5)) +
	ylab(expression(italic(P['50']))) +
	theme(axis.title.x=element_blank(),legend.title = element_blank())
ggsave(p,file='figures/gelada_p50_nolog.pdf',useDingbats=FALSE,height=5)

p = ggplot(within(p50.summary,{taxon=factor(taxon,levels=c('human','baboon','gelada')) }),aes(set,p50.mean,fill=taxon)) +
	geom_bar(stat='identity',position='dodge') +
	geom_errorbar(aes(ymin=lb,ymax=ub),position=position_dodge2(width = 0.25, padding = 0.75)) +
	coord_flip() +
	theme_classic(base_size=24) +
	scale_fill_manual(values=c('#396caf','#81c77f','#beaed4')[3:1]) +
	scale_y_continuous(breaks=seq(0,15,5)) +
	ylab(expression(italic(P['50']))) +
	theme(axis.title.y=element_blank(),legend.title = element_blank(),legend.position='none')
ggsave(p,file='figures/gelada_p50_nolog_flip.pdf',useDingbats=FALSE,height=4)

p = ggplot(p50,aes(set,p50,fill=taxon)) +
	geom_boxplot() +
	theme_classic() +
	scale_fill_brewer(palette='Dark2') +
	scale_y_continuous(trans='log10',limits=c(1e0,max(p50$p50)),breaks=seq(5,20,5))

