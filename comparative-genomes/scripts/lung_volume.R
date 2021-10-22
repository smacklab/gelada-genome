#!/usr/bin/env Rscript

babdat = readRDS('checkpoints/babdat.rds')

babdat.vars = readRDS('checkpoints/babdat_vars.rds')

# subset(babdat.vars,grepl('girth',label))

# wtk = weight
# digimi = min girth
# digimx = max girth
# sex
# tx
# age_est

babdat.full = subset(babdat,select=c('ads','tx','sex','age_est','wtk','digimi','digimx'))
babdat.full = babdat.full[complete.cases(babdat.full[c('tx','sex','age_est','digimx')]),]

gelada.morph = read.delim('data/gelada_morphometrics.txt')

library(lubridate)

gelada.morph$Date = mdy(gelada.morph$Date)

individuals = read.csv('data/individuals.csv')

individuals$Date.of.Birth = ymd(individuals$Date.of.Birth)

gelada.morph$Code = toupper(substr(gsub(' ','',gelada.morph$Name),1,3))

dobs = subset(individuals,select=c('Code','Date.of.Birth'))

gelada.morph = merge(gelada.morph,dobs,by='Code',all.x=TRUE,all.y=FALSE)

gelada.morph = gelada.morph[order(gelada.morph$ID),]

gelada.morph$age_est = with(gelada.morph,as.integer(difftime(Date,Date.of.Birth,'days'))/30)

gelada.morph$Chest = apply(subset(gelada.morph,select=grepl('Chest[1-3]',names(gelada.morph))),1,function(x) {
	if (is.na(x['Chest1'])) {
		NA
	} else if (!is.na(x['Chest2'])) {
		mean(sort(c(x['Chest1'],x['Chest2'],x['Chest3']))[which.min(abs(diff(sort(c(x['Chest1'],x['Chest2'],x['Chest3']))))):(which.min(abs(diff(sort(c(x['Chest1'],x['Chest2'],x['Chest3'])))))+1)],na.rm=TRUE)
	} else {
		mean(c(x['Chest1'],x['Chest2']),na.rm=TRUE)
	}
})

gelada.morph$Waist = apply(subset(gelada.morph,select=grepl('Waist[1-3]',names(gelada.morph))),1,function(x) {
	if (is.na(x['Waist1'])) {
		NA
	} else if (!is.na(x['Waist2'])) {
		mean(sort(c(x['Waist1'],x['Waist2'],x['Waist3']))[which.min(abs(diff(sort(c(x['Waist1'],x['Waist2'],x['Waist3']))))):(which.min(abs(diff(sort(c(x['Waist1'],x['Waist2'],x['Waist3'])))))+1)],na.rm=TRUE)
	} else {
		mean(c(x['Waist1'],x['Waist2']),na.rm=TRUE)
	}
})

gelada.morph$taxon = 'gelada'

baboon.final = subset(babdat.full,tx %in% c('Hamadryas','Guinea (Chicago)','Anubis (Awash)','Yellow (Mikumi)','Moremi chacmas','Zambian  Kindas') & age_est >= 72,select=c('ads','tx','sex','wtk','digimi','digimx'))
gelada.final = subset(gelada.morph,Age == 'Adult',select=c('ID','taxon','Sex','Mass','Waist','Chest'))

names(baboon.final) = names(gelada.final) = c('id','taxon','sex','mass','waist','chest')

baboon.final$sex = substr(baboon.final$sex,1,1)

baboon.final$taxon = factor(baboon.final$taxon,levels=c('Gelada','Hamadryas','Guinea (Chicago)','Anubis (Awash)','Yellow (Mikumi)','Moremi chacmas','Zambian  Kindas'))
levels(baboon.final$taxon) = c('gelada','hamadryas','Guinea','anubis','yellow','chacma','Kinda')
gelada.final$taxon = factor(gelada.final$taxon,levels=levels(baboon.final$taxon))

papionin.final = rbind(baboon.final,gelada.final)

papionin.final = papionin.final[complete.cases(papionin.final[c('taxon','sex','chest')]),]

# Drop outliers
papionin.final = subset(papionin.final,!(!is.na(id) & id %in% c('DRT_2017_001','DRT_2019_001')))

# For now, drop chacmas
papionin.final = droplevels(subset(papionin.final,taxon != 'chacma'))

papionin.final$sex = factor(papionin.final$sex,levels=c('F','M'),labels=c('female','male'))

library(ggplot2)
library(ggbeeswarm)
library(RColorBrewer)

p = ggplot(papionin.final,aes(x=taxon,y=chest/waist)) +
	geom_boxplot(aes(alpha=sex),width=0.1,size=0.5,outlier.colour=NA,position=position_dodge(width = 0.8),color='#cccccc',show.legend=FALSE) +
	geom_quasirandom(aes(color=sex),dodge.width=0.8) +
	scale_color_brewer(palette='Dark2') +
	scale_alpha_manual(values=c(1,1)) +
	guides(alpha=guide_legend(show=FALSE),color=guide_legend(override.aes=list(size=2))) +
	theme_classic(base_size=24) +
	theme(axis.title.x=element_blank(),legend.title = element_blank(),axis.text.x=element_text(angle = -45, hjust = 0, vjust = 1)) +
	ylab('Chest:waist ratio')
ggsave(p,file='figures/gelada_chest_waist_ratio.pdf',height=6,useDingbats=FALSE)

p = ggplot(papionin.final,aes(x=taxon,y=chest/mass)) +
	geom_boxplot(aes(alpha=sex),width=0.1,size=0.5,outlier.colour=NA,position=position_dodge(width = 0.8),color='#cccccc',show.legend=FALSE) +
	geom_quasirandom(aes(color=sex),dodge.width=0.8) +
	scale_color_brewer(palette='Dark2') +
	scale_alpha_manual(values=c(1,1)) +
	guides(alpha=guide_legend(show=FALSE),color=guide_legend(override.aes=list(size=2))) +
	theme_classic(base_size=24) +
	theme(axis.title.x=element_blank(),legend.title = element_blank(),axis.text.x=element_text(angle = -45, hjust = 0, vjust = 1)) +
	ylab('Chest / body mass (cm/kg)')
ggsave(p,file='figures/gelada_chest_weight_ratio.pdf',height=6,useDingbats=FALSE)

# Custom colors
taxon.colors = c('#386cb0','#beaed4','#fdc086','#7fc97f','#ffff99','#f0027f')

p = ggplot(papionin.final,aes(x=mass,y=chest/waist)) +
#	geom_boxplot(aes(alpha=sex),width=0.1,size=0.5,outlier.colour=NA,position=position_dodge(width = 0.8),color='#cccccc',show.legend=FALSE) +
#	geom_quasirandom(aes(color=sex),dodge.width=0.8) +
	geom_point(aes(color=taxon,shape=sex)) +
#	scale_color_brewer(palette='Accent') +
	scale_color_manual(values=taxon.colors) +
	scale_alpha_manual(values=c(1,1)) +
	guides(alpha=guide_legend(show=FALSE),color=guide_legend(override.aes=list(size=2))) +
	theme_classic(base_size=24) +
	theme(legend.title = element_blank()) +
	ylab('Chest:waist ratio') + xlab('Body mass')
ggsave(p,file='figures/gelada_chest_waist_ratio_vs_weight.pdf',height=6,useDingbats=FALSE)


p = ggplot(papionin.final,aes(x=mass,y=chest/mass)) +
#	geom_boxplot(aes(alpha=sex),width=0.1,size=0.5,outlier.colour=NA,position=position_dodge(width = 0.8),color='#cccccc',show.legend=FALSE) +
#	geom_quasirandom(aes(color=sex),dodge.width=0.8) +
	geom_point(aes(color=taxon,shape=sex)) +
#	scale_color_brewer(palette='Accent') +
	scale_color_manual(values=taxon.colors) +
	scale_alpha_manual(values=c(1,1)) +
	guides(alpha=guide_legend(show=FALSE),color=guide_legend(override.aes=list(size=2))) +
	theme_classic(base_size=24) +
	theme(legend.title = element_blank()) +
	ylab('Chest / body mass (cm/kg)') + xlab('Body mass')
ggsave(p,file='figures/gelada_chest_weight_ratio_vs_weight.pdf',height=6,useDingbats=FALSE)


papionin.final = within(papionin.final,{
	chest_waist = chest/waist
	chest_mass = chest/mass
	genus = c('Theropithecus','Papio')[as.integer(factor(taxon == 'gelada',levels=c('TRUE','FALSE')))]
})

summary(lm(chest_waist~taxon + sex,data=papionin.final))
summary(lm(chest_waist~genus + genus:taxon + sex,data=papionin.final))
summary(lm(chest_mass~taxon + sex,data=papionin.final))
summary(lm(chest_mass~genus + genus:taxon + sex,data=papionin.final))

# Significant when we drop Kindas
summary(lm(chest_mass~genus + genus:taxon + sex,data=droplevels(subset(papionin.final,taxon!='Kinda'))))

# Run all of the above with mass as a covariate

summary(lm(chest_waist~taxon + sex + mass,data=papionin.final))
summary(lm(chest_waist~genus + genus:taxon + sex + mass,data=papionin.final))
summary(lm(chest_mass~taxon + sex + mass,data=papionin.final))
summary(lm(chest_mass~genus + genus:taxon + sex + mass,data=papionin.final))


res = lm(chest_waist~genus + genus:taxon + sex + mass,data=papionin.final)
coef(summary(res))['genusTheropithecus','Estimate']
pt(coef(summary(res))['genusTheropithecus','t value'],df=res$df.residual,lower.tail=FALSE)
res = lm(chest_mass~genus + genus:taxon + sex + mass,data=papionin.final)
coef(summary(res))['genusTheropithecus','Estimate']
pt(coef(summary(res))['genusTheropithecus','t value'],df=res$df.residual,lower.tail=FALSE)



# Final model
res = lm(chest~genus + genus:taxon + sex + mass + waist,data=papionin.final)
coef(summary(res))['genusTheropithecus','Estimate']
pt(coef(summary(res))['genusTheropithecus','t value'],df=res$df.residual,lower.tail=FALSE)

res = lm(chest~genus + genus:taxon + sex + waist,data=papionin.final)
coef(summary(res))['genusTheropithecus','Estimate']
pt(coef(summary(res))['genusTheropithecus','t value'],df=res$df.residual,lower.tail=FALSE)

res = lm(chest~genus + genus:taxon + sex + mass,data=papionin.final)
coef(summary(res))['genusTheropithecus','Estimate']
pt(coef(summary(res))['genusTheropithecus','t value'],df=res$df.residual,lower.tail=FALSE)



p1 = ggplot(aes(mass,chest),data=papionin.final) +
	geom_point(aes(color=taxon,shape=genus)) +
#	geom_smooth(aes(color=taxon,size=genus),method=lm,se=FALSE,show.legend=FALSE) +
	geom_smooth(data=subset(papionin.final,genus=='Theropithecus'),aes(mass,chest),method=lm,color=taxon.colors[1],se=FALSE,show.legend=FALSE) +
	geom_smooth(data=subset(papionin.final,genus=='Papio'),aes(mass,chest),method=lm,color='#000000',linetype=2,size=0.5,se=FALSE,show.legend=FALSE) +
	scale_color_manual(values=taxon.colors) +
	scale_shape_manual(values=c(21,19)) +
#	scale_size_manual(values=c(0.5,1)) +
	scale_y_continuous(limits=c(40,80),breaks=seq(40,80,10)) +
	theme_classic(base_size=30) +
	theme(legend.title = element_blank(),legend.position='none') +
	guides(shape=FALSE,color=guide_legend(override.aes=list(size=2))) +
	ylab('Chest circ. (cm)') +
	xlab('Body weight (kg)')

p2 = ggplot(aes(waist,chest),data=papionin.final) +
	geom_point(aes(color=taxon,shape=genus)) +
#	geom_smooth(aes(color=taxon,size=genus),method=lm,se=FALSE,show.legend=FALSE) +
	geom_smooth(data=subset(papionin.final,genus=='Theropithecus'),aes(waist,chest),method=lm,color=taxon.colors[1],se=FALSE,show.legend=FALSE) +
	geom_smooth(data=subset(papionin.final,genus=='Papio'),aes(waist,chest),method=lm,color='#000000',linetype=2,size=0.5,se=FALSE,show.legend=FALSE) +
	scale_color_manual(values=taxon.colors) +
	scale_shape_manual(values=c(21,19)) +
#	scale_size_manual(values=c(0.5,1)) +
	scale_y_continuous(limits=c(40,80),breaks=seq(40,80,10)) +
	theme_classic(base_size=30) +
	theme(
		axis.text.y = element_blank(),
		axis.title.y = element_blank(),
		axis.ticks.y = element_blank(),
		legend.title = element_blank(),
		legend.position='none'
	) +
	guides(shape=FALSE,color=guide_legend(override.aes=list(size=2))) +
	ylab('Chest circ. (cm)') +
	xlab('Waist circ. (cm)')

p3 = ggplot(aes(NA,NA),data=papionin.final) +
	geom_point(aes(color=taxon,shape=genus),alpha=0) +
#	geom_smooth(aes(color=taxon,size=genus),method=lm,se=FALSE,show.legend=FALSE) +
	scale_color_manual(values=taxon.colors) +
	scale_shape_manual(values=c(21,19)) +
#	scale_size_manual(values=c(0.5,1)) +
	coord_fixed(ratio=100) +
	theme_classic(base_size=30) +
	theme(
		axis.title=element_blank(),
		axis.text=element_blank(),
		axis.ticks=element_blank(),
		legend.title = element_blank(),
		axis.line=element_blank()
	) +
	guides(shape=FALSE,color=guide_legend(override.aes=list(size=2,alpha=1)))

library(egg)

pdf(file='figures/gelada_chest_combined.pdf',useDingbats=FALSE,height=7,width=11)
	ggarrange(p1,p2,p3,ncol=3,nrow=1,widths=c(1.25,1,1),newpage=FALSE)
dev.off()