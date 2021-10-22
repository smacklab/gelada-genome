#!/usr/bin/env Rscript

gelada.hb = read.delim('data/gelada_hb_concentrations.txt')

library(lubridate)

gelada.hb$Date = mdy(gelada.hb$Date)

individuals = read.csv('data/individuals.csv')

individuals$Date.of.Birth = ymd(individuals$Date.of.Birth)

gelada.hb$Code = toupper(substr(gsub(' ','',gelada.hb$Name),1,3))

dobs = subset(individuals,select=c('Code','Date.of.Birth'))

gelada.hb = merge(gelada.hb,dobs,by='Code',all.x=TRUE,all.y=FALSE)

gelada.hb = gelada.hb[order(gelada.hb$ID),]

gelada.hb$age_est = with(gelada.hb,as.integer(difftime(Date,Date.of.Birth,'days'))/30)

gelada.hb$taxon = 'gelada'

gelada.hb$Hematocrit = with(gelada.hb,Hemoglobin * 2.94)


# Obtained from ISIS data
zoo.hb = data.frame(
	mean=12.8,
	sd = 1.5,
	min = 7.4,
	max=15.2,
	n=42
)

zoo.crit = data.frame(
	mean=39.1,
	sd=4.2,
	min=26,
	max=49,
	n=65
)

# Function for comparing summary stats of two populations and calculating joint summary stats
combine.samples = function(l1,l2) {
	n1 = l1$n
	n2 = l2$n
	x1 = l1$mean
	x2 = l2$mean
	s1 = l1$sd
	s2 = l2$sd

	x0 = (x1*n1+x2*n2)/(n1+n2)
	s0 = sqrt((n1*(s1^2+(x1-x0)^2) + n2*(s2^2+(x2-x0)^2))/(n1+n2))

	data.frame(
		mean=x0,
		sd=s0,
		min=min(c(l1$min,l2$min)),
		max=max(c(l1$max,l2$max)),
		n=n1+n2
	)
}

# Obtained from 10.1111/j.1600-0684.1999.tb00085.x
hamadryas.hb.males = data.frame(
	mean=12.63,
	sd=1.08,
	min=10.47,
	max=14.78,
	n=407
)
hamadryas.hb.females = data.frame(
	mean=12.52,
	sd=0.97,
	min=10.58,
	max=14.45,
	n=616
)

hamadryas.hb = combine.samples(hamadryas.hb.males,hamadryas.hb.females)

hamadryas.crit.males = data.frame(
	mean=40,
	sd=3,
	min=33,
	max=46,
	n=407
)
hamadryas.crit.females = data.frame(
	mean=40,
	sd=3,
	min=33,
	max=46,
	n=617
)

hamadryas.crit = combine.samples(hamadryas.crit.males,hamadryas.crit.females)

simiens.hb = with(gelada.hb,data.frame(
	mean=mean(Hemoglobin,na.rm=TRUE),
	sd=sd(Hemoglobin,na.rm=TRUE),
	min=min(Hemoglobin,na.rm=TRUE),
	max=max(Hemoglobin,na.rm=TRUE),
	n=sum(!is.na(Hemoglobin))
))

simiens.crit = with(gelada.hb,data.frame(
	mean=mean(Hematocrit,na.rm=TRUE),
	sd=sd(Hematocrit,na.rm=TRUE),
	min=min(Hematocrit,na.rm=TRUE),
	max=max(Hematocrit,na.rm=TRUE),
	n=sum(!is.na(Hematocrit))
))

# Function for calculating T statistic from two populations
t.compare = function(l1,l2) {
	n1 = l1$n
	n2 = l2$n
	x1 = l1$mean
	x2 = l2$mean
	s1 = l1$sd
	s2 = l2$sd

	se = sqrt( (s1^2/n1) + (s2^2/n2) )
	DF = ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
	T = (x1-x2)/se

# 	x0 = (x1*n1+x2*n2)/(n1+n2)
# 	s0 = sqrt((n1*(s1^2+(x1-x0)^2) + n2*(s2^2+(x2-x0)^2))/(n1+n2))
# 	T = (x1-x2)/(s0/sqrt(n1+n2))

#	sp = sqrt((s1^2+s2^2)/2)
#
#	T = (x1-x2)/(sp*sqrt(2/(n1+n2)))
#	DF = (n1+n2) - 2
#	pt(q=abs(T),df=DF,lower.tail=FALSE)

	out = c(x1-x2, se, T, 2*pt(-abs(T),DF)) 
	names(out) = c('Difference of means', 'Std Error', 't', 'p-value')
	out
}

# Data from 10.1002/(SICI)1096-8644(199807)106:3<385::AID-AJPA10>3.0.CO;2-X

tibet.male = data.frame(
	mean=15.6,
	sd=0.17,
	min=11.9,
	max=20,
	n=75
)
aymara.male = data.frame(
	mean=19.1,
	sd=0.18,
	min=15,
	max=25.1,
	n=91
)
tibet.female = data.frame(
	mean=14.2,
	sd=0.14,
	min=12.4,
	max=17.6,
	n=61
)
aymara.female = data.frame(
	mean=17.8,
	sd=0.23,
	min=11.8,
	max=23.6,
	n=83
)

tibet.hb = combine.samples(tibet.male,tibet.female)
aymara.hb = combine.samples(aymara.male,aymara.female)

# Data from 10.1002/ajh.25194
amhara.male = data.frame(
	mean=16.7,
	sd=1.1,
	min=NA,
	max=NA,
	n=99
)
oromo.male = data.frame(
	mean=18.7,
	sd=1.7,
	min=NA,
	max=NA,
	n=66
)
amhara.female = data.frame(
	mean=15.5,
	sd=1.2,
	min=NA,
	max=NA,
	n=33
)
oromo.female = data.frame(
	mean=16.8,
	sd=1.6,
	min=NA,
	max=NA,
	n=35
)
amhara.hb = combine.samples(amhara.male,amhara.female)
oromo.hb = combine.samples(oromo.male,oromo.female)

# Data from 10.1182/blood-2005-07-3046
san.diego = read.delim('data/usa_hb_concentrations.txt')

san.diego$min = san.diego$max = NA

us.hb = Reduce(combine.samples,split(san.diego,1:nrow(san.diego)))

t.compare(hamadryas.hb,simiens.hb)
t.compare(zoo.hb,simiens.hb)

t.compare(hamadryas.crit,simiens.crit)
t.compare(zoo.crit,simiens.crit)

comparative.hb = rbind(
	data.frame(simiens.hb,pop='geladas (wild)',alt='high'),
	data.frame(zoo.hb,pop='geladas (zoo)',alt='low'),
	data.frame(hamadryas.hb,pop='baboons',alt='low'),
	data.frame(amhara.hb,pop='humans (Ethiopian Amhara)',alt='high'),
	data.frame(oromo.hb,pop='humans (Ethiopian Oromo)',alt='high'),
	data.frame(tibet.hb,pop='humans (Tibetans)',alt='high'),
	data.frame(aymara.hb,pop='humans (Andeans)',alt='high'),
	data.frame(us.hb,pop='humans (USA)',alt='low')
)

comparative.hb$pop = factor(comparative.hb$pop,
	levels=c('geladas (wild)','geladas (zoo)','baboons','humans (USA)','humans (Tibetans)','humans (Ethiopian Amhara)','humans (Ethiopian Oromo)','humans (Andeans)'))

library(ggplot2)

p = ggplot(aes(pop,mean,fill=alt),data=comparative.hb) +
	geom_bar(stat='identity') +
	geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2) +
	scale_fill_brewer(palette='Dark2',name='Altitude') +
	theme_classic(base_size=24) +
	theme(axis.title.x=element_blank(),axis.text.x=element_text(angle = -45, hjust = 0, vjust = 1)) +
	ylab('Hb concentration (g/dL)')


p = ggplot(aes(pop,mean),data=comparative.hb) +
	geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2,color='#000000',show.legend=FALSE) +
	geom_point(aes(color=alt),size=10) +
	geom_text(aes(label=format(round(mean,1))),color='#ffffff') +
	scale_color_brewer(palette='Dark2',name='Altitude') +
	scale_y_continuous(limits=c(10,20),breaks=seq(10,20,5)) +
	coord_flip() +
	theme_classic(base_size=24) +
	theme(axis.title.y=element_blank()) +
	guides() +
	ylab('Hb concentration (g/dL)')
ggsave(p,file='figures/hb_concentrations.pdf',height=5,width=10)





hb.gradient.raw = read.delim('data/elevation_hb_concentrations.txt')

hb.gradient = subset(hb.gradient.raw,Sex %in% c('m','f','m-adu','f-adu','f-preg'))

# hb.gradient$Category[hb.gradient$Category %in% c('Sherpa','Tibetan')] = 'Tibetan'
# hb.gradient$Category[hb.gradient$Category %in% c('Han Chinese','Asia')] = 'Asia'

hb.gradient.split = split(hb.gradient,list(hb.gradient$Category,hb.gradient$Elevation))
hb.gradient.split = hb.gradient.split[as.logical(unlist(lapply(hb.gradient.split,nrow)))]

hb.gradient.all = subset(do.call(rbind,lapply(hb.gradient.split,function(x) {
	out = unique(x[c('Human','Category','Elevation')])
	out$N = with(x,sum((!is.na(Hb)) * N,na.rm=TRUE))
	if (!out$N) {
		out$Hb = NA
	} else {
		out$Hb = with(x,sum(Hb * N,na.rm=TRUE)/sum((!is.na(Hb)) * N,na.rm=TRUE))
	}
	out
})),N>0)


hb.gradient.all$Category = factor(hb.gradient.all$Category,
	levels=c(
		'Gorilla','Hylobates','Nomascus','Pan','Pongo',
		'Tibetan','Han Chinese','Sherpa','Asia',
		'Ethiopia','Africa',
		'Quechua','Central/South America',
		'USA'
	),
	labels=c(
		'Gorilla','Hylobates','Nomascus','Pan','Pongo',
		'Tibetan','Han Chinese','Sherpa','Asia (other)',
		'Ethiopia','Africa (other)',
		'Quechua','Central/South America (other)',
		'USA'
	)
)



p = ggplot() +
	geom_smooth(data=droplevels(subset(hb.gradient.all,Human=='Human' & !Category %in% c('Sherpa','Tibetan'))),aes(Elevation/1000,Hb),method=lm,color='black',se=TRUE,size=0.5,linetype=1) +
	geom_smooth(data=droplevels(subset(hb.gradient.all,Human=='Human' & Category %in% c('Sherpa','Tibetan'))),aes(Elevation/1000,Hb),method=lm,color='black',fill='#cccccc',se=TRUE,size=0.5,linetype=2) +
	geom_point(data=droplevels(subset(hb.gradient.all,Human=='Human')),aes(Elevation/1000,Hb,color=Category),shape=19,size=2) +
#	geom_point(data=subset(hb.gradient.all,Human=='NHP'),aes(Elevation/1000,Hb,shape=Category),color='black',size=4) +
#	geom_smooth(data=subset(hb.gradient.all,Human=='NHP'),aes(Elevation/1000,Hb),method=lm,color='black',se=TRUE,size=0.5,linetype=3) +
	geom_point(data=simiens.hb,aes(3250/1000,mean),color='black',size=2) +
	geom_errorbar(data=simiens.hb,aes(x=3250/1000,ymin=mean-sd,ymax=mean+sd),color='black',size=0.5,width=0.2) +
	geom_point(data=zoo.hb,aes(100/1000,mean),color='black',size=2) +
	geom_errorbar(data=zoo.hb,aes(x=100/1000,ymin=mean-sd,ymax=mean+sd),color='black',size=0.5,width=0.2) +
#	geom_point(data=hamadryas.hb,aes(0/1000,mean),color='black',size=2) +
#	geom_errorbar(data=hamadryas.hb,aes(x=0/1000,ymin=mean-sd,ymax=mean+sd),color='black',size=0.5,width=0.2) +
#	coord_flip(xlim=c(0,6)) +
#	scale_color_brewer(palette='Accent') +
	coord_cartesian(xlim=c(0,6)) +
	scale_color_manual(values=c('#1f78b4','#a6cee3','#33a02c','#b2df8a','#e31a1c','#fb9a99','#6a3d9a','#cab2d6','#666666')) +
	scale_shape_manual(values=c('g','h','n','c','o')) +
	scale_x_continuous(limits=c(-0.1,6),breaks=seq(0,6,1)) +
	scale_y_continuous(limits=c(10,20),breaks=seq(10,20,2)) +
	theme_classic(base_size=20) +
	theme(legend.title=element_blank()) +
	xlab('Elevation (km)') + ylab('Hb concentration (g/dl)')
ggsave(p,file='figures/hb_gradient.pdf',useDingbats=FALSE,height=5,width=9)