#pacakges used (need to be installed first)
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(ggtext)

#define output directory 
plot_dir<-"../figures/"

#define some plotting parameters
text_size=9
#colour scheme
cp=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#e6ab02")

###################
### Functions of invasion conditions for feminisers and masculinisers
###################

feminiser_invasion_condition<-function(theta, delta){
	2*(1-theta*delta)-1
}
masculiniser_invasion_condition<-function(theta, delta){
	2*(1-theta*delta)*((1)/(1-theta))-1
}

#This function returns the masculiniser invasion condition outside the region where gynodioecy evolves
masc_into_gynodioecy<-function(k, theta, delta){
	ifelse(feminiser_invasion_condition(theta, delta)-k<0,
		(1+(1-2*delta)*theta)/(k-(1-2*delta)*theta),
		masculiniser_invasion_condition(theta, delta))
}

#This function returns the masculiniser invasion condition but only applies where gynodioecy evolves. 
masc_into_gynodioecy2<-function(k, theta, delta){
	ifelse(feminiser_invasion_condition(theta, delta)-k<0,
		(1+(1-2*delta)*theta)/(k-(1-2*delta)*theta),
		NA)
}

#This function is used when the double mutant is sterile and both mutations are dominant
masc_into_gynodioecy_dom_dom<-function(k, theta, delta){
	ifelse(feminiser_invasion_condition(theta, delta)-k<0,
		((1+k)+2*(1-delta)*theta)/((1+k)-2*(1-delta)*theta),
		masculiniser_invasion_condition(theta, delta))
}

#This function is used to plot the contour for invasion in terms of inbreeding depression (delta) for a specified k and theta. Used to add lines to the invasion plot
feminiser_invasion_condition_delta<-function(k, theta){
	(1-k)/(2*theta)
}


#Not used in paper, equilibrium frequencies of females and males
gynodioecy_female_freq<-function(k, theta, delta){
	ifelse(feminiser_invasion_condition(theta, delta)-k<0,
	(k+2*theta*delta-1)/(2*(k+theta*delta)),
	NA)
}
androdioecy_male_freq<-function(KK, theta, delta){
	ifelse(masculiniser_invasion_condition(theta, delta)-KK<0,
	((1+KK)*(1-theta)-2*(1-delta*theta))/(2*(KK)*(1-delta*theta)),
	NA)	
}


# Get Rcrit for the case of a recessive feminiser and a dominant masculiniser with double mutant sterility. 
acoef<-function(k, theta, delta){((1-theta)*(k+delta*theta))/((1+k)^2)}
bcoef<-function(k, theta, delta){(-2+2*k+theta+3*delta*theta)/(2*(1+k))}
ccoef<-function(k, theta, delta){(1-k-2*delta*theta)/(2*(k+delta*theta))}
Y1<-function(k, theta, delta){
	(-bcoef(k, theta, delta)+sqrt((bcoef(k, theta, delta))^2-4*acoef(k, theta, delta)*ccoef(k, theta, delta)))/(2*acoef(k, theta, delta))
}
Y2<-function(k, theta, delta){
	(-bcoef(k, theta, delta)-sqrt((bcoef(k, theta, delta))^2-4*acoef(k, theta, delta)*ccoef(k, theta, delta)))/(2*acoef(k, theta, delta))
}
Z<-function(k, theta, delta){
	(-1+k+2*delta*theta)/(2*(k+delta*theta))
}
aparam<-function(k, theta, delta,KK){
	(((1-Y1(k, theta, delta)-Z(k, theta, delta))+(Y1(k, theta, delta)/2))*(1-theta)*(1+KK))/((1-Z(k, theta, delta))*(1+k))
}
bparam<-function(k, theta, delta,KK){
	(((Y1(k, theta, delta)/2)*(1-theta)+Z(k, theta, delta)*(1+k))*(1+KK))/((1-Z(k, theta, delta))*(1+k))
}
Rcrit<-function(k, theta, delta,KK){
	(aparam(k, theta, delta,KK)+bparam(k, theta, delta,KK)-1)/(aparam(k, theta, delta,KK)+bparam(k, theta, delta,KK)-((aparam(k, theta, delta,KK)*bparam(k, theta, delta,KK))/(1-aparam(k, theta, delta,KK))))
}


###################
###The following functions are used to conduct simulations. 
###################

#function to convert p (frequency of A2 allele), the inbreeding coefficient (usually called F), and F2 and M2 (frequency of the feminising and masculinising mutations). 
startgen_strings <- read_file("../results/recursions/TwoLocusStartGen.txt")
inputStartGen<-str_remove(startgen_strings, "InputForm\\[\\{") %>% str_remove("\\}\\]$") %>% str_split(pattern=", ") %>% unlist()
startgen<-function(F2, M2){
	as.numeric(lapply(inputStartGen, function(x)(eval(parse(text=x)))))
}

#functions to get genotype from haplotype number
FGeno<-function(x){ceiling((x-1)/1)%%2+1}
MGeno<-function(x){ceiling((x-2)/2)%%2+1}

#A genotype matrix is convenient to convert between diploid genotype frequencies and allele frequencies. 
#This creates "gen", which has the upper triangular laid out according to the 10 diploid genotype frequencies. Only the upper triangular part is needed because we assume no parent-of-origin effects. 
create_geno_mat<-function(size){
	genoMat <- matrix(NA, size, size)
	i.low <- which(lower.tri(genoMat, diag = TRUE), arr.ind=TRUE)
	genoMat[i.low] <- 1:((size*size/2)+size/2)
	gen=t(genoMat)
	gen
}

#functions to get allele frequencies from diploid genotype frequency vectors
freqF2<-function(x){
	geno<-create_geno_mat(4)
	F2gens<-c(geno[sapply(1:4, FGeno)==2,], geno[,sapply(1:4, FGeno)==2])
	F2gens<-sort(F2gens[!is.na(F2gens)])
	sum(x[F2gens]/2)} 
freqM2<-function(x){
	geno<-create_geno_mat(4)
	M2gens<-c(geno[sapply(1:4, MGeno)==2,], geno[,sapply(1:4, MGeno)==2])
	M2gens<-sort(M2gens[!is.na(M2gens)])
	sum(x[M2gens]/2)} 

#function to run simulation
runSim<-function(sta, end=0.0001, Rec, k, KK, theta, delta, F2dom, M2dom, F2epidom, M2epidom, min_generations=100, max_generations=5000){
	#read in recursions outputted from Mathematica
	recursions_strings <- read_file(paste0("../results/recursions/TwoLocusRecursions.txt"))
	inputRecursions<-str_remove(recursions_strings, "InputForm\\[\\{") %>% str_remove("\\}\\]$") %>% str_split(pattern=", ") %>% unlist()
	recursions<-function(pop, Rec, k, KK, theta, delta, F2dom, M2dom, F2epidom, M2epidom){
		as.numeric(lapply(inputRecursions, function(x)(eval(parse(text=x)))))
	}
	#initialise population with starting frequencies
	pop<-list()
	pop[[1]]<-sta
	#record the difference in allele frequencies (and two generational difference) to determine simulation end point
	diffF2=c(0,0)
	diffM2=c(0,0)
	doublediffF2=abs(diffF2[2])-abs(diffF2[1])
	doublediffM2=abs(diffM2[2])-abs(diffM2[1])
	generation=1
	while(!(
		diffF2[2]<end && doublediffF2<=0 && generation>min_generations+1 && 
		diffM2[2]<end && doublediffM2<=0 && generation>min_generations+1 && 
		generation<max_generations+1)){
			#get next generation genotype frequencies using recursions
		generation=generation+1
		pop[[generation]]<-recursions(pop[[generation-1]], Rec, k, KK, theta, delta, F2dom, M2dom, F2epidom, M2epidom)
		diffF2=c(diffF2[2],freqF2(pop[[generation]])-freqF2(pop[[generation-1]]))
		doublediffF2=abs(diffF2[2])-abs(diffF2[1])
		diffM2=c(diffM2[2],freqM2(pop[[generation]])-freqM2(pop[[generation-1]]))
		doublediffM2=abs(diffM2[2])-abs(diffM2[1])
	}
	popmat<-do.call(rbind, pop)
	data.frame(generation=seq(1,nrow(popmat)), popmat, freqF2=apply(popmat,1,freqF2), freqM2=apply(popmat,1,freqM2), Rec=Rec, k=k, KK=KK, theta=theta, delta=delta, F2dom=F2dom, M2dom=M2dom, F2epidom=F2epidom, M2epidom=M2epidom, min_generations=min_generations, max_generations=max_generations, stringsAsFactors=FALSE)
}

#run the simulation for introducing the masculiniser at a specified frequency, then introducing the masculiniser. 
runSim_Fem_Masc<-function(startFreqF2, startFreqM2, end=0.0001, Rec, k, KK, theta, delta, F2dom, M2dom, F2epidom, M2epidom, min_generations=100, max_generations=5000){
	#first introduce the feminising allele 
	tblF2<-runSim(startgen(startFreqF2,0), end=end, Rec=Rec, k=k, KK=KK, theta=theta, delta=delta, F2dom=F2dom, M2dom=M2dom, F2epidom=F2epidom, M2epidom=M2epidom, min_generations=min_generations, max_generations=max_generations) %>% addSexFreq()

	#we then take the genotype frequencies at the end of this simulation and use this as the starting frequencies for introducing the masculinising allele. 
	starting_frequencies<-as.numeric(tblF2[nrow(tblF2),grep("X[0-9]", colnames(tblF2))])

	#Introduces the masculinising allele in linkage equilibrium with the feminising allele
	tblM2<-runSim(introduce_M2(starting_frequencies, startFreqM2), end=end, Rec, k=k, KK=KK, theta=theta, delta=delta, F2dom=F2dom, M2dom=M2dom, F2epidom=F2epidom, M2epidom=M2epidom, min_generations=min_generations, max_generations=max_generations) %>% addSexFreq()
	
	#output the result of both invasions. Also add the frequency of the sexes and label the tables. 
	mutate(rbind(
		mutate(tblF2, invasion="F2invasion"), 
		mutate(tblM2, invasion="M2invasion")), 
		startFreqF2= startFreqF2, startFreqM2=startFreqM2)
}

#function to introduce M2 allele given a vector of starting haplotype frequencies where M1 is fixed. 
introduce_M2<-function(starting_frequencies, M2freq){
	#genotype matrix
	geno<-create_geno_mat(4)
	#get those parts that correspond to M1 homozygotes and put in the corresponding frequencies from the starting frequencies vector
	M1hom<-c(geno[sapply(1:4, MGeno)==1,sapply(1:4, MGeno)==1])
	M1homfreqs<-starting_frequencies[M1hom]
	#Now remove M2freq from these M1 homozygote frequencies
	geno2<-matrix(nrow=4, ncol=4)
	geno2[sapply(1:4, MGeno)==1, sapply(1:4, MGeno)==1] <- M1homfreqs*((1-M2freq)^2)

	#the following is to add heterozgygotes, we must transpose to get ij and ji heterozygotes. 
	temp<-matrix(nrow=2,ncol=2)
	temp[1:2,1:2]<-M1homfreqs*(M2freq*(1-M2freq))
	temp[is.na(temp)]<-0
	
	geno2[sapply(1:4, MGeno)==1, sapply(1:4, MGeno)==2]<-temp+t(temp) #add temp and the transpose
	geno2[sapply(1:4, MGeno)==2, sapply(1:4, MGeno)==2]<-M1homfreqs*((M2freq)^2)

	newstart<-c(t(geno2)) #must transpose to get the genotypes in the correct order
	newstart<-newstart[!is.na(newstart)]
	newstart
}

#adds the sex frequencies to the genotype frequencies table according to the dominance and epistasis between mutations
addSexFreq<-function(df){
	Fgen1=sapply(1:4, FGeno)==1
	Mgen1=sapply(1:4, MGeno)==1
	
	geno<-create_geno_mat(4)

	frequencies<-df[,grep("X[0-9]", colnames(df))]
	
	if(df$F2dom[1]==1){ Feffect<-matrix(!as.logical(outer(Fgen1,Fgen1)),4)}
	if(df$M2dom[1]==1){ Meffect<-matrix(!as.logical(outer(Mgen1,Mgen1)),4)}

	if(df$F2dom[1]==0){ Feffect<-matrix(as.logical(outer(!Fgen1,!Fgen1)),4)}
	if(df$M2dom[1]==0){ Meffect<-matrix(as.logical(outer(!Mgen1,!Mgen1)),4)}
	
	cosex<-sort(c(geno[!Feffect & !Meffect]))
	cosex<-cosex[!duplicated(cosex)]
	females<-sort(c( geno[Feffect & (!Meffect)] ))
	females<-females[!duplicated(females)]
	males<-sort(c( geno[Meffect & (!Feffect)] ))
	males<-males[!duplicated(males)]
	sterile<-c()
	
	if(df$M2epidom[1]==1){
		males<-c(males,sort(c( geno[Feffect & Meffect] )))
		males<-sort(males[!duplicated(males)])}
	if(df$F2epidom[1]==1){
		females<-c(females,sort(c( geno[Feffect & Meffect] )))
		females<-sort(females[!duplicated(females)])}
	if(df$F2epidom[1]==0 && df$M2epidom[1]==0){
		sterile<-c(sterile,sort(c( geno[Feffect & Meffect] )))
		sterile<-sort(sterile[!duplicated(sterile)])}

	df$cosex<-apply(data.frame(frequencies[, cosex]),1,sum)
	df$females<-apply(data.frame(frequencies[, females]),1,sum)
	df$males<-apply(data.frame(frequencies[, males]),1,sum)
	df$sterile<-apply(data.frame(frequencies[, sterile]),1,sum)
	
	df
}

###################
### These are some functions used for plotting
###################

basic_theme<-theme_minimal() + 
	theme(
		legend.title=element_text(size=text_size),
	    panel.grid.minor = element_blank(),
	    text=element_text(size=text_size),
    	axis.title= element_text(size=text_size), 
    	axis.text=element_text(size=text_size),
	    strip.text = element_text(size=text_size),
	    legend.box.margin=margin(0,0,0,0), 		
		legend.margin=margin(0,0.2,0,0.2, unit="cm")) 

basic_theme_border<-basic_theme +
	theme(panel.border = element_rect(colour = "black", fill=NA)) 

basic_contour_plt<-function(dat, xvar="theta", yvar="delta", zvar="k", steps=seq(0,2,0.2), begin=0, end=1){
	ggplot(dat, aes(x=.data[[xvar]], y= .data[[yvar]], z=.data[[zvar]])) +
	geom_contour_filled(breaks=c(-1,steps)) +
	scale_x_continuous(name=expression("selfing rate,"~theta), limits=c(0,1)) +
	scale_y_continuous(name=expression("inbreeding depression,"~delta), limits=c(0,1)) +
	scale_fill_viridis_d(name=zvar, labels=steps, end=end) +
	basic_theme_border
}

#produce a plot of sex and allele frequencies from the output of a single simulation
sex_and_allele_plot<-function(df){
	#first we select the relevant columns to plot the sex frequencies
	sexfreqs<-gather(df[,c("generation", "invasion", "males", "females", "cosex","sterile")], key="sex", value="frequency", -c("generation", "invasion"))
	sexfreqs$sex<-factor(sexfreqs$sex, levels=c("cosex", "females", "males","sterile"))
	#then we use another table to plot the allele frequencies 
	allelefreqs<-gather(df[,c("generation", "invasion", "freqF2", "freqM2")], key="allele", value="frequency", -c("generation", "invasion"))
	
	ggplot() +
	#plot the sex frequencies using geom_area
	geom_area(data= sexfreqs, aes(x=generation, y=frequency, group=sex, fill=sex)) +
	#allele frequencies using geom_line
	geom_line(data=allelefreqs, aes(x=generation, y=frequency, group=allele, colour=allele), size=1) +
	#split by invasion
	facet_wrap(~invasion, scales="free") +
	#colours
	scale_fill_manual(values=c(cp[1], "grey80", "grey30", cp[4])) +
	scale_colour_manual(values=c(cp[2], cp[3])) +
	basic_theme_border

}

invasion_sim_plot<-function(df, type="feminiser"){
	
	#first we select the relevant columns to plot the sex frequencies
	sexfreqs<-gather(df[,c("generation", "males", "females", "cosex", "F2dom", "M2dom", "F2epidom", "M2epidom")], key="sex", value="frequency", -c("generation", "F2dom", "M2dom", "F2epidom", "M2epidom"))
	sexfreqs$sex<-factor(sexfreqs$sex, levels=c("cosex", "females", "males"))
	#then we use another table to plot the allele frequencies 
	allelefreqs_tmp<-gather(df[,c("generation", "F2dom", "M2dom", "F2epidom", "M2epidom", "freqF2", "freqM2")], key="allele", value="frequency", -c("generation", "F2dom", "M2dom", "F2epidom", "M2epidom"))
	
	#There are two things that change between the plots for the feminiser and the masculiniser invasions. 
	#First, the feminiser invasion plots do not include the masculiniser
	#Second, the colour of the text for the panels and legend 
	if(type=="feminiser"){
	allelefreqs<-mutate(allelefreqs_tmp, 
	frequency=ifelse(allele=="freqF2", frequency, NA),
	double_dominance=factor(ifelse(F2dom==1,
		"<b style='color:#d95f02;'>dominant </b><br>",
		"<b style='color:#d95f02;'>recessive</b><br>")),
	allele=factor(ifelse(allele=="freqF2",
		"<b style='color:#d95f02;'>feminiser</b>",
		"<b style='color:#7570b3;'>masculiniser</b>"
		), levels=c("<b style='color:#d95f02;'>feminiser</b>","<b style='color:#7570b3;'>masculiniser</b>")))
	} else {
	allelefreqs<-mutate(allelefreqs_tmp, double_dominance=factor(ifelse(F2dom==1,
	ifelse(M2dom==1,
		"<b style='color:#d95f02;'>dominant </b><br><b style='color: #7570b3;'>dominant </b>",
		"<b style='color:#d95f02;'>dominant</b><br><b style='color: #7570b3;'>recessive</b>"),
	ifelse(M2dom==1,
		"<b style='color:#d95f02;'>recessive</b><br><b style='color: #7570b3;'>dominant </b>",
		"<b style='color:#d95f02;'>recessive</b><br><b style='color: #7570b3;'>recessive</b>"))),
	allele=factor(ifelse(allele=="freqF2",
		"<b style='color:#d95f02;'>feminiser</b>",
		"<b style='color:#7570b3;'>masculiniser</b>"
		), levels=c("<b style='color:#d95f02;'>feminiser</b>","<b style='color:#7570b3;'>masculiniser</b>")))
	sexfreqs<-mutate(sexfreqs, double_dominance=factor(ifelse(F2dom==1,
	ifelse(M2dom==1,
		"<b style='color:#d95f02;'>dominant </b><br><b style='color: #7570b3;'>dominant </b>",
		"<b style='color:#d95f02;'>dominant</b><br><b style='color: #7570b3;'>recessive</b>"),
	ifelse(M2dom==1,
		"<b style='color:#d95f02;'>recessive</b><br><b style='color: #7570b3;'>dominant </b>",
		"<b style='color:#d95f02;'>recessive</b><br><b style='color: #7570b3;'>recessive</b>"))))	
	}
	
	ggplot() +
	#plot the sex frequencies using geom_area
	geom_area(data= sexfreqs, aes(x=generation, y=frequency, group=sex, fill=sex)) +
	#allele frequencies using geom_line
	geom_line(data=allelefreqs, aes(x=generation, y=frequency, group=allele, colour=allele), size=1) +
	facet_wrap(~ double_dominance, scales="free_x", nrow=1) +
	#colours
	scale_fill_manual(values=c(cp[1], "grey80", "grey30")) +
	scale_colour_manual(values=c(cp[2], cp[3])) +
	basic_theme_border +
	theme(strip.text=element_markdown(), legend.text=element_markdown(size=text_size), legend.position="top", legend.title=element_blank()) +
	guides(fill=guide_legend(order=2),
			colour=guide_legend(order=1))

}
