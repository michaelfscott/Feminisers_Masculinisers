#load functions and libraries used
source("functions.R")

#Define parameters 
startFreqF2=1/1000 #startFreqF2, frequency that the the feminising allele is introduced at
startFreqM2=1/1000 #startFreqM2, frequency that the the masculinising allele is introduced at (assumed to be in linkage equilibrium with the feminising allele after the feminising allele has reach an equilibrium frequency)
end=0.0001 #end, maximum change in allele frequency at which to end simulation, i.e., equilibrium has been reached
Rec=0.5 #Rec, recombination rate between feminising and masculinising alleles
k=0.5 #k, relative increase in female fitness of females compared to cosexuals
KK=1 #KK, relative increase in male fitness of males compared to cosexuals
theta=0.7 #theta, selfing rate among cosexuals
delta=0.8 #delta, inbreeding depression for selfed offspring
F2dom=0 #F2dom, dominance of the feminising allele (valid values 0 or 1)
M2dom=1 #M2dom, dominance of the masculinising allele (valid values 0 or 1)
F2epidom=0 #F2epidom, epistatic dominance of the feminising allele (valid values 0 or 1)
M2epidom=1 #M2epidom, epistatic dominance of the masculinising allele (valid values 0 or 1)
min_generations=100 #min_generations, minimum number of generations to run before the simulation ends
max_generations=5000 #max_generations, maximum number of generations before the simulation ends

#run the simulation for the invasion of a feminiser, which will reach equilibrium before a masculiniser is introduced
res=runSim_Fem_Masc(startFreqF2, startFreqM2, end, Rec, k, KK, theta, delta, F2dom, M2dom, F2epidom, M2epidom, min_generations, max_generations)

#produce plot
plt<-sex_and_allele_plot(res)

#save plot
ggsave(filename=paste0(plot_dir,"example.pdf"), plot=plt, width=6, height=3)

###to produce the simulations used for figure 2, uncomment the following 
#dom_comb<-rbind(data.frame(expand.grid(F2dom=c(0,1), M2dom=c(0,1)), F2epidom=0, M2epidom=1), data.frame(expand.grid(F2dom=c(0,1), M2dom=c(0,1)), F2epidom=1, M2epidom=0)) #first get the eight combinations of dominance and double mutant phenotype
#res=do.call(rbind,lapply(1:nrow(dom_comb), function(x)(runSim_Fem_Masc(startFreqF2, startFreqM2, end, Rec, k, KK, theta, delta, F2dom=dom_comb[x,1], M2dom=dom_comb[x,2], F2epidom=dom_comb[x,3], M2epidom=dom_comb[x,4], min_generations, max_generations)))) #run simulation for each combination and combine output tables
#write_tsv(res, file="../results/simulations/fig2_simulation.tsv") #write output