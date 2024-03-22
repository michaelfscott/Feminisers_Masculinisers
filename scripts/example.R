#load functions and libraries used
source("functions.R")

#Define parameters 
startFreqF2=0.01 #startFreqF2, frequency that the the feminising allele is introduced at
startFreqM2=0.01 #startFreqM2, frequency that the the masculinising allele is introduced at (assumed to be in linkage equilibrium with the feminising allele after the feminising allele has reach an equilibrium frequency)
end=0.001 #end, maximum change in allele frequency at which to end simulation, i.e., equilibrium has been reached
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