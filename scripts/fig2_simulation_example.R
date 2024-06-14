source("functions.R")

#this produces a plot from a set of simulations for different dominance and double mutant phenotype combinations. Here we read in the results of these simulations
data<-read_tsv("../results/simulations/fig2_simulation.tsv")
#These simulations can be reproduced using the example.R script. 

#First we produce plots for the initial invasion by the feminiser. 
dfFem<-filter(data, invasion=="F2invasion", (F2dom==1 & M2dom==1) | (F2dom==0 & M2dom==1), F2epidom==1)
domFemInv<-invasion_sim_plot(filter(dfFem, F2dom==1), type="feminiser")
recFemInv<-invasion_sim_plot(filter(dfFem, F2dom==0), type="feminiser")

#Then, we produce the plots for invasion by the masculiniser
dfMasc<-filter(data, invasion=="M2invasion")
doubleMutMaleInv<-invasion_sim_plot(filter(dfMasc, F2epidom==0, M2epidom==1), type="masculiniser")
doubleMutFemaleInv<-invasion_sim_plot(filter(dfMasc, F2epidom==1, M2epidom==0), type="masculiniser")

#Using a shared legend
legend_grb<-get_legend(doubleMutFemaleInv, position = "top")

#arranging the panels
fig2<-ggarrange(domFemInv, doubleMutMaleInv , recFemInv, doubleMutFemaleInv, ncol=2, nrow=2, common.legend = TRUE, legend="top", labels=c("a","c","b","d"), widths=c(1,4), legend.grob=legend_grb)

#outputting the figure in different formats. 
output_name<-"fig2"
output_formats=c(".pdf", ".png", ".eps", ".jpeg", ".tiff")
lapply(output_formats, function(x)(ggsave(filename=paste0(plot_dir, output_name,x), plot=fig2, width=7.5, height=3.75)))
