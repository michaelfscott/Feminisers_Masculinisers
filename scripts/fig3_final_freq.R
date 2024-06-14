source("functions.R")

#read in the data, which is the final genotype frequency for different parameter combinations
tbl<-read_tsv("../results/simulations/fig3_final_frequencies.tsv")
#re-arrange table to give the allele frequencies
plt_dat<-gather(tbl[,c("feminiser_dominance", "masculiniser_dominance", "feminiser", "masculiniser", "epistatic_dominance", "fem_inv", "masc_inv", "generation")], key="allele", value="frequency", -c("feminiser_dominance", "masculiniser_dominance", "epistatic_dominance", "fem_inv", "masc_inv", "generation"))  %>%
	mutate(allele=factor(ifelse(allele=="feminiser", #change the names to write them in colour on plot
	"<b style='color:#d95f02;'>feminiser</b>",
	"<b style='color:#7570b3;'>masculiniser</b>"
	), levels=c("<b style='color:#d95f02;'>feminiser</b>","<b style='color:#7570b3;'>masculiniser</b>")))

#produce plot
fig3<-plt_dat %>%
	ggplot(aes(x=paste0(feminiser_dominance, "\n", masculiniser_dominance), y=frequency, colour=allele)) +
	geom_point(position = position_jitter(seed = 1234), alpha=0.25, pch=16) +
	facet_wrap(~ epistatic_dominance) +
	ylim(0,1.0000001) +
	ylab("final allele frequency") +
	scale_color_manual(values=cp[c(2,3)]) +
	basic_theme_border +
	guides(col = guide_legend(override.aes = list(shape = 15, size = 5, alpha=1))) +
	scale_x_discrete(labels=c(
		"<b style='color:#d95f02;'>dominant </b><br><b style='color: #7570b3;'>dominant </b>",
		"<b style='color:#d95f02;'>dominant</b><br><b style='color: #7570b3;'>recessive</b>",
		"<b style='color:#d95f02;'>recessive</b><br><b style='color: #7570b3;'>dominant </b>",
		"<b style='color:#d95f02;'>recessive</b><br><b style='color: #7570b3;'>recessive</b>")) +
	theme(axis.title.x=element_blank(), axis.text.x= element_markdown(), strip.text=element_markdown(), 
	legend.position=c(0.475,0.05), legend.justification=c(1,0), legend.title=element_blank(), legend.box.background=element_rect(fill="white"), legend.box.margin=margin(0,0,0,0), legend.margin=margin(-0.05,0.2,0,0, unit="cm"), legend.text=element_markdown()) 

#outputting the figure in different formats. 
output_name<-"fig3"
output_formats=c(".pdf", ".png", ".jpeg", ".tiff")
lapply(output_formats, function(x)(ggsave(filename=paste0(plot_dir, output_name,x), plot= fig3, width=6, height=2.5)))


