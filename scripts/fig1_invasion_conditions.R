#load functions and libraries used
source("functions.R")

### This figure shows the regions for invasion by feminisers and masculinisers (into cosexuality and gynodioecy). 

#First, create a parameter grid of selfing and inbreeding depression parameters, theta and delta.
theta_delta_parameters<-expand.grid(theta=seq(0,1,0.001), delta=seq(0,1,0.001))

#now get the compensation required for invasion by feminisers and masculinisers (into cosexuality and gynodioecy) across these parameters. 
feminiser_inv_data<-mutate(theta_delta_parameters, minimum_compensation=feminiser_invasion_condition(theta, delta), type="feminiser into cosex, condition (1)")
masc_inv_data<-mutate(theta_delta_parameters, minimum_compensation=masculiniser_invasion_condition(theta, delta), type="masculiniser into cosex, condition (2)") %>% filter(!is.na(minimum_compensation))
masc_into_gyno_data<-mutate(theta_delta_parameters, minimum_compensation=masc_into_gynodioecy2(0.5, theta, delta), type="masculiniser into gynodioecy, condition (3)") %>% filter(!is.na(minimum_compensation)) 
invasion_data<-rbind(feminiser_inv_data, masc_inv_data, masc_into_gyno_data)

#brks gives the contour lines to be used
brks=seq(0,2,0.2)

#we are also going to add a line for the invasion by k_f=0.5 on two panels, which requires the following data. 
tmp<-seq(0,1,0.01)
line_data1<-tibble(theta=tmp, K=NA, k=NA) %>% 
	mutate(delta= feminiser_invasion_condition_delta(0.5, theta), type="feminiser into cosex, condition (1)", minimum_compensation=NA) %>%
	filter(0<delta & delta<1)
line_data2<-tibble(theta=tmp, K=NA, k=NA) %>% mutate(delta= feminiser_invasion_condition_delta(0.5, theta), type="masculiniser into gynodioecy, condition (3)", minimum_compensation=NA) %>%
	filter(0<delta & delta<1)
line_data=rbind(line_data1, line_data2)

#we also add a point to show the parameters used for figure 2, here we create the data to plot this point. 
point_data<-tibble(theta=0.7, delta=0.8, K=NA, minimum_compensation=NA, type="masculiniser into gynodioecy, condition (3)")


#produce plot
fig1<-ggplot(invasion_data, aes(x=theta, y=delta, z=minimum_compensation)) + #set up main data and aesthetics
	geom_contour_filled(breaks=c(-1,brks)) + #add contours from data
	geom_point(data=point_data, aes(x=theta, y=delta), colour="white") + #add point of parameters used for figure 2
	geom_line(data=line_data, aes(x=theta, y=delta), colour="black", linetype="dashed") + #add line for invasion by k_f=0.5
	scale_x_continuous(name=expression("selfing rate,"~theta), limits=c(0,1), labels = function(x) ifelse(x == 0, "0", x)) + #set axis name and remove trailing zeroes
	scale_y_continuous(name=expression("inbreeding depression,"~delta), limits=c(0,1), labels = function(x) ifelse(x == 0, "0", x)) + #set axis name and remove trailing zeroes
	scale_fill_viridis_d(name="minimum compensation required for invasion", labels=brks) + #colour scheme for contours, name and labels
	facet_wrap(~type) + #"type" used for the different invasion types, displayed across panels
	basic_theme_border + #theme 
	theme(legend.position=c(0.995,0.01), legend.justification=c(1,0), legend.direction="vertical", legend.text.align=0, legend.spacing.x = unit(0, 'cm'), legend.spacing.y = unit(0.1, 'cm'), legend.box.background = element_rect(color = "black", fill="white")) + #adjustments to legend location
	guides(fill = guide_legend(nrow = 1, label.position="bottom", label.hjust=0.5, title.hjust=0.5)) #adjustments to legend

#save in different file formats
output_name<-"fig1"
output_formats=c(".pdf", ".png", ".eps", ".jpeg", ".tiff")
lapply(output_formats, function(x)(ggsave(filename=paste0(plot_dir, output_name,x), plot=fig1, height=3.6, width=9.45)))

