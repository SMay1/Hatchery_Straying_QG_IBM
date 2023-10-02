##Visualization for Hathcery_Straying_QG_IBM


#### 1: Load Data and Packages ####

#Set Working Directory to Source File Location

library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
options(dplyr.group_by.inform = FALSE)
library(ggplot2)
#library(sjPlot)
#library(AICcmodavg)
#library(mgcv)
#library(MASS)
library(ggpubr)
#library(TruncatedNormal)
#library(rockchalk)
#devtools::install_github("zeehio/facetscales")
library(facetscales)
library(RColorBrewer)
library(purrr)
library(furrr)
library(TruncatedNormal)


#Script to load simulation iterations, calculate summary statistics, and plot data.
#Summary statistics data provided. Individual model outputs are not provided, as they are very large.

###Because the output is very large, and because we are only using the summary statistics for each iteration
###Easier to calculate summary stats for each iter, and then combine after.
###If, in the future, one wants to load in all together - here are two different methods to do so:

#gen_data_all <- list.files(path = "../data/iters/",pattern = "*.csv") %>% 
#  map_df(~read_csv(paste("../data/iters/",.,sep="")))

#file_names = list.files(path = "../data/iters/",pattern="*.csv",full.names = T)
#gen_data_all<-do.call(rbind, lapply(file_names, read.csv, header = T))


#summary stats loop:
#Reads each output file, calculates summary stats, adds them to "summary_stats_all"


plan(multisession, workers = availableCores())

file_names <- list.files(path = "output/Experiment_4/", pattern = "*.csv", full.names = TRUE, include.dirs = FALSE)

summary_stats_all <- future_map_dfr(file_names, ~ {
  gen_data_i <- read.csv(.x)
  summary_stats_i <- gen_data_i %>%
    group_by(G, beta, rho, var_entry_day, var_RLS, PopSize, iter, var_env_entry, var_env_RLS, OMEGA, pHOS, MU_entry_day_wild, MU_entry_day_hatchery, MU_RLS_hatchery, origin, F0_hatchery, F1_hatchery, n_mates_lambda) %>%
    dplyr::summarise(
      Nc = n(),
      meanK = mean(RS_obs),
      varK = var(RS_obs),
      S = sum(RS_obs) / 2,
      Ne = ((meanK * Nc) - 2) / (meanK - 1 + (varK / meanK)),
      Ne_N = Ne / Nc,
      Nm = sum(sex == "m"),
      Nf = sum(sex == "f"),
      sex_ratio = Nm / Nf,
      mean_entry_day = mean(entry_day),
      sd_entry_day = sd(entry_day),
      mean_pHOA = mean(pHOA),
      mean_RLS = mean(RLS),
      sd_RLS = sd(RLS),
      Nb = ((meanK * Nc) - 2) / (meanK - 1 + (varK / meanK)),
      breeders = sum(RS_obs > 0),
      mean_pHOA = mean(pHOA),
      .groups = "keep"  # Set .groups argument to "keep" to avoid the grouping message
    )
  
  summary_stats_i
})


plan(sequential)

#save(list = "summary_stats_all",file = "output/Experiment_1/Experiment_1_sumstats.Rdata")
#save(list = "summary_stats_all",file = "output/Experiment_2-3/Experiment_2-3_sumstats.Rdata")
#save(list = "summary_stats_all",file = "output/Experiment_4/Experiment_4_sumstats.Rdata")
load(file = "output/Experiment_1/Experiment_1_sumstats.Rdata")
load("output/Experiment_2-3_sumstats.Rdata")
load(file = "output/Experiment_4/Experiment_4_sumstats.Rdata")


#summary_stats_all<-read.csv("summary_stats_all.csv")

#pHOS_palette <- brewer.pal(name="OrRd",n=9)[c(3,5,6,8)]
pHOS_palette <- viridis::viridis(4,direction = -1)
return_palette <- brewer.pal(name="Blues",n=9)[c(8,6,4,2)]





###Figure1 Schematic

sim_data<-data.frame(origin = c(rep("wild",100000),rep("hatchery",120000)),
                     mu_return_day = c(rep(10,100000),rep(c(10.1,15,20,25),each=30000)),
                     return_day = numeric(220000)) # make sure it's numeric

sim_data$return_day[1:100000]<-rtnorm(100000,mu=10,sd=sqrt(20),lb=0,ub=45) #generate 10,000 numbers
sim_data$return_day[100001:130000]<-rtnorm(30000,mu=10,sd=sqrt(20),lb=0,ub=45) # generate 3,000 numbers
sim_data$return_day[130001:160000]<-rtnorm(30000,mu=15,sd=sqrt(20),lb=0,ub=45) # generate 3,000 numbers
sim_data$return_day[160001:190000]<-rtnorm(30000,mu=20,sd=sqrt(20),lb=0,ub=45) # generate 3,000 numbers
sim_data$return_day[190001:220000]<-rtnorm(30000,mu=25,sd=sqrt(20),lb=0,ub=45) # generate 3,000 numbers

smoothed_data <- sim_data %>% 
  mutate(entry = as.factor(mu_return_day)) %>% 
  group_by(bin = cut_width(return_day, 1), entry) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(entry) %>%
  mutate(count_smooth = rollmean(count, k = 4, fill = NA)) %>% 
  mutate(count_smooth=count_smooth/100)

# plot the smoothed data
Fig1b<-ggplot(smoothed_data, aes(x = as.numeric(bin), y = count_smooth, color = entry)) +
  geom_line(lwd=1.5) +
  geom_vline(xintercept=c(10),lty=2,show.legend=F) +
  labs(x="Return Day", y="Count",color = "Return") +
  theme_classic()+
  #scale_color_brewer(type="seq",palette = "Blues",direction = -1)+
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.position = "none")+
  scale_color_manual(values=return_palette[c(1,1,2,3,4)])



ggsave(Fig1b,dpi=500,filename = "Figures/Figure1_hist.png",units = "in",height = 3,width=4)  


#Experiment 1 Plots
load("output/Experiment_1_sumstats.Rdata")

FigS2<-summary_stats_all %>%
  group_by(G,pHOS,MU_entry_day_hatchery,origin,iter) %>% 
  filter(!is.na(meanK)) %>% 
  dplyr::summarize(meanK=mean(meanK)) %>% ungroup() %>% group_by(G,pHOS,MU_entry_day_hatchery,iter) %>% 
  dplyr::summarize(RRS=meanK[origin=="Hatchery"]/meanK[origin=="Wild"]) %>%
  group_by(G,pHOS,MU_entry_day_hatchery) %>%
  dplyr::summarize(sdRRS = sd(RRS),n=n(),RRS=mean(RRS),error = qnorm(0.975)*sdRRS/sqrt(n)) %>% 
  mutate(RRS = ifelse(is.na(RRS),0,RRS)) %>% 
  ggplot(aes(x=G,y=RRS,color=as.factor(MU_entry_day_hatchery-10)))+
  geom_ribbon(aes(fill=as.factor(MU_entry_day_hatchery),ymin=RRS-error,ymax=RRS+error,group=MU_entry_day_hatchery),alpha=0.4,show.legend = F)+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = c(0.25,0.5,0.75,1),lty = 2)+
  theme_classic()+
  labs(x="Generation",color= "Hatchery Mismatch\nFrom Local Run\nTiming (days)",y = "Relative Reproductive Success of Hatchery Fish")+
  scale_color_manual(values=return_palette)+
  scale_fill_manual(values=return_palette)+
  scale_y_continuous(breaks = c(0.25,0.5,0.75,1))

ggsave(FigS2,dpi=300,filename = "Figures/FigureS2.png",units = "in",height = 4,width=5)  


Fig1c<- summary_stats_all %>% filter(G == 25,!is.na(meanK)) %>% 
  group_by(pHOS,MU_entry_day_hatchery,iter) %>% 
  dplyr::summarize(RRS=meanK[origin=="Hatchery"]/meanK[origin=="Wild"]) %>%
  group_by(pHOS,MU_entry_day_hatchery) %>% 
  dplyr::summarize(meanRRS=mean(RRS),sdRRS = sd(RRS), error = qnorm(0.975)*sdRRS/sqrt(n())) %>% 
  ggplot(aes(x = MU_entry_day_hatchery - 10, y=meanRRS,fill=as.factor(MU_entry_day_hatchery-10)))+
  geom_col(color="black",show.legend = F)+
  geom_errorbar(aes(ymin = meanRRS - error, ymax = meanRRS + error),width =1.5)+
  theme_classic()+
  scale_color_manual(values=return_palette)+
  scale_fill_manual(values=return_palette)+
  labs(x="Hatchery Mismatch From\nLocal Run Timing (days)", y = "Relative Reproductive Success\nof Hatchery-Origin Fish")+
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.position = "none")

Figure1<-ggarrange(Fig1b,Fig1c,align = "h")

ggsave(Figure1,dpi=500,filename = "Figures/Figure1_hist.png",units = "in",height = 3.5,width=7)  


summary_stats_all %>% filter(G == 25,!is.na(meanK)) %>% 
  group_by(pHOS,MU_entry_day_hatchery,iter) %>% 
  dplyr::summarize(RRS=meanK[origin=="Hatchery"]/meanK[origin=="Wild"]) %>%
  group_by(pHOS,MU_entry_day_hatchery) %>% 
  dplyr::summarize(meanRRS=mean(RRS),sdRRS = sd(RRS), error = qnorm(0.975)*sdRRS/sqrt(n()),
                   Lower = meanRRS - error, Upper = meanRRS + error)

summary_stats_all %>% filter(G == 25,!is.na(meanK)) %>% 
  group_by(pHOS,MU_entry_day_hatchery,origin) %>% 
  dplyr::summarize(mean_entry = mean(mean_entry_day),sd_entry = sd(mean_entry_day),var = var(mean_entry_day))

#Experiment 2 Plots
load("output/Experiment_2-3_sumstats.Rdata")

Nc_plot_MU_hatchery<-summary_stats_all %>%
  filter(pHOS == 0.1, 
         origin=="Wild") %>%  
  group_by(G,MU_entry_day_hatchery,pHOS) %>% dplyr::summarize(sd=sd(Nc),Nc=mean(Nc),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=Nc,color=as.factor(MU_entry_day_hatchery)))+
  geom_ribbon(aes(fill=as.factor(MU_entry_day_hatchery),x=G,ymin=Nc-error,ymax=Nc+error),alpha=0.4,show.legend = F)+
  geom_point()+
  geom_vline(xintercept=25,lty=2)+
  theme_classic()+
  labs(x="Generation",color= "MU_entry_day_hatchery")+
  scale_color_brewer(type="seq",palette = "Blues")+
  scale_fill_brewer(type="seq",palette="Blues")+
  facet_grid(.~MU_entry_day_hatchery)
Nc_plot_MU_hatchery


Fig2<-summary_stats_all %>%
  filter(#pHOS == 0.1, 
    origin=="Wild") %>%  
  group_by(G,MU_entry_day_hatchery,pHOS) %>% dplyr::summarize(sd=sd(Nc),Nc=mean(Nc),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=Nc,color=as.factor(pHOS)))+
  geom_ribbon(aes(fill=as.factor(pHOS),x=G,ymin=Nc-error,ymax=Nc+error),alpha=0.4,show.legend = F)+
  geom_point(size = 0.8)+
  geom_vline(xintercept=25,lty=2)+
  theme_classic()+
  ggtitle("",subtitle = "Hatchery Mismatch From Local Run Timing (days)")+
  labs(x="Generation",color= "pHOS",y="Natural-Origin Population Size")+
  #scale_color_brewer(type="seq",palette = "Dark2")+
  #scale_fill_brewer(type="seq",palette="Dark2")+
  scale_color_manual(values=pHOS_palette)+
  scale_fill_manual(values = pHOS_palette)+
  facet_grid(.~as.factor(MU_entry_day_hatchery-10))+
  guides(colour = guide_legend(override.aes = list(size=2)))
Fig2
ggsave(Fig2,dpi=500,filename = "Figures/Figure2.png",units = "in",height = 4,width=7)  

summary_stats_all %>%
  filter(G == 50, 
    origin=="Wild") %>%  
  group_by(G,MU_entry_day_hatchery,pHOS) %>% dplyr::summarize(sd=sd(Nc),Nc=mean(Nc),error=qnorm(0.975)*sd/sqrt(n()),
                                                              lower = Nc - error, upper = Nc + error)


Figure_S3<-summary_stats_all %>% 
  filter(origin == "Wild",G==25) %>% 
  group_by(pHOS,MU_entry_day_hatchery) %>% dplyr::summarize(sd=sd(Nc),Nc=mean(Nc),error=qnorm(0.975)*Nc/sqrt(n())) %>%  
  ggplot(aes(x=pHOS,y=Nc,color=as.factor(pHOS)))+
  geom_line(color="gray")+
  geom_errorbar(aes(ymin=Nc-error,ymax=Nc+error),width=0.03)+
  geom_point()+
  facet_grid(.~MU_entry_day_hatchery)+
  theme_classic()+
  labs(y="Population Size (Nc)",y = "pHOS",color="pHOS")+
  scale_color_manual(values = pHOS_palette)
ggsave(Figure_S3,dpi=500,filename = "Figures/FigureS3.png",units = "in",height = 4,width=7)  

#Experiment 3:
load("output/Experiment_2-3_sumstats.Rdata")

Fig3<-summary_stats_all %>%
  filter(origin == "Wild") %>% 
  group_by(G,pHOS,MU_entry_day_hatchery) %>% dplyr::summarize(sd=sd(mean_pHOA),pHOA=mean(mean_pHOA),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=pHOA,color=as.factor(MU_entry_day_hatchery-10)))+
  geom_vline(xintercept = 25,lty = 2)+
  geom_ribbon(aes(fill=as.factor(MU_entry_day_hatchery-10),x=G,ymin=pHOA-error,ymax=pHOA+error,group=MU_entry_day_hatchery),alpha=0.4,show.legend = F)+
  geom_point(size = 1)+
  #geom_line()+
  theme_classic()+
  labs(x="Generation",color= "Hatchery Mismatch\nFrom Local Run\nTiming (days)",
       y = "Proportion of Hatchery-Origin\nAncestry in Natural-Origin Fish")+
  scale_color_manual(values = return_palette[c(2,3,4)])+
  scale_fill_manual(values = return_palette[c(2,3,4)])+
  facet_grid(.~pHOS)+
  ggtitle("",subtitle = "Proportion of Hatchery-Origin Spawners Added Each Generation")+
  guides(colour = guide_legend(override.aes = list(size=2)))
ggsave(Fig3,dpi=300,filename = "Figures/Figure3.png",units = "in",height = 4,width=7)

summary_stats_all %>%
  filter(origin == "Wild",G==26) %>% 
  group_by(G,pHOS,MU_entry_day_hatchery) %>% dplyr::summarize(sd=sd(mean_pHOA),
                                                              pHOA=mean(mean_pHOA),
                                                              error=qnorm(0.975)*sd/sqrt(n()),
                                                              lower = pHOA - error,
                                                              upper = pHOA + error)

summary_stats_all %>%
  filter(origin == "Wild",G==25) %>% 
  group_by(G,pHOS,MU_entry_day_hatchery) %>% dplyr::summarize(sd=sd(mean_pHOA),pHOA=mean(mean_pHOA),error=qnorm(0.975)*sd/sqrt(n()),meanK=mean(meanK,na.rm=T)) %>% 
  ggplot(aes(x=pHOS,y=pHOA,color=factor(MU_entry_day_hatchery)))+
  geom_errorbar(aes(ymin=pHOA - error,ymax = pHOA+error,width = 0.1))+
  geom_point(size=3)+
  facet_grid(.~MU_entry_day_hatchery)+
  theme_classic()+
  labs(x="Mean Entry Day of Hatchery Fish",color= "Mean Hatchery Entry")

#Experiment 4 Plots
load("output/Experiment_4_sumstats.Rdata")

FigS4<-summary_stats_all %>%
  filter(origin=="Wild") %>%  
  group_by(G,MU_entry_day_wild,pHOS) %>% dplyr::summarize(sd=sd(mean_entry_day),mean_entry_day=mean(mean_entry_day),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=mean_entry_day,color=as.factor(MU_entry_day_wild)))+
  geom_vline(xintercept=25,lty=2)+
  geom_hline(yintercept=13,lty=2)+
  geom_ribbon(inherit.aes=F,aes(fill=as.factor(MU_entry_day_wild),x=G,ymin=mean_entry_day-error,ymax=mean_entry_day+error),alpha=0.4,show.legend = F)+
  geom_point()+
  geom_line()+
  theme_classic()+
  facet_grid(.~pHOS)+
  labs(x="Generation",color= "Optimal Return Day", y = "Mean Return Day")+
  scale_color_brewer(type="seq",palette = "PRGn")+
  scale_fill_brewer(type="seq",palette="PRGn")+
  ggtitle("",subtitle = "Proportion of Hatchery-Origin Spawners Added Each Generation")
FigS4
ggsave(FigS4,dpi=300,filename = "Figures/FigureS4.png",units = "in",height = 4,width=7)

FigS4<-summary_stats_all %>%
  filter(origin=="Wild", G == 26) %>%  
  group_by(G,MU_entry_day_wild,pHOS) %>% dplyr::summarize(sd=sd(mean_entry_day),mean_entry_day=mean(mean_entry_day),error=qnorm(0.975)*sd/sqrt(n()),
                                                          lower = mean_entry_day - error, upper = mean_entry_day + error)
  

#Running average:
summary_stats_all %>%
  filter(origin == "Wild") %>%  
  group_by(G, MU_entry_day_wild, pHOS) %>%
  dplyr::summarize(
    sd = sd(mean_entry_day),
    mean_entry_day = mean(mean_entry_day),
    error = qnorm(0.975) * sd / sqrt(n())
  ) %>% 
  ggplot(aes(x = G, y = mean_entry_day)) +
  geom_point(aes(color = as.factor(MU_entry_day_wild))) +
  geom_smooth(aes(color = as.factor(MU_entry_day_wild)), method = "loess", se = FALSE) +  # Add running average
  theme_classic() +
  facet_grid(. ~ pHOS) +
  labs(x = "Generation", color = "Optimal Return Day", y = "Mean Return Day")

CV_portf <- summary_stats_all %>%
  filter(origin=="Wild", G == 25) %>%  
  group_by(MU_entry_day_wild,pHOS) %>% dplyr::summarize(sd=sd(mean_entry_day),mean_entry_day=mean(mean_entry_day),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ungroup() %>% group_by(pHOS) %>% 
  dplyr::summarise(CV = 100*sd(mean_entry_day)/mean(mean_entry_day))  
CV_portf

Fig4<-summary_stats_all %>%
  filter(origin=="Wild", G == 25) %>%  
  group_by(MU_entry_day_wild,pHOS) %>% dplyr::summarize(sd=sd(mean_entry_day),mean_entry_day=mean(mean_entry_day),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=MU_entry_day_wild,y=mean_entry_day,color=as.factor(pHOS)))+
  geom_hline(yintercept = 13, lty = 2)+
  geom_errorbar(aes(ymin=mean_entry_day-error,ymax=mean_entry_day+error),width=0)+
  geom_point()+
  geom_line()+
  facet_grid(.~pHOS)+
  geom_text(data=CV_portf,aes(label=paste("CV =",round(CV,1))),x = 10,y=20,inherit.aes = F)+
  labs(x="Optimal Return Day For Different Populations",
       color= "pHOS",
       y = "Mean Return Day After 25 Generations")+
  ggtitle("",subtitle = "Proportion of Hatchery-Origin Spawners Added Each Generation")+
  scale_color_manual(values=pHOS_palette)+
  scale_fill_manual(values = pHOS_palette)+
  theme_classic()+
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)#,
    #strip.background = element_rect(fill = "transparent", colour = NA)
  )
Fig4
ggsave(Fig4,dpi=500,filename = "Figures/Figure5.png",units = "in",height = 5,width=7,bg = "transparent")

summary_stats_all %>%
  filter(origin=="Wild", G == 25, MU_entry_day_wild == 20) %>%  
  group_by(MU_entry_day_wild,pHOS) %>% dplyr::summarize(sd=sd(mean_entry_day),mean_entry_day=mean(mean_entry_day),error=qnorm(0.975)*sd/sqrt(n()),
                                                        lower = mean_entry_day - error, upper = mean_entry_day + error)
  


###Density Dependence and number of mates supplemental:
load("output/Density_and_N_Mates_sumstats.Rdata")

Figure_S1<-summary_stats_all %>%
  filter(origin == "Wild") %>% 
  group_by(G,pHOS,n_mates_lambda,beta) %>% dplyr::summarize(sd=sd(Nc),Nc=mean(Nc),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=Nc,color=as.factor(beta)))+
  geom_vline(xintercept = 25, lty = 2)+
  geom_ribbon(aes(fill=as.factor(beta),x=G,ymin=Nc-error,ymax=Nc+error),alpha=0.4,show.legend = F)+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x="Generation",color= expression(beta), y = "Natural-Origin Population Size")+
  scale_color_brewer(type="seq",palette = "Purples")+
  scale_fill_brewer(type="seq",palette="Purples")+
  facet_grid(pHOS~n_mates_lambda,scales="free_y")+
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "pHOS", breaks = NULL, labels = NULL))+
  scale_x_continuous(sec.axis = sec_axis(~ . , name = expression(lambda["mates"]), breaks = NULL, labels = NULL))
Figure_S1
ggsave(Figure_S1,dpi=500,filename = "Figures/Figure_S1.png",units = "in",height = 5,width=7)




#Testing the effect of n_mates_lambda:
a<-summary_stats_all %>% 
  filter(origin == "Wild",G==10) %>% mutate(n_mates_lambda = ifelse(is.na(n_mates_lambda),Inf,n_mates_lambda)) %>% 
  group_by(n_mates_lambda) %>% dplyr::summarize(sd=sd(Nc),Nc=mean(Nc),error=qnorm(0.975)*Nc/sqrt(n())) %>%  
  ggplot(aes(x=n_mates_lambda,y=Nc))+
  geom_errorbar(aes(ymin=Nc-error,ymax=Nc+error),width=0.03)+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(y="Population Size (Nc)\nafter 10 generations",y = "Mean and variance in Number of Mates")
a

b<-summary_stats_all %>% mutate(n_mates_lambda = ifelse(is.na(n_mates_lambda),Inf,n_mates_lambda)) %>% 
  filter(origin=="Wild") %>%  
  group_by(G,n_mates_lambda) %>% dplyr::summarize(sd=sd(mean_entry_day),mean_entry_day=mean(mean_entry_day),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=mean_entry_day,color=as.factor(n_mates_lambda)))+
  geom_ribbon(aes(fill=as.factor(n_mates_lambda),ymin=mean_entry_day-error,ymax=mean_entry_day+error),alpha=0.4,show.legend = F)+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x="Generation",color= "n_mates_lambda")
b  

c<-summary_stats_all %>% mutate(n_mates_lambda = ifelse(is.na(n_mates_lambda),Inf,n_mates_lambda)) %>% 
  filter(origin=="Wild") %>%  
  group_by(G,n_mates_lambda) %>% dplyr::summarize(sd=sd(mean_RLS),mean_RLS=mean(mean_RLS),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=mean_RLS,color=as.factor(n_mates_lambda)))+
  geom_ribbon(aes(fill=as.factor(n_mates_lambda),ymin=mean_RLS-error,ymax=mean_RLS+error),alpha=0.4,show.legend = F)+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x="Generation",color= "n_mates_lambda")
c

d<-summary_stats_all %>% mutate(n_mates_lambda = ifelse(is.na(n_mates_lambda),Inf,n_mates_lambda)) %>% 
  filter(origin=="Wild") %>%  
  group_by(G,n_mates_lambda) %>% dplyr::summarize(sd=sd(Nc),mean_Nc=mean(Nc),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=mean_Nc,color=as.factor(n_mates_lambda)))+
  geom_ribbon(aes(fill=as.factor(n_mates_lambda),ymin=mean_Nc-error,ymax=mean_Nc+error),alpha=0.4,show.legend = F)+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x="Generation",color= "n_mates_lambda")
d

e<-summary_stats_all %>% mutate(n_mates_lambda = ifelse(is.na(n_mates_lambda),Inf,n_mates_lambda)) %>% 
  filter(origin=="Wild") %>%  
  group_by(G,n_mates_lambda) %>% dplyr::summarize(sd=sd(Ne,na.rm=T),mean_Ne=mean(Ne,na.rm=T),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=mean_Ne,color=as.factor(n_mates_lambda)))+
  geom_ribbon(aes(fill=as.factor(n_mates_lambda),ymin=mean_Ne-error,ymax=mean_Ne+error),alpha=0.4,show.legend = F)+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x="Generation",color= "n_mates_lambda")
e  
f<-summary_stats_all %>% mutate(n_mates_lambda = ifelse(is.na(n_mates_lambda),Inf,n_mates_lambda)) %>% 
  filter(origin=="Wild") %>%  
  group_by(G,n_mates_lambda) %>% dplyr::summarize(sd=sd(Ne_N,na.rm=T),mean_Ne=mean(Ne_N,na.rm=T),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=mean_Ne,color=as.factor(n_mates_lambda)))+
  geom_ribbon(aes(fill=as.factor(n_mates_lambda),ymin=mean_Ne-error,ymax=mean_Ne+error),alpha=0.4,show.legend = F)+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x="Generation",color= "n_mates_lambda",y = "Ne/N")
f  

ggarrange(a,b,c,d,e,f,common.legend = T)



##Single iteration example distribution of RS and number of mates:

output_example<-read.csv("output/portfolios/gen_data_2002.csv")

output_example<-output_example %>% filter(G==1)
output_example %>% ggplot(aes(x=RS_obs))+geom_histogram(binwidth = 1)


dams<-unique(output_example$dam)
sires<-unique(output_example$sire)

dams<-data.frame(id=dams,parent="dam",n_mates=NA)
sires<-data.frame(id=sires,parent="sire",n_mates=NA)

for(i in 1:nrow(dams)){
  dams$n_mates[i]<-length(unique(output_example$sire[which(output_example$dam==dams$id[i])]))
}

for(i in 1:nrow(sires)){
  sires$n_mates[i]<-length(unique(output_example$dam[which(output_example$sire==sires$id[i])]))
}
rbind(sires,dams) %>% ggplot(aes(x=n_mates))+geom_histogram(binwidth = 1)+facet_grid(.~parent)









