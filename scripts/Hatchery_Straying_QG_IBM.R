#title: "Quantitative genetic model to examine the effects of hatchery straying on wild population recruitment and viability"
#authors: "Samuel May..."
#date: "09/08/2022"




#### Load Data and Packages ####

library(dplyr)
library(TruncatedNormal)
library(rockchalk)
library(purrr)
library(foreach)
library(doParallel)


seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

#### Description of Function Arguments: ####

#iteration = number of iterations to run
#rho = phenotypic correlation between entry day and RLS
#variance_entry_day = phenotypic variance in entry day
#variance_RLS = phenotypic variance in RLS
#assortative = should mating occur randomly or assortatively? (T or F) 
#PopSize = allow the population size to vary or remain constant at Nc_initial? (variable or constant)... 
#G = number of generations to run for
#Nc_initial = number of individuals in the starting population
#OMEGA = shape of the fitness surface - approx 1-4 phenotypic standard deviations
#var_env_entry = Interannual variability in optimal entry day
#var_env_RLS = Interannual variability in optimal lifespan
#PATH = the directory to which output files should be written

#New arguments:
#pHOS: proportion of hatchery-origin spawners to be added to the wild population
#MU_entry_day_hatchery: average day of hatchery entry
#MU_RLS_hatchery: average hatchery lifespan
#season_length: number of days in spawning season
#F0_hatchery : In the first generation, should pHOS be additive or subtractive:
#in other words, among simulations do we want to compare F0 with the same sized wild population 
#or the same size total population (hatchery + wild)
#F1_hatchery : pHOS is constant, do we want to add or subtract from wild pop - or constant (same as F0 number)
#beta: strength of density dependence
#hatchery_gens : how many generations to run the hatchery for?
#n_mates_lambda : mean and variance of a Poisson distribution, from which the maximum number of possible mates for each putative female parent is drawn.


##All files currently in directory:
file_names = list.files(path = "output/Density_Dependence/",pattern="*.csv",full.names = T,include.dirs = F)


#### Model Initialization ####

Assortative_Mating_IBM<-function(iteration = 1,
                                 rho = -0.3,
                                 MU_optimum_entry = 10,
                                 MU_entry_day = 10,
                                 MU_RLS = 7,
                                 variance_entry_day = 20,
                                 variance_RLS = 10,
                                 PopSize = "variable",
                                 G = 50,
                                 OMEGA = 2,
                                 var_env_entry = 20,
                                 var_env_RLS = 1,
                                 Nc_initial = 500,
                                 PATH = "output/n_mates_test2/",
                                 pHOS = 0.2,
                                 MU_entry_day_hatchery = 20,
                                 MU_RLS_hatchery = 7,
                                 season_length = 45,
                                 F0_hatchery = "constant_wild", #constant_wild or "constant_total"
                                 F1_hatchery = "add",#add, subtract, or constant
                                 beta = -1e-5,
                                 hatchery_gens = 25,
                                 n_mates_lambda = 2){ 
  
  
  ###HPC keeps crashing for some reason, so let's just make sure when we are running something, that we don't put
  #computational effort into re-running things that already have output files:
  file_names = list.files(path = PATH,pattern="*.csv",full.names = T,include.dirs = F)
  output_file_name<-paste(PATH,"gen_data_",n_mates_lambda,"_",beta,"_",MU_entry_day,"_",MU_entry_day_hatchery,"_",pHOS,"_",
                          iteration,".csv",sep="")
  if(output_file_name%in%file_names){break}
  
  #How many hatchery fish should be initialized?
  if(F0_hatchery == "constant_wild"){
    Nc_hatchery <- round((pHOS*Nc_initial)/(1-pHOS))}
  if(F0_hatchery == "constant_total"){
    Nc_hatchery <- Nc_initial*pHOS
    Nc_initial <- Nc_initial - Nc_hatchery
  }
  
  Nc_total = Nc_initial+Nc_hatchery
  
  pHOS_initial <- pHOS #pHOS may change later - we want to use this in our outputs
  #Make Sigma (The phenotpyic VCV matrix, or "P matrix"):
  #Sigma = Phenotypic Matrix = Additive Genetic VCV matrix + Residual (environmental) VCV matrix
  
  Sigma<-rockchalk::lazyCov(Rho = rho, Sd = c(sqrt(variance_RLS),sqrt(variance_entry_day)))
  rownames(Sigma)<-c("RLS","entry_day")
  colnames(Sigma)<-c("RLS","entry_day")
  
  Corr<-cov2cor(Sigma)
  
  #Initiate the starting population (F0) from Sigma and population trait means. 
  F0<-data.frame(id=1:Nc_total,
                 entry_day = NA,
                 RLS = NA,
                 origin=c(rep("Wild",Nc_initial),rep("Hatchery",Nc_hatchery)),
                 sex=rep(c("m","f"),length.out=Nc_total),
                 year=rep(2022,Nc_total)) #year is irrelevant and redundant with Generation... but could be used if desired.
  u = c(MU_RLS,MU_entry_day) 
  u_hatchery = c(MU_RLS_hatchery,MU_entry_day_hatchery)
  
  #Draw entry days for Nc_initial individuals
  F0$entry_day[F0$origin=="Wild"]<-TruncatedNormal::rtnorm(n=Nc_initial,mu=u[2],sd=sqrt(Sigma[2,2]),lb=0,ub=season_length) #30 day spawning season
  
  #Draw entry days for Nc_hatchery individuals
  if(pHOS>0){
    F0$entry_day[F0$origin=="Hatchery"]<-TruncatedNormal::rtnorm(n=Nc_hatchery,mu=u_hatchery[2],sd=sqrt(Sigma[2,2]),lb=0,ub=season_length)} #30 day spawning season
  
  #Draw RLS values for those individuals, conditional upon entry day such that the traits will be correlated
  for(i in 1:nrow(F0)){
    if(F0$origin[i]=="Wild"){
      cond_mean = u[1] + Corr[1,2]*sqrt(Sigma[1,1])/sqrt(Sigma[2,2]) * (F0[i,"entry_day"] - u[2])
      cond_sd = sqrt(Sigma[1,1]) * sqrt(1 - Corr[1,2])}
    if(F0$origin[i]=="Hatchery"){
      cond_mean = u_hatchery[1] + Corr[1,2]*sqrt(Sigma[1,1])/sqrt(Sigma[2,2]) * (F0[i,"entry_day"] - u_hatchery[2])
      cond_sd = sqrt(Sigma[1,1]) * sqrt(1 - Corr[1,2])}
    
    
    #RLS cannot be greater than the length of the spawning season (30) - entry day
    F0$RLS[i]<-TruncatedNormal::rtnorm(1,mu=cond_mean,sd=cond_sd,lb=0,ub=season_length-F0[i,"entry_day"])
  }
  
  #Round these trait values to the nearest day
  F0$entry_day<-round(F0$entry_day)
  F0$RLS<-round(F0$RLS)
  
  
  #Initialize the data to keep track of for each individual within each generation (gen_data)
  #The gen_data object will be be the eventual output of the model
  #Add F0 information to this data frame
  gen_data <- matrix(nrow=Nc_total,ncol=9)
  colnames(gen_data) <- c("ID", "G", "sire", "dam", "sex", "entry_day", "RLS","RS_exp","RS_obs") 
  gen_data <- as_tibble(gen_data)
  gen_data$ID<-1:nrow(F0)
  gen_data$G<-0
  gen_data$sex<-F0$sex
  gen_data$entry_day<-F0$entry_day
  gen_data$RLS<-F0$RLS
  gen_data$year<-F0$year
  gen_data$origin<-F0$origin
  gen_data<-gen_data %>% mutate(pHOA = ifelse(origin=="Wild",0,1)) #prop. hatchery origin ancestry
  
  
  for (g in 0:G){ #Begin looping through generations
    
    #Set parents of current generation
    parents<-gen_data%>%filter(G==g)
    
    #### Fitness Estimation Module ####
    
    #Fitness surface from Lande and Arnold, determined by trait optima (theta values) and omega
    #Here, OMEGA is the number of standard deviations for the multivariate normal distribution
    #Which determines the strength of selection. We use mean theta values equal to the initial
    #trait mean in the wild population (MU_entry_day, MU_RLS) 
    
    #Draw a new theta value for each trait, for each generation. 
    #environmental variance (var_env) determines how much theta can vary from year to year
    THETA_entry_day <- TruncatedNormal::rtnorm(n=1, mu = MU_optimum_entry, lb=1,ub=season_length, sd = sqrt(var_env_entry))
    THETA_RLS <- TruncatedNormal::rtnorm(n=1, mu = MU_RLS, lb=1, ub=season_length-THETA_entry_day, sd = sqrt(var_env_RLS))
    
    
    #Fitness estimation of a two-trait surface, from Lande and Arnold
    parents$RS_exp<-exp((-1/(2*(1-(rho^2))*OMEGA))*(((parents$entry_day - (THETA_entry_day+0))/sqrt(variance_entry_day))^2 + 
                                                      ((parents$RLS - (THETA_RLS+0))/sqrt(variance_RLS))^2 - 
                                                      2*rho*((parents$entry_day - (THETA_entry_day+0))/sqrt(variance_entry_day))*((parents$RLS - 
                                                                                                                                     (THETA_RLS+0))/sqrt(variance_RLS))))
    
    
    #for each individual, make a list of how many days they were there
    parents<-parents %>% mutate(days_present = seq2(entry_day,entry_day+RLS))
    
    #count how many individuals were in the stream on each day:
    total_vs_day<-data.frame(day = seq(0,season_length),Nc_day=NA)
    for(i in 0:nrow(total_vs_day)){
      total_vs_day$Nc_day[i] <- sum(unlist(parents$days_present)%in%total_vs_day$day[i])}
    
    #How many fish were in the stream on the day a fish entered:
    parents$Nc_day<-total_vs_day$Nc_day[match(parents$entry_day,total_vs_day$day)]
    #Use this Nc_day in a Ricker function:
    parents$RS_exp <- parents$RS_exp * exp(beta * parents$Nc_day)
    
    #Draw Expected RS values from a Poisson distribution, using the 0-1 scale generated above as probabilities
    parents$RS_exp<-qpois(parents$RS_exp,lambda=2)
    
    #Size of the next generation is the sum of RS_exp values of all females in the populations,
    #or we can hold the population size constant (Nc_initial is the carrying capacity)
    if(PopSize=="variable"){
      Nc<-round(parents%>%filter(sex=='f')%>%
                  #group_by(spawning_location)%>%
                  dplyr::summarize(Nc=sum(RS_exp))%>%
                  dplyr::select(Nc)%>%as.matrix%>%c())
      
      if(F1_hatchery=="subtract"){
        Nc_hatchery<-round(pHOS*Nc)
        Nc <- Nc-Nc_hatchery}
    }
    
    if(PopSize=="constant"){
      Nc<-Nc_initial
    }
    
    if(Nc<50){break} #If the population crashes below 50, assume extinction.
    if(Nc>4000){Nc<-4000}
    
    #### Reproduction Module ####
    # Here we draw parents for the number of offspring in the next generation
    
    #days_df provides one row for each day that each individual was available to mate  
    days_df<-as.data.frame(matrix(nrow=0,ncol=4)) 
    colnames(days_df)<-c('day','ID','sex')
    
    for (i in 1:nrow(parents)){
      ind_i<-parents$ID[i]	
      days<-seq(parents$entry_day[i],parents$entry_day[i]+parents$RLS[i])
      sex<-parents$sex[i]
      temp_df<-data.frame(day=days,ID=rep(ind_i,length(days)), sex=rep(sex,length(days)))
      days_df<-rbind(days_df,temp_df)
    }
    
    
    #The following is a slight code modification from the original MS that is optimized so it runs faster:
    #instead of making a matrix of Nc x Nc, we just skip that and make a data frame of possible parent-pairs. 
    
    # subset data for only male individuals
    male_df <- days_df %>% filter(sex == 'm')
    
    # subset data for only female individuals
    female_df <- days_df %>% filter(sex == 'f')
    
    # create pairwise data frame of overlapping males and females
    possible_mates <- male_df %>% 
      inner_join(female_df, by = 'day') %>% 
      select(ID_1 = ID.x, ID_2 = ID.y, sex_1 = sex.x, sex_2 = sex.y) %>% 
      distinct()
    
    
    
    
    ######################################
    #To better replicate wild mating dynamics, here we limit the number of mates per female
    #We assume a Poisson distribution with a mean and variance equal to ~1 (input parameter n_mates_lambda)
    #We will test a range of values to see how this affects model behavior
    
    possible_mates_filtered <- possible_mates
    
    if(is.na(n_mates_lambda)==F){
      #Get the female parents:
      female_parents<-parents %>% filter(sex=="f")
      
      #Assign each female parent a maximum number of mates drawn from possible mates in days_matrix
      female_parents$n_mates<-rpois(n = nrow(female_parents),lambda=n_mates_lambda)
      
      #make a new data frame that's a filtered version of possible_mates, where each female can't match with a male more than n_mates:
      
      possible_mates_filtered <- possible_mates %>%
        group_split(ID_2) %>%
        map_dfr(~ {
          filtered_df <- .
          n_mates <- female_parents$n_mates[female_parents$ID == filtered_df$ID_2[1]]
          if (nrow(filtered_df) > n_mates) {
            filtered_df <- filtered_df %>% sample_n(n_mates,)
          }
          return(filtered_df)
        })
    }
    #This doesn't use purrr, but is slower:
    #possible_mates_filtered <- lapply(split(possible_mates, possible_mates$ID_2), function(df) {
    #  filtered_df <- df
    #  n_mates <- female_parents$n_mates[female_parents$ID == filtered_df$ID_2[1]]
    #  if (nrow(filtered_df) > n_mates) {
    #    filtered_df <- filtered_df[sample(nrow(filtered_df), n_mates), ]
    #  }
    #  return(filtered_df)
    #})
    #
    #possible_mates_filtered <- do.call(rbind, possible_mates_filtered)
    #}
    ######################################
    
    #In original script, we used all parents - here we only use possible_mates_filtered
    #Assign Parent Weights (W[i,j]) with overlap weight * fitness weight of dad * fitness weight of mum:
    
    #relative reproductive success values are not sex-specific.
    parents<-parents%>%mutate(Nc=Nc)
    parents<-parents%>%mutate(meanKexp=mean(RS_exp),fitness_weight=RS_exp/meanKexp)
    
    
    #Add fitness weights (Relative reproductive success values)
    possible_mates_filtered$dam_fitness<-parents$fitness_weight[match(possible_mates_filtered$ID_2,parents$ID)]
    possible_mates_filtered$sire_fitness<-parents$fitness_weight[match(possible_mates_filtered$ID_1,parents$ID)]
    
    
    possible_mates_filtered<-possible_mates_filtered %>% mutate(weight=dam_fitness*sire_fitness)
    
    #Draw parents using weights
    
    #Initialize new data frame for offspring
    offspring<-data.frame(ID=seq(max(parents$ID)+1,max(parents$ID)+sum(Nc)),
                          dam=NA,
                          sire=NA,
                          pHOA=NA,
                          gen=g+1)
    
    #Parents:
    offspring[,c(3,2)]<-possible_mates_filtered[sample(x=1:nrow(possible_mates_filtered),
                                                       size=Nc,replace=T,
                                                       prob=possible_mates_filtered$weight),1:2]
    
    
    #### Inheritance Module ####
    #Here we assign traits to offspring
    
    #Get the trait values of their parents
    offspring$dam_return<-gen_data$entry_day[match(offspring$dam,gen_data$ID)]
    offspring$sire_return<-gen_data$entry_day[match(offspring$sire,gen_data$ID)]
    offspring$midparent_return<-rowMeans(offspring%>%dplyr::select(dam_return,sire_return))
    
    offspring$dam_RLS<-gen_data$RLS[match(offspring$dam,gen_data$ID)]
    offspring$sire_RLS<-gen_data$RLS[match(offspring$sire,gen_data$ID)]
    offspring$midparent_RLS<-rowMeans(offspring%>%dplyr::select(dam_RLS,sire_RLS))
    
    offspring$dam_pHOA<-gen_data$pHOA[match(offspring$dam,gen_data$ID)]
    offspring$sire_pHOA<-gen_data$pHOA[match(offspring$sire,gen_data$ID)]
    offspring<-offspring %>% mutate(pHOA=(dam_pHOA+sire_pHOA)/2)
    
    #Inheritance of correlated traits
    for(i in 1:nrow(offspring)){
      
      #Means equal midparent values
      u = c(offspring$midparent_RLS[i],offspring$midparent_return[i]) #midparent RLS , midparents entry_day
      
      #Draw entry day from truncated normal from 0 to 30
      offspring[i,"entry_day"]<-TruncatedNormal::rtnorm(n=1,mu=u[2],sd=sqrt(Sigma[2,2]),lb=0,ub=season_length) #30 day spawning season
      
      #RLS is correlated with entry day, so drawn from conditional mean and sd.
      cond_mean = u[1] + (Corr[1,2] * sqrt(Sigma[1,1])/sqrt(Sigma[2,2]) * (offspring[i,"entry_day"] - u[2]))
      cond_sd = sqrt(Sigma[1,1]) * sqrt(1 - Corr[1,2])
      
      #Individuals can't live past the last day of the season.
      offspring[i,"RLS"]<-TruncatedNormal::rtnorm(n=1,mu=cond_mean,sd=cond_sd,lb=0,ub=season_length-offspring[i,"entry_day"]) 
    }
    
    offspring$RLS<-round(offspring$RLS)
    offspring$entry_day<-round(offspring$entry_day)
    
    #Assuming a sex ratio of 1:1 and randomly draw sexes. 
    #Thus traits are completely independent of sex
    offspring$sex<-sample(c('m','f'),size = nrow(offspring),replace = T) 
    
    #Assign parents RS_obs
    parents$RS_obs<-0
    sires_RS_obs<-offspring%>%group_by(sire)%>%dplyr::summarize(RS_obs=n(),.groups="drop")%>%dplyr::rename(parent=sire) %>% as.data.frame()
    dams_RS_obs<-offspring%>%group_by(dam)%>%dplyr::summarize(RS_obs=n(),.groups="drop")%>%dplyr::rename(parent=dam)%>%as.data.frame()
    parents_RS_obs<-rbind(sires_RS_obs,dams_RS_obs)
    
    parents$RS_obs[match(parents_RS_obs$parent,parents$ID)]<-parents_RS_obs$RS_obs
    
    
    
    ############......Finish parents in gen_data################
    
    gen_data$RS_exp[which(gen_data$G==g)]<-parents$RS_exp
    gen_data$RS_obs[which(gen_data$G==g)]<-parents$RS_obs
    
    ############......Add offspring as next generation in gen_data################
    
    
    offspring$RS_obs<-NA
    offspring$RS_exp<-NA
    offspring$year<-max(gen_data$year)+4 #arbitrary generation length of 4 years.
    offspring$origin<-"Wild"
    
    gen_data_g<-offspring%>%dplyr::select(ID,gen,sire,dam,sex,entry_day,RLS,RS_exp,RS_obs,year,origin,pHOA)
    colnames(gen_data_g)<-colnames(gen_data)
    
    gen_data<-rbind(gen_data,gen_data_g)
    
    ############...Add hatchery fish for next generation in gen_data#########
    
    #If we want to "turn off the tap" and remove hatchery fish after so long, after how many generations?
    if(g>hatchery_gens){
      pHOS <- 0}
    
    if(pHOS>0){
      if(F1_hatchery=="add"){
        Nc_hatchery<-round((pHOS*nrow(offspring))/(1-pHOS))}
      
      
      if(F1_hatchery=="constant"){
        Nc_hatchery<-Nc_hatchery}
      
      hatchery_next_gen<-data.frame(ID=(max(gen_data$ID)+1):(max(gen_data$ID)+Nc_hatchery),
                                    G = max(gen_data$G),
                                    sire = NA, dam = NA, sex = rep(c("m","f"),length.out=Nc_hatchery),
                                    entry_day = NA, RLS = NA,
                                    RS_exp = NA, RS_obs=NA,
                                    year=max(gen_data$year), #year is irrelevant and redundant with Generation... but could be used if desired.
                                    origin=rep("Hatchery",Nc_hatchery),
                                    pHOA=1)
      u_hatchery = c(MU_RLS_hatchery,MU_entry_day_hatchery)
      
      for(i in 1:nrow(hatchery_next_gen)){  
        hatchery_next_gen$entry_day<-TruncatedNormal::rtnorm(n=Nc_hatchery,mu=u_hatchery[2],sd=sqrt(Sigma[2,2]),lb=0,ub=season_length)
        
        cond_mean = u_hatchery[1] + Corr[1,2]*sqrt(Sigma[1,1])/sqrt(Sigma[2,2]) * (hatchery_next_gen[i,"entry_day"] - u_hatchery[2])
        cond_sd = sqrt(Sigma[1,1]) * sqrt(1 - Corr[1,2])
        
        
        # RLS cannot be greater than the length of the spawning season (45) - entry day
        hatchery_next_gen$RLS[i]<-TruncatedNormal::rtnorm(1,mu=cond_mean,sd=cond_sd,lb=0,ub=season_length-hatchery_next_gen[i,"entry_day"])
      }
      
      #Round these trait values to the nearest day
      hatchery_next_gen$entry_day<-round(hatchery_next_gen$entry_day)
      hatchery_next_gen$RLS<-round(hatchery_next_gen$RLS)
      
      gen_data<-rbind(gen_data,hatchery_next_gen)
    }
    
  } #end G Loop
  
  #Add input parameters to gen data
  gen_data$rho<-rho
  gen_data$beta<-beta
  gen_data$var_entry_day<-variance_entry_day
  gen_data$var_RLS<-variance_RLS
  gen_data$PopSize<-PopSize
  gen_data$iter<-iteration
  gen_data$var_env_entry = var_env_entry
  gen_data$var_env_RLS = var_env_RLS
  gen_data$OMEGA<- OMEGA
  gen_data$pHOS<-pHOS_initial
  gen_data$Nc_initial<-Nc_initial
  gen_data$MU_entry_day_wild<-MU_entry_day
  gen_data$MU_entry_day_hatchery<-MU_entry_day_hatchery
  gen_data$MU_RLS_wild<-MU_RLS
  gen_data$MU_RLS_hatchery<-MU_RLS_hatchery
  gen_data$F1_hatchery<-F1_hatchery
  gen_data$F0_hatchery <- F0_hatchery
  gen_data$n_mates_lambda<-n_mates_lambda
  
  
  
  #At the end of each iteration, export gen_data to a .csv
  #The title of this .csv contains all unique input parameters,
  #so that files can easily read in and concatenated with rbind()
  
  #write.csv(gen_data,file=paste(PATH,pHOS,MU_entry_day_hatchery,OMEGA,rho,variance_entry_day,variance_RLS,
  write.csv(gen_data,file=paste(PATH,"gen_data_",n_mates_lambda,"_",beta,"_",MU_entry_day,"_",MU_entry_day_hatchery,"_",pHOS_initial,"_",
                                iteration,".csv",sep=""))
  
  
} #end function




#Register a cluster of cores to run the function on in parallel:
clust<-makeCluster(8)
registerDoParallel(cl = clust, cores = 8)

#registerDoParallel(20)

#for each input parameter to vary, define a new parameter with desired values
#Here the new parameter is just the input parameter with a "z" on the end.

#Experiment 1:
file_names = list.files(path = "output/Experiment_1",pattern="*.csv",full.names = T,include.dirs = F)
foreach(iterz = 1:1000, .packages = c("dplyr","rockchalk","TruncatedNormal","doParallel","purrr"),
        .export = "file_names") %:% 
  foreach(MU_entry_day_hatcheryz = c(10,15,20,25)) %dopar% {
    Assortative_Mating_IBM(PATH = "output/Experiment_1/",
                           iteration = iterz,
                           pHOS = 0.2,
                           MU_entry_day_hatchery = MU_entry_day_hatcheryz)}
#stopImplicitCluster()

#Experiment 2-3:
file_names = list.files(path = "output/Experiment_2-3",pattern="*.csv",full.names = T,include.dirs = F)
foreach(iterz = 101:1000, .packages = c("dplyr","rockchalk","TruncatedNormal","purrr"),
        .export = "file_names") %:%
  foreach(pHOSz = c(0, 0.01,0.1,0.3)) %:%
  foreach(MU_entry_day_hatcheryz = c(15,20,25)) %dopar% {
    Assortative_Mating_IBM(PATH = "output/Experiment_2-3/",
                           iteration = iterz,
                           pHOS = pHOSz,
                           MU_entry_day_hatchery = MU_entry_day_hatcheryz)}

#Experiment 4
file_names = list.files(path = "output/Experiment_4",pattern="*.csv",full.names = T,include.dirs = F)
foreach(iterz = 1:1000, .packages = c("dplyr","rockchalk","TruncatedNormal","purrr"),
        .export = "file_names") %:%
  foreach(pHOSz = c(0,0.1,0.2,0.3)) %:%
  foreach(MU_entry_dayz = c(5,10,13,15,17,20)) %:%
  foreach(MU_entry_day_hatcheryz = 13) %dopar% {
    Assortative_Mating_IBM(PATH = "output/Experiment_4/",
                           iteration = iterz,
                           pHOS = pHOSz,
                           MU_entry_day = MU_entry_dayz,
                           MU_optimum_entry = MU_entry_dayz,
                           MU_entry_day_hatchery = MU_entry_day_hatcheryz)}


#Testing density dependence
foreach(iterz = 1:100, .packages = c("dplyr","rockchalk","TruncatedNormal","purrr")) %:%
  foreach(betaz = c(-0.001,-0.0001,-0.00001)) %:% 
  foreach(n_matesz = c(1,2,3,NA)) %:%
  foreach(pHOSz = c(0,0.3)) %dopar% {
    Assortative_Mating_IBM(G = 50,hatchery_gens = 25,
                           PATH = "output/Density_Dependence/",
                           iteration = iterz,
                           n_mates_lambda = n_matesz,
                           beta = betaz,
                           pHOS = pHOSz)}



stopImplicitCluster() #closes unused processor connections after running foreach


