library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(stringr)

### work directory
sdmoutput <- ""

setwd(sdmoutput)

# Species
plant_keys <- read.csv("../data_outputs/plant_keys.csv")
# plant_keys <- read.csv("C:/Users/claud/Desktop/UCADOCTORADO/REPOSITORIOS/chapter_3/data_outputs/plant_keys.csv")
plant_keys <- plant_keys[,"plant_sp"]

myRespName <- stringr::str_replace_all(plant_keys, " ", "") #eliminamos espacio para evitar problemas al guardar archvos
myRespName <- stringr::str_sub(myRespName, start=1,end=12) 

comb <- 12 # numero de combinaciones probadas

runROC <- NA
# Prepare our ROC evaluation
for(i in 1:length(plant_keys)) {  
  for (j in 1:comb){
   if(file.exists(paste(myRespName[i],"/",myRespName[i],"_mod_", j,"_evalDF_ROC.csv",sep = ""))){
     eval_ROC <- read.csv(file=paste(myRespName[i],"/",myRespName[i],"_mod_", j,"_evalDF_ROC.csv",sep = "")) 
     names(eval_ROC)[1] <- 'model'
     eval_ROC$species <- myRespName[i]
     eval_ROC$comb <- j
     eval_ROC <- pivot_longer(eval_ROC,cols =starts_with("RUN"),names_to="run",values_to ="roc")
     
     runROC <- rbind(runROC,eval_ROC)
   }

  }
}
  
runROC <- runROC[2:nrow(runROC),]
  
ggplot(runROC, aes(x=model, y=roc)) +
    geom_boxplot() +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    ggtitle("") +
    xlab("") +
    facet_wrap(~species)

runTSS <- NA
for(i in 1:length(plant_keys)) {  
  for (j in 1:comb){
    if(file.exists(paste(myRespName[i],"/",myRespName[i],"_mod_", j,"_evalDF_TSS.csv",sep = ""))){
      eval_TSS <- read.csv(file=paste(myRespName[i],"/",myRespName[i],"_mod_", j,"_evalDF_TSS.csv",sep = "")) 
      names(eval_TSS)[1] <- 'model'
      eval_TSS$species <- myRespName[i]
      eval_TSS$comb <- j
      eval_TSS <- pivot_longer(eval_TSS,cols =starts_with("RUN"),names_to="run",values_to ="TSS")
      
      runTSS <- rbind(runTSS,eval_TSS)
    }
    
  }
}

runTSS <- runTSS[2:nrow(runTSS),]

ggplot(runTSS, aes(x=as.factor(comb), y=TSS,color=model)) +
  geom_boxplot() +
  geom_jitter(size=0.4, alpha=0.9) +
  ggtitle("") +
  xlab("") +
  facet_wrap(~species)

eval_global <- merge(runROC,runTSS)

eval_global_mean <- eval_global %>%
  group_by(model,species) %>%
  summarise(roc=mean(roc,na.rm=T),tss=mean(TSS,na.rm=T))

# Plot
library(ggpubr)
ggscatterhist(eval_global_mean,x="roc",y="tss", color="model",
                        size = 3, alpha = 0.6,
                        palette = c("#00AFBB", "#E7B800", "#FC4E07","light green", "yellow", "red"),
                        margin.plot = "boxplot",
                        margin.params = list(fill = "model", size = 0.2),
                        label="species",font.label = c(10, "plain"),
                        ggtheme=theme_bw(base_size = 14),
                        xlab="ROC",ylab="TSS")





# Prepare our variable importance evaluation

imp_all <- NA

for(i in 1:length(plant_keys)) {  
  for (j in 1:comb){
  if(file.exists(paste(myRespName[i],"/",myRespName[i],"_mod_", j,"_variables_importance.csv",sep = ""))){
  
      var_imp <-read.csv(file=paste(myRespName[i],"/",myRespName[i],"_mod_", j,"_variables_importance.csv",sep = ""),header=FALSE) 
      
      seleccion <- seq(1,nrow(var_imp),15)
      
      for (k in seleccion){
        run <- var_imp[k,3] %>%
          str_sub(2,100) %>% #eliminamos espacio inicial 
          str_split_fixed(pattern="_",n=4)
        ini =k+2
        fin =k+13
        imp <- var_imp[ini:fin,1] 
        imp <- gsub("\\s+"," ",imp) # eliminamos espacios extra)
        imp <- as.data.frame(str_split_fixed(imp, pattern=" ",n=4))
        imp <- cbind(imp,run) %>%
          pivot_longer(cols =c("V2","V3","V4"),names_to="rand",values_to ="imp")
        names(imp) <- c("variable","species","PA","RUN","model","RUN_IMP","imp")
        imp$comb <- j
        imp_all <- rbind(imp_all,imp)
      }
  }
  }
  print(paste("species done: ", myRespName[i]))
  
}
  
imp_all <- imp_all[2:nrow(imp_all),]

imp_all$imp <- as.numeric(imp_all$imp)

# Save file
write.csv(imp_all, file=paste("variables_importance_stand.csv",sep = "")) 


# By RUN iteration
imp_all_mean <- imp_all %>%
  group_by(RUN,species,variable) %>%
  summarise(imp=mean(imp,na.rm=T))

  ggplot(imp_all_mean, aes(x=reorder(variable,-imp,FUN=median), y=imp, fill=variable)) +
    geom_boxplot() +
    ggtitle("") +
    xlab("")+
    facet_wrap(~ species)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  # By RUN of importance
  imp_all_mean <- imp_all %>%
    group_by(RUN_IMP,species,variable) %>%
    summarise(imp=mean(imp,na.rm=T))
  
  ggplot(imp_all_mean, aes(x=reorder(variable,-imp,FUN=median), y=imp, fill=variable)) +
    geom_boxplot() +
    ggtitle("") +
    xlab("")+
    facet_wrap(~ species)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
# By model
  imp_all_mean <- imp_all %>%
    group_by(model,species,variable) %>%
    summarise(imp=mean(imp,na.rm=T))
  
  ggplot(imp_all_mean, aes(x=reorder(variable,-imp,FUN=median), y=imp, fill=variable)) +
    geom_boxplot() +
    ggtitle("") +
    xlab("")+
    facet_wrap(~ species)+
    theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


  
  ggplot(imp_all_mean,aes(x=reorder(variable,-imp,FUN=median,), y=imp,fill=variable))+
  
    geom_violin(width=1.4) +
    
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    
    theme_bw(base_size = 16) + theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    
    xlab("") +   ylab("Variable importance (%)") +
      facet_wrap(~ species)





