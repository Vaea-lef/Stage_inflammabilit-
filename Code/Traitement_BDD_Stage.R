
#### TRAITEMENT DONNEE STAGE INFLAMMABILITE ####

############# Traitement des données feuille ###########################


##  importation BDD_feuille en format CSV
setwd("C:/IRD/Stage_inflammabilit-") #définition du répertoire de travail
BDD_leaf<-read.csv2("Data/BDD_Feuille.csv", header = TRUE) #importation de la base

# pour éventuellement visualiser la base importée
View(BDD_leaf) 

## Manipulation des données : calcul de Leaf Moisture Content (LMC) à t0 et t24,Leaf Dry Matter Content (LDMC) et Specific leaf area (SLA) 

#calcul contenu en eau "frais" (LMC_t0)
LMC_t0<-((BDD_leaf$Masse_F_T0-BDD_leaf$Masse_F_Tsec)/BDD_leaf$Masse_F_Tsec)*100
LMC_t0<-round(LMC_t0,2)
LMC_t0

#calcul contenu en eau "avant brûlage" (LMC_t24) 
LMC_t24<-((BDD_leaf$Masse_F_T24-BDD_leaf$Masse_F_Tsec)/BDD_leaf$Masse_F_Tsec)*100
LMC_t24<-round(LMC_t24,2)
LMC_t24

#calcul matière sèche (LDMC) (mg/g)
LDMC<-((BDD_leaf$Masse_F_Tsec * 1000)/BDD_leaf$Masse_F_T0)
LDMC<-round(LDMC,2)
LDMC

#calcul surface foliaire spécifique (SLA) (mm²/g)
SLA<-((BDD_leaf$Surface_F*100)/BDD_leaf$Masse_F_Tsec)
SLA<-round(SLA,2)
SLA

#Leaf Thickness (LT) = épaisseur de la feuille (mm)
LT<-((1/SLA)*LDMC)
LT<-round(LT,4)
LT

#Perte en eau des feuilles (PEF)(%)
PEF<-(1-(BDD_leaf$Masse_F_T24/BDD_leaf$Masse_F_T0))*100
PEF<-round(PEF,2)
PEF

# création d'une nouvelle base de données "feuille" calculée avec ajout des colonnes calculées
BDD_leaf_calcule<-data.frame(BDD_leaf,LMC_t0,LMC_t24,PEF,LDMC,SLA,LT)
View(BDD_leaf_calcule) #pour voir la BDD finale

###### Export de la BDD_Leaf_F CSV
write.csv2(BDD_leaf_F,"Data/BDD_leaf_calcule.csv")






###################  traitement des données de Tige  #############################

##  importation BDD_traits CSV
BDD_traits<-read.csv2("Data/BDD_Traits.csv")

# pour éventuellement visualiser la base importée
View(BDD_traits) 

## Manipulation des données : calcul du Volume et de la densité de la tige ainsi que Twig Moisture Content (TMC) et Twig Dry Matter Content (TDMC)

#calcul volume de la tige "TWig volume" (TV) (cm3)
TV<-(pi*((((BDD_traits$diametre_T_1_mm + BDD_traits$diametre_T_2_mm)/2)/2)/10)* BDD_traits$Longeur_T_cm)
TV<-round(TV,2)
TV

#calcul densitée de la tige "Twig density" (TD)(g/cm3)
TD<-(BDD_traits$Masse_T_T0/TV)
TD<-round(TD,2)
TD

#calcul TMC frais
TMC_t0<-((BDD_traits$Masse_T_T0-BDD_traits$Masse_T_Tsec)/BDD_traits$Masse_T_Tsec)*100
TMC_t0<-round(TMC_t0,2)
TMC_t0

#calcul contenu en eau tige "avant brûlage" (TMC)
TMC_t24<-((BDD_traits$Masse_T_T24-BDD_traits$Masse_T_Tsec)/BDD_traits$Masse_T_Tsec)*100
TMC_t24<-round(TMC_t24,2)
TMC_t24


#calcul matière sèche tige (TDMC) (mg/g)
TDMC<-((BDD_traits$Masse_T_Tsec * 1000)/BDD_traits$Masse_T_T0)
TDMC<-round(TDMC,2)
TDMC

#Perte en eau de la Tige (PET)(%)
PET<-(1-(BDD_traits$Masse_T_T24/BDD_traits$Masse_T_T0))*100
PET<-round(PET,2)
PET

# création d'une nouvelle base de données calculée avec ajout des colonnes
BDD_traits_calcule<-data.frame(BDD_traits,TV,TD,TMC_t0,TMC_t24,PET,TDMC)
View(BDD_traits_calcule) #pour voir la BDD finale

###### Export de la BDD traits calculée CSV

# Exporter les données avec les colonnes calculées CSV
write.csv2(BDD_traits_calcule,"Data/BDD_traits_calcule.csv")




############# Traitement des données Inflammabilité ###########################


BDD_infla<-read.csv2("Data/BDD_Inflammabilite.csv")
View(BDD_infla)

#calcul du volume des échantillons (m3)
SV<-((2/3)*pi*BDD_infla$longeur*((BDD_infla$largeur/2)*(BDD_infla$hauteur/2)))/1000000
SV<-round(SV,3)
SV

#calcul de la SD density (kg/m3)
SD<-((BDD_infla$masse/1000)/SV)
SD<-round(SD,2)
SD

#calcul de Burning Time (BT) (s) temps de combustion
BT<-(BDD_infla$temps_total - (120 + BDD_infla$DI))
BT<-round(BT)
BT

# création d'une nouvelle base de données calculée avec ajout des colonnes
BDD_infla_calcule<-data.frame(BDD_infla,BT,SV,SD)
View(BDD_infla_calcule) #pour voir la BDD finale

###### Export de la base infla calculée CSV
# Exporter les données avec les colonnes calculées CSV
write.csv2(BDD_infla_calcule,"Data/BDD_infla_calcule.csv")





######### Assembler les bases de données calculée en une seule base #########

##  importation BDD_esp CSV
BDD_esp<-read.csv2("Data/BDD_Espece.csv")
View(BDD_esp)
#importation de la base échantillonnage
BDD_ech<-read.csv2("Data/BDD_Echantillonnage.csv")
View(BDD_ech)

#assemblage : BDD esp et ech
data1<-merge(BDD_ech,BDD_esp,"ID_espece","ID_espece",all.x=T)
View(data1)

#assemblage data1 avec BDD inflammabilité calculée
data2<-merge(data1,BDD_infla_calcule,"ID_echantillon","ID_echantillon",all.x=T)
View(data2)

#assemblage data2 avec BDD traits
data3<-merge(data2,BDD_traits_calcule,"ID_echantillon","ID_echantillon",all.x=T)
View(data3)

#assemblage data2 avec BDD Feuille
BDD_finale<-merge(data3,BDD_leaf_calcule,"ID_echantillon","ID_echantillon",all.x=T)
View(BDD_finale)

#export de la BDD finale 
write.csv2(BDD_finale,"Data/BDD_finale.csv")


############# Base à l'échelle de l'espèce #######

#création de table avec moyenne, sd, min et max pour chaque variable en fonction du nom de l'espèce



moy_DI<-aggregate(DI~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_DI
moy_MT<-aggregate(MT~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_MT
moy_BB<-aggregate(BB~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_BB
moy_BT<-aggregate(BT~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_BT
moy_SV<-aggregate(SV~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_SV
moy_SD<-aggregate(SD~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_SD
moy_TMC_t0<-aggregate(TMC_t0~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_TMC_t0
moy_TMC_t24<-aggregate(TMC_t24~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_TMC_t24
moy_PET<-aggregate(PET~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_PET
moy_TDMC<-aggregate(TDMC~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_TDMC
moy_TD<-aggregate(TD~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_TD
moy_Gmin<-aggregate(Gmin~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_Gmin
moy_LMC_t0<-aggregate(LMC_t0~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_LMC_t0
moy_LMC_t24<-aggregate(LMC_t24~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_LMC_t24
moy_PEF<-aggregate(PEF~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_PEF
moy_LDMC<-aggregate(LDMC~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_LDMC
moy_Surface_F<-aggregate(Surface_F~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_Surface_F
moy_SLA<-aggregate(SLA~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_SLA
moy_LT<-aggregate(LT~Nom_scientifique, data = BDD_finale, FUN = function(x) c(moy=mean(x),sd=sd(x),min=min(x),max=max(x)))
moy_LT

##assemblage des table pour faire une seule BDD
data4<-merge(moy_DI,moy_MT,"Nom_scientifique","Nom_scientifique",all.x=T)
data5<-merge(data4,moy_BB,"Nom_scientifique","Nom_scientifique",all.x=T)
data6<-merge(data5,moy_BT,"Nom_scientifique","Nom_scientifique",all.x=T)
data7<-merge(data6,moy_SV,"Nom_scientifique","Nom_scientifique",all.x=T)
data8<-merge(data7,moy_SD,"Nom_scientifique","Nom_scientifique",all.x=T)
data9<-merge(data8,moy_TMC_t0,"Nom_scientifique","Nom_scientifique",all.x=T)
data10<-merge(data9,moy_TMC_t24,"Nom_scientifique","Nom_scientifique",all.x=T)
data11<-merge(data10,moy_PET,"Nom_scientifique","Nom_scientifique",all.x=T)
data12<-merge(data11,moy_TDMC,"Nom_scientifique","Nom_scientifique",all.x=T)
data13<-merge(data12,moy_Gmin,"Nom_scientifique","Nom_scientifique",all.x=T)
data14<-merge(data13,moy_LMC_t0,"Nom_scientifique","Nom_scientifique",all.x=T)
data15<-merge(data14,moy_LMC_t24,"Nom_scientifique","Nom_scientifique",all.x=T)
data16<-merge(data15,moy_PEF,"Nom_scientifique","Nom_scientifique",all.x=T)
data17<-merge(data16,moy_LDMC,"Nom_scientifique","Nom_scientifique",all.x=T)
data18<-merge(data17,moy_Surface_F,"Nom_scientifique","Nom_scientifique",all.x=T)
data19<-merge(data18,moy_SLA,"Nom_scientifique","Nom_scientifique",all.x=T)
BDD_finale_esp<-merge(data19,moy_LT,"Nom_scientifique","Nom_scientifique",all.x=T)
View(BDD_finale_esp)

data19

# Arrondir toutes les colonnes numériques à 2 décimales
ncol(BDD_finale_esp)
BDD_finale_esp[] <- round((BDD_finale_esp[,2:73]), 2)
BDD_finale_esp
########################################

# Exemple de structure de ton tableau

# Initialiser les espèces uniques
variables <- BDD_finale[, sapply(BDD_finale, is.numeric)]
variables

temp<-BDD_finale[,c("MT","BT","DI")]

moy<-aggregate(temp[,1]~Nom_scientifique, data = BDD_finale, FUN = mean)
moy[,2]<-round(moy[,2],2)
colnames(moy)[2]<-colnames(temp)[1]

ET<-aggregate(temp[,1]~Nom_scientifique, data = BDD_finale, FUN = sd)
ET[,2]<-round(ET[,2],2)
colnames(ET)[2]<-colnames(temp)[1]

for (i in 2:ncol(temp)){
  temp2<-aggregate(temp[,i]~Nom_scientifique, data = BDD_finale, FUN = mean)
  temp2[,2]<-round(temp2[,2],2)
  colnames(temp2)[2]<-colnames(temp)[i]
  moy<-merge(x=moy,y=temp2,by.x="Nom_scientifique",by.y="Nom_scientifique",all.x=T,all.y=T)
  temp2<-aggregate(temp[,i]~Nom_scientifique, data = BDD_finale, FUN = sd)
  temp2[,2]<-round(temp2[,2],2)
  colnames(temp2)[2]<-colnames(temp)[i]
  ET<-merge(x=ET,y=temp2,by.x="Nom_scientifique",by.y="Nom_scientifique",all.x=T,all.y=T)
}

DATA<-list(moy,ET)

resultats








# Créer une liste vide pour stocker les résultats
resultats <- data.frame()

# Boucle for pour chaque espèce
for (e in especes) {
  sous_ensemble <- data[data$espece == e, ]  
  
  stats <- data.frame(
    espece = e,
    moyenne_temp = mean(sous_ensemble$temperature, na.rm = TRUE),
    sd_temp = sd(sous_ensemble$temperature, na.rm = TRUE),
    min_temp = min(sous_ensemble$temperature, na.rm = TRUE),
    max_temp = max(sous_ensemble$temperature, na.rm = TRUE),
    
    moyenne_temps = mean(sous_ensemble$temps, na.rm = TRUE),
    sd_temps = sd(sous_ensemble$temps, na.rm = TRUE),
    min_temps = min(sous_ensemble$temps, na.rm = TRUE),
    max_temps = max(sous_ensemble$temps, na.rm = TRUE),
    
    moyenne_masse = mean(sous_ensemble$masse, na.rm = TRUE),
    sd_masse = sd(sous_ensemble$masse, na.rm = TRUE),
    min_masse = min(sous_ensemble$masse, na.rm = TRUE),
    max_masse = max(sous_ensemble$masse, na.rm = TRUE)
  )
  
  resultats <- rbind(resultats, stats)
}

# Afficher le tableau final
print(resultats)

ncol(resultats)

#export de la BDD finale 
write.csv2(BDD_finale_esp,"Data/BDD_finale_esp.csv")




#arrondir les valeurs
BDD_finale_esp<-round(c(2:50),2)
View(BDD_finale_esp)
names(BDD_finale_esp)[c(2:8)] = c("moy_DI", "sd_DI","min_DI","max_DI","moy_MT","sd_MT","min_MT")





# Écart-type
ecarts<-aggregate(cbind(DI,MT,BB,BT,SV,SD,TMC_t0,TMC_t24,PET,TDMC,TD,Gmin,LMC_t0,LMC_t24,PEF,LDMC,Surface_F,SLA,LT)~Nom_scientifique, data = BDD_finale, FUN = sd)
ecarts



