
#### TRAITEMENT DONNEE STAGE INFLAMMABILITE ####

############# Traitement des données feuille ###########################


##  importation BDD_feuille en format CSV
setwd("C:/IRD/Stage_inflammabilit-") #définition du répertoire de travail
BDD_leaf<-read.csv2("Data/BDD_Feuille.csv", header = TRUE) #importation de la base

# pour éventuellement visualiser la base importée
BDD_leaf

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
BDD_leaf_calcule #pour voir la BDD finale

###### Export de la BDD_Leaf_F CSV
write.csv2(BDD_leaf_F,"Data/BDD_leaf_calcule.csv")






###################  traitement des données de Tige  #############################

##  importation BDD_traits CSV
BDD_traits<-read.csv2("Data/BDD_Traits.csv")

# pour éventuellement visualiser la base importée
BDD_traits 

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

#calcul du diamètre de la tige
TDIA<-((BDD_traits$diametre_T_1_mm + BDD_traits$diametre_T_2_mm)/2)
TDIA

# création d'une nouvelle base de données calculée avec ajout des colonnes
BDD_traits_calcule<-data.frame(BDD_traits,TV,TD,TMC_t0,TMC_t24,TDMC,TDIA)
BDD_traits_calcule #pour voir la BDD finale

###### Export de la BDD traits calculée CSV

# Exporter les données avec les colonnes calculées CSV
write.csv2(BDD_traits_calcule,"Data/BDD_traits_calcule.csv")




############# Traitement des données Inflammabilité ###########################


BDD_infla<-read.csv2("Data/BDD_Inflammabilite.csv")
BDD_infla

# Exemple avec trois ID
ids_a_remplacer <- c("02_5", "13_6", "30_5","39_6","42_5")

# Boucle simple pour remplacer ligne par ligne
for (id in ids_a_remplacer) {
  lignes <- which(BDD_infla$ID_echantillon == id)
  BDD_infla[lignes, c("DI", "BB", "temps_total", "MT","DI_test","BB_test")] <- NA
}
View(BDD_infla)


#calcul de la SD density (kg/m3)
SD<-((BDD_infla$masse/1000)/(((2/3)*pi*BDD_infla$longeur*((BDD_infla$largeur/2)*(BDD_infla$hauteur/2)))/1000000))
SD<-round(SD,2)
SD

#calcul de Burning Time (BT) (s) temps de combustion
BT<-(BDD_infla$temps_total - (120 + BDD_infla$DI))
BT<-round(BT)
BT


#calcul de l'inverse de DI
BDD_infla$DI<-(10-BDD_infla$DI)
BDD_infla$DI


# création d'une nouvelle base de données calculée avec ajout des colonnes
BDD_infla_calcule<-data.frame(BDD_infla,BT,SD)
BDD_infla_calcule #pour voir la BDD finale

###### Export de la base infla calculée CSV
# Exporter les données avec les colonnes calculées CSV
write.csv2(BDD_infla_calcule,"Data/BDD_infla_calcule.csv")





######### Assembler les bases de données en une seule base avec TOUTES les infos  #########
print(BDD_esp$Nom.scientifique)

##  importation BDD_esp CSV
BDD_esp<-read.csv2("Data/BDD_Espece.csv")
names(BDD_esp)[which(names(BDD_esp) == "Nom.scientifique")] <- "Nom_scientifique"

BDD_esp
#importation de la base échantillonnage
BDD_ech<-read.csv2("Data/BDD_Echantillonnage.csv")
BDD_ech

#assemblage : BDD esp et ech
data1 <- merge(BDD_ech, BDD_esp, by.x = "ID_espece", by.y = "ID.espece", all.x = TRUE)
data1

#assemblage data1 avec BDD inflammabilité calculée
data2<-merge(data1,BDD_infla_calcule,"ID_echantillon","ID_echantillon",all.x=T)
data2

#assemblage data2 avec BDD traits
data3<-merge(data2,BDD_traits_calcule,"ID_echantillon","ID_echantillon",all.x=T)
data3

#assemblage data2 avec BDD Feuille
BDD_finale<-merge(data3,BDD_leaf_calcule,"ID_echantillon","ID_echantillon",all.x=T)
BDD_finale

#export de la BDD finale 
write.csv2(BDD_finale,"Data/BDD_finale.csv")


############# création d'une BDD avec seulement les infos pour les analyses ############

BDD_ana_ech<-subset(BDD_finale, select=c(Nom_scientifique, Milieu_recolte,ID_espece,ID_echantillon,ID_Feuille,DI,BT,MT,BB,Nb_ramifications,SV,SD,TMC_t0,TMC_t24,TDMC,TD,TDIA,Gmin,LMC_t0,LMC_t24,PEF,LDMC,Surface_F,SLA,LT))
head(BDD_ana_ech)
dim(BDD_ana_ech)
#export de la BDD 
write.csv2(BDD_ana_ech,"Data/BDD_ana_ech.csv")





############# Base à l'échelle de l'échantillon #######################

#création de table avec moyenne et sd pour chaque variable en fonction du nom de l'espèce
temp<-BDD_ana_ech[,6:25] ###sélection des colonnes comprenant les variables pour les intégrer dans la boucle
temp

#création d'une bdd d'origine pour moyenne (sert pour merge)
BDD_moy_ech <- aggregate(temp[,1] ~ ID_echantillon + Nom_scientifique + ID_espece + Milieu_recolte, data = BDD_ana_ech, FUN = mean, na.rm = TRUE)
BDD_moy_ech[,5] <- round(BDD_moy_ech[,5], 2)
colnames(BDD_moy_ech)[5] <- colnames(temp)[1]

#création d'un bdd d'origine pour sd (sert pour merge)
BDD_sd_ech <- aggregate(temp[,1] ~ ID_echantillon + Nom_scientifique + ID_espece+ Milieu_recolte, data = BDD_ana_ech, FUN = sd, na.rm = TRUE)
BDD_sd_ech[,5] <- round(BDD_sd_ech[,5], 2)
colnames(BDD_sd_ech)[5] <- colnames(temp)[1]

# Boucle pour les calcul des moyennes et écart-types
for (i in 2:ncol(temp)) {
  # Moyenne
  temp_moy <- aggregate(temp[, i] ~ ID_echantillon + Nom_scientifique + ID_espece+ Milieu_recolte,data = BDD_ana_ech, FUN = mean, na.rm = TRUE)
  temp_moy[,5] <- round(temp_moy[,5], 2)
  colnames(temp_moy)[5] <- colnames(temp)[i]
  BDD_moy_ech <- merge(BDD_moy_ech, temp_moy, by = c("ID_echantillon", "Nom_scientifique", "ID_espece","Milieu_recolte"), all = TRUE)
  
  # Écart-type
  temp_sd <- aggregate(temp[, i] ~ ID_echantillon + Nom_scientifique + ID_espece+ Milieu_recolte,data = BDD_ana_ech, FUN = sd, na.rm = TRUE)
  temp_sd[,5] <- round(temp_sd[,5], 2)
  colnames(temp_sd)[5] <- colnames(temp)[i]
  BDD_sd_ech <- merge(BDD_sd_ech, temp_sd,by = c("ID_echantillon", "Nom_scientifique", "ID_espece","Milieu_recolte"),all = TRUE)
  }

# Résultats
head(BDD_moy_ech)
head(BDD_sd_ech)

# Export des données
write.csv2(BDD_moy_ech, "Data/BDD_moy_ech.csv", row.names = FALSE)
write.csv2(BDD_sd_ech, "Data/BDD_sd_ech.csv", row.names = FALSE)




############# Base à l'échelle de l'espèce ######################################

#création d'un bdd d'origine pour moyenne (sert pour merge)
BDD_moy_esp <- aggregate(temp[,1] ~ Nom_scientifique + ID_espece + Milieu_recolte, data = BDD_finale, FUN = mean, na.rm = TRUE)
BDD_moy_esp[,4] <- round(BDD_moy_esp[,4], 2)
colnames(BDD_moy_esp)[4] <- colnames(temp)[1]

#création d'un bdd d'origine pour sd (sert pour merge)
BDD_sd_esp <- aggregate(temp[,1] ~ Nom_scientifique + ID_espece + Milieu_recolte, data = BDD_finale, FUN = sd, na.rm = TRUE)
BDD_sd_esp[,4] <- round(BDD_sd_esp[,4], 2)
colnames(BDD_sd_esp)[4] <- colnames(temp)[1]

#Boucle pour les calcul des moyennes et écart-types
for (i in 2:ncol(temp)) {
  
  # Moyenne
  temp_moy_esp <- aggregate(temp[, i] ~ Nom_scientifique + ID_espece+ Milieu_recolte, data = BDD_finale, FUN = mean, na.rm = TRUE)
  temp_moy_esp[,4] <- round(temp_moy_esp[,4], 2)
  colnames(temp_moy_esp)[4] <- colnames(temp)[i]
  BDD_moy_esp <- merge(BDD_moy_esp, temp_moy_esp, by = c("Nom_scientifique", "ID_espece", "Milieu_recolte"), all = TRUE)
  
  # Ecart-type
  temp_sd_esp <- aggregate(temp[, i] ~ Nom_scientifique + ID_espece+ Milieu_recolte, data = BDD_finale, FUN = sd, na.rm = TRUE)
  temp_sd_esp[,4] <- round(temp_sd_esp[,4], 2)
  colnames(temp_sd_esp)[4] <- colnames(temp)[i]
  BDD_sd_esp <- merge(BDD_sd_esp, temp_sd_esp, by = c("Nom_scientifique", "ID_espece" , "Milieu_recolte"), all = TRUE)
}

#Résultat
head(BDD_moy_esp)
head(BDD_sd_esp)

#export de la BDD 
write.csv2(BDD_moy_esp,"Data/BDD_moy_esp.csv")
write.csv2(BDD_sd_esp, "Data/BDD_sd_ech.csv")





################ NOMBRE D'ESPECES PAR MILIEU ###########
library(dplyr)
library(tidyr)
library(ggplot2)

# Comptage des espèces par milieu
df_plot <- aggregate(Nom_scientifique ~ Milieu_recolte, data = BDD_esp, FUN = length)
colnames(df_plot) <- c("Milieu", "Nb_especes") 

# Tri du plus grand au plus petit
df_plot <- df_plot[order(df_plot$Nb_especes, decreasing = TRUE), ]

# Diagramme
ggplot(df_plot, aes(x = reorder(Milieu, -Nb_especes), y = Nb_especes)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  labs(x = "Milieu", y = "Nombre d'espèces", title = "Nombre d'espèces par milieu") +
  theme_minimal()


################ NOMBRE D'ESPECES PAR statut ###########

# comptage des espèces par milieu
df_statut <- aggregate(BDD_esp$Nom_scientifique ~ BDD_esp$Statut, FUN = length)
colnames(df_statut)[1] <- "Statut"
colnames(df_statut)[2] <- "Nb_especes"

# tri du plus grand au plus petit
df_statut <- df_statut[order(df_statut$Nb_especes, decreasing = TRUE), ]
names(df_statut)

# Étape 4 : diagramme
ggplot(df_statut, aes(x = reorder(Statut, -Nb_especes), y = Nb_especes)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  labs(x = "Milieu", y = "Nombre d'espèces", title = "Nombre d'espèces par milieu") 

