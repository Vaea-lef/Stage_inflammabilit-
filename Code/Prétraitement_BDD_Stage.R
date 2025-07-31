
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
SLA
SLA<-round(SLA,2)
SLA

#Leaf Thickness (LT) = épaisseur de la feuille (mm)
LT<-(((1/SLA)*LDMC)*1000)
LT
LT<-round(LT,2)
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
BDD_infla[BDD_infla$ID_echantillon == "21_2", "MT"] <- 150
BDD_infla[BDD_infla$ID_echantillon == "21_2", "temps_total"] <- 130

#calcul de la SD density (kg/m3)
SD<-((BDD_infla$masse/1000)/(((2/3)*pi*BDD_infla$longeur*((BDD_infla$largeur/2)*(BDD_infla$hauteur/2)))/1000000))
SD<-round(SD,2)
SD

#calcul de Burning Time (BT) (s) temps de combustion
BT<-(BDD_infla$temps_total - (120 + BDD_infla$DI))
BT<-round(BT)
BT

#calcul du BT avec extra case
BT_test<-(BDD_infla$temps_total - (120 + BDD_infla$DI_test))
BT_test<-round(BT_test)
BT_test

#calcul de l'inverse de DI_test
BDD_infla$DI<-(10-BDD_infla$DI)
BDD_infla$DI



#calcul de l'inverse de DI_test
score_DI<-(10-BDD_infla$DI_test)
score_DI

#calcul fréquence ignition 
FI <- ifelse(BDD_infla$MT > 150, 1, 0)
FI

#calcul VPD
#VPD = VPsaturated − VPair, with VPsaturated = 0.61 08 × exp((17.27 × T°)/(T° + 237.2)) and VPair = HR/(100 × VPsaturated)
VPsat<-0.6108*exp((17.27*BDD_infla$T_ambiante)/(BDD_infla$T_ambiante+237.2))
VPair<- (BDD_infla$Humidite/100)*VPsat
VPD<-VPsat-VPair
summary(VPD)

# création d'une nouvelle base de données calculée avec ajout des colonnes
BDD_infla_calcule<-data.frame(BDD_infla,score_DI,BT, BT_test,FI,SD,VPD)
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

colnames(BDD_finale)


############# création d'une BDD avec seulement les infos pour les analyses ############
BDD_ana_ech<-subset(BDD_finale, select=c(Nom_scientifique, Milieu_recolte,ID_espece,ID_echantillon,ID_Feuille,T_ambiante,Vent,Humidite,DI,DI_test,score_DI,BT,BT_test,MT,BB,BB_test,FI,Nb_ramifications,SD,TMC_t0,TMC_t24,TDMC,TD,TDIA,Gmin,LMC_t0,LMC_t24,PEF,LDMC,Surface_F,SLA,LT,VPD))
head(BDD_ana_ech)
dim(BDD_ana_ech)
#export de la BDD 
write.csv2(BDD_ana_ech,"Data/BDD_ana_ech.csv")










############# Base à l'échelle de l'échantillon #######################

#création de table avec moyenne et sd pour chaque variable en fonction du nom de l'espèce
temp<-BDD_ana_ech[,6:33] ###sélection des colonnes comprenant les variables pour les intégrer dans la boucle
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



# Vecteur des ID à supprimer
ids_a_supprimer <- c("02_5", "13_6", "30_5", "39_6", "42_5")

# Trouver les indices à supprimer
lignes_a_supprimer <- c()
for (id in ids_a_supprimer) {
  lignes_a_supprimer <- c(lignes_a_supprimer, which(BDD_moy_ech$ID_echantillon == id))
}

# Supprimer ces lignes dans l'objet d'origine
BDD_moy_ech <- BDD_moy_ech[-lignes_a_supprimer, ]

# Résultats
head(BDD_moy_ech)
head(BDD_sd_ech)

# Export des données
write.csv2(BDD_moy_ech, "Data/BDD_moy_ech.csv", row.names = FALSE)
write.csv2(BDD_sd_ech, "Data/BDD_sd_ech.csv", row.names = FALSE)

View(BDD_moy_ech)


####### suppression des échantillons pour MT
BDD_moy_ech1 <- BDD_moy_ech[BDD_moy_ech$MT != 150, ]
View(BDD_moy_ech1)

####### suppression des échantillons pour reste composates
BDD_moy_ech2 <- BDD_moy_ech[BDD_moy_ech$DI_test != 10, ]
View(BDD_moy_ech2)

# Export des données
write.csv2(BDD_moy_ech1, "Data/BDD_moy_ech1.csv", row.names = FALSE)
write.csv2(BDD_moy_ech2, "Data/BDD_moy_ech2.csv", row.names = FALSE)
dim(BDD_moy_ech)















############# Base à l'échelle de l'espèce (complète) ######################################
#création de table avec moyenne et sd pour chaque variable en fonction du nom de l'espèce
tem1<-BDD_moy_ech[,5:32] ###sélection des colonnes comprenant les variables pour les intégrer dans la boucle
tem1
#création d'un bdd d'origine pour moyenne (sert pour merge)
BDD_moy_esp <- aggregate(tem1[,1] ~ Nom_scientifique + ID_espece + Milieu_recolte, data = BDD_moy_ech, FUN = mean, na.rm = TRUE)
BDD_moy_esp[,4] <- round(BDD_moy_esp[,4], 2)
colnames(BDD_moy_esp)[4] <- colnames(tem1)[1]

#création d'un bdd d'origine pour sd (sert pour merge)
BDD_sd_esp <- aggregate(tem1[,1] ~ Nom_scientifique + ID_espece + Milieu_recolte, data = BDD_moy_ech, FUN = sd, na.rm = TRUE)
BDD_sd_esp[,4] <- round(BDD_sd_esp[,4], 2)
colnames(BDD_sd_esp)[4] <- colnames(tem1)[1]

#Boucle pour les calcul des moyennes et écart-types
for (i in 2:ncol(tem1)) {
  
  # Moyenne
  tem1_moy_esp <- aggregate(tem1[, i] ~ Nom_scientifique + ID_espece+ Milieu_recolte, data = BDD_moy_ech, FUN = mean, na.rm = TRUE)
  tem1_moy_esp[,4] <- round(tem1_moy_esp[,4], 2)
  colnames(tem1_moy_esp)[4] <- colnames(tem1)[i]
  BDD_moy_esp <- merge(BDD_moy_esp, tem1_moy_esp, by = c("Nom_scientifique", "ID_espece", "Milieu_recolte"), all = TRUE)
  
  # Ecart-type
  tem1_sd_esp <- aggregate(tem1[, i] ~ Nom_scientifique + ID_espece+ Milieu_recolte, data = BDD_moy_ech, FUN = sd, na.rm = TRUE)
  tem1_sd_esp[,4] <- round(tem1_sd_esp[,4], 2)
  colnames(tem1_sd_esp)[4] <- colnames(tem1)[i]
  BDD_sd_esp <- merge(BDD_sd_esp, tem1_sd_esp, by = c("Nom_scientifique", "ID_espece" , "Milieu_recolte"), all = TRUE)
}


BDD_FI_esp <- aggregate(FI ~ Nom_scientifique + ID_espece + Milieu_recolte, 
                        data = BDD_moy_ech, 
                        FUN = sum, 
                        na.rm = TRUE)
colnames(BDD_FI_esp)[4] <- "Nb_FI"

BDD_FI_esp


BDD_moy_esp <- merge(BDD_moy_esp, BDD_FI_esp, 
                     by = c("Nom_scientifique", "ID_espece", "Milieu_recolte"), 
                     all.x = TRUE)

#Résultat
head(BDD_moy_esp)
head(BDD_sd_esp)

View(BDD_sd_esp)

#export de la BDD 
write.csv2(BDD_moy_esp,"Data/BDD_moy_esp.csv")
write.csv2(BDD_sd_esp, "Data/BDD_sd_esp.csv")




############# Base à l'échelle de l'espèce (nettoyée MT) ######################################
#création de table avec moyenne et sd pour chaque variable en fonction du nom de l'espèce
tem2<-BDD_moy_ech1[,5:32] ###sélection des colonnes comprenant les variables pour les intégrer dans la boucle
tem2
#création d'un bdd d'origine pour moyenne (sert pour merge)
BDD_moy_esp1 <- aggregate(tem2[,1] ~ Nom_scientifique + ID_espece + Milieu_recolte, data = BDD_moy_ech1, FUN = mean, na.rm = TRUE)
BDD_moy_esp1[,4] <- round(BDD_moy_esp1[,4], 2)
colnames(BDD_moy_esp1)[4] <- colnames(tem2)[1]

#création d'un bdd d'origine pour sd (sert pour merge)
BDD_sd_esp <- aggregate(tem2[,1] ~ Nom_scientifique + ID_espece + Milieu_recolte, data = BDD_moy_ech1, FUN = sd, na.rm = TRUE)
BDD_sd_esp[,4] <- round(BDD_sd_esp[,4], 2)
colnames(BDD_sd_esp)[4] <- colnames(tem2)[1]

#Boucle pour les calcul des moyennes et écart-types
for (i in 2:ncol(tem2)) {
  
  # Moyenne
  tem2_moy_esp <- aggregate(tem2[, i] ~ Nom_scientifique + ID_espece+ Milieu_recolte, data = BDD_moy_ech1, FUN = mean, na.rm = TRUE)
  tem2_moy_esp[,4] <- round(tem2_moy_esp[,4], 2)
  colnames(tem2_moy_esp)[4] <- colnames(tem2)[i]
  BDD_moy_esp1 <- merge(BDD_moy_esp1, tem2_moy_esp, by = c("Nom_scientifique", "ID_espece", "Milieu_recolte"), all = TRUE)
  
  # Ecart-type
  tem2_sd_esp <- aggregate(tem2[, i] ~ Nom_scientifique + ID_espece+ Milieu_recolte, data = BDD_moy_ech1, FUN = sd, na.rm = TRUE)
  tem2_sd_esp[,4] <- round(tem2_sd_esp[,4], 2)
  colnames(tem2_sd_esp)[4] <- colnames(tem2)[i]
  BDD_sd_esp <- merge(BDD_sd_esp, tem2_sd_esp, by = c("Nom_scientifique", "ID_espece" , "Milieu_recolte"), all = TRUE)
}


BDD_FI_esp <- aggregate(FI ~ Nom_scientifique + ID_espece + Milieu_recolte, 
                        data = BDD_moy_ech11, 
                        FUN = sum, 
                        na.rm = TRUE)
colnames(BDD_FI_esp)[4] <- "Nb_FI"

BDD_FI_esp


BDD_moy_esp1 <- merge(BDD_moy_esp1, BDD_FI_esp, 
                     by = c("Nom_scientifique", "ID_espece", "Milieu_recolte"), 
                     all.x = TRUE)

#Résultat
head(BDD_moy_esp1)
View(BDD_moy_esp1)

#export de la BDD 
write.csv2(BDD_moy_esp1,"Data/BDD_moy_esp1.csv")



############# Base à l'échelle de l'espèce (nettoyée reste) ######################################
#création de table avec moyenne et sd pour chaque variable en fonction du nom de l'espèce
tem2<-BDD_moy_ech2[,5:32] ###sélection des colonnes comprenant les variables pour les intégrer dans la boucle
tem2
#création d'un bdd d'origine pour moyenne (sert pour merge)
BDD_moy_esp2 <- aggregate(tem2[,1] ~ Nom_scientifique + ID_espece + Milieu_recolte, data = BDD_moy_ech2, FUN = mean, na.rm = TRUE)
BDD_moy_esp2[,4] <- round(BDD_moy_esp2[,4], 2)
colnames(BDD_moy_esp2)[4] <- colnames(tem2)[1]

#création d'un bdd d'origine pour sd (sert pour merge)
BDD_sd_esp <- aggregate(tem2[,1] ~ Nom_scientifique + ID_espece + Milieu_recolte, data = BDD_moy_ech2, FUN = sd, na.rm = TRUE)
BDD_sd_esp[,4] <- round(BDD_sd_esp[,4], 2)
colnames(BDD_sd_esp)[4] <- colnames(tem2)[1]

#Boucle pour les calcul des moyennes et écart-types
for (i in 2:ncol(tem2)) {
  
  # Moyenne
  tem2_moy_esp <- aggregate(tem2[, i] ~ Nom_scientifique + ID_espece+ Milieu_recolte, data = BDD_moy_ech2, FUN = mean, na.rm = TRUE)
  tem2_moy_esp[,4] <- round(tem2_moy_esp[,4], 2)
  colnames(tem2_moy_esp)[4] <- colnames(tem2)[i]
  BDD_moy_esp2 <- merge(BDD_moy_esp2, tem2_moy_esp, by = c("Nom_scientifique", "ID_espece", "Milieu_recolte"), all = TRUE)
  
  # Ecart-type
  tem2_sd_esp <- aggregate(tem2[, i] ~ Nom_scientifique + ID_espece+ Milieu_recolte, data = BDD_moy_ech2, FUN = sd, na.rm = TRUE)
  tem2_sd_esp[,4] <- round(tem2_sd_esp[,4], 2)
  colnames(tem2_sd_esp)[4] <- colnames(tem2)[i]
  BDD_sd_esp <- merge(BDD_sd_esp, tem2_sd_esp, by = c("Nom_scientifique", "ID_espece" , "Milieu_recolte"), all = TRUE)
}


BDD_FI_esp <- aggregate(FI ~ Nom_scientifique + ID_espece + Milieu_recolte, 
                        data = BDD_moy_ech2, 
                        FUN = sum, 
                        na.rm = TRUE)
colnames(BDD_FI_esp)[4] <- "Nb_FI"

BDD_FI_esp


BDD_moy_esp2 <- merge(BDD_moy_esp2, BDD_FI_esp, 
                      by = c("Nom_scientifique", "ID_espece", "Milieu_recolte"), 
                      all.x = TRUE)

#Résultat
head(BDD_moy_esp2)
View(BDD_moy_esp2)

#export de la BDD 
write.csv2(BDD_moy_esp2,"Data/BDD_moy_esp2.csv")




















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






# Vecteur des ID à supprimer
ids_a_supprimer <- c("02_5", "13_6", "30_5", "39_6", "42_5")
# Initialisation d'un vecteur d'indices à supprimer
lignes_a_supprimer <- c()
lignes_a_supprimer <- unique(lignes_a_supprimer)  # au cas où il y aurait des doublons


# Boucle pour trouver les indices à supprimer
for (id in ids_a_supprimer) {
  lignes_a_supprimer <- c(lignes_a_supprimer, which(BDD_infla$ID_echantillon == id))
}

# Suppression des lignes
BDD_infla <- BDD_infla[-lignes_a_supprimer, ]
BDD_infla <- BDD_infla[BDD_infla$MT != 150, ]



ID_a_supprimer <- c("02_5", "13_6", "30_5", "39_6", "42_5")

BDD_infla <- BDD_infla[
  !(BDD_infla$ID_echantillon == "02_5" |
      BDD_infla$ID_echantillon == "13_6" |
      BDD_infla$ID_echantillon == "30_5" |
      BDD_infla$ID_echantillon == "39_6" |
      BDD_infla$ID_echantillon == "42_5"), ]










########### coef de variation ################
BDD_coef<-merge(BDD_moy_esp,BDD_sd_esp,"ID_espece","ID_espece",all.x=T)
View(BDD_coef)

BDD_coef$coefMT<-(BDD_coef$MT.y/BDD_coef$MT.x)*100
BDD_coef$coefBB<-(BDD_coef$BB_test.y/BDD_coef$BB_test.x)*100
BDD_coef$coefBT<-(BDD_coef$BT_test.y/BDD_coef$BT_test.x)*100
BDD_coef$coefDI<-(BDD_coef$DI_test.y/BDD_coef$DI_test.x)*100


BDD_coef$coefSLA<-(BDD_coef$SLA.y/BDD_coef$SLA.x)*100
BDD_coef$coefLA<-(BDD_coef$Surface_F.y/BDD_coef$Surface_F.x)*100
BDD_coef$coefLDMC<-(BDD_coef$LDMC.y/BDD_coef$LDMC.x)*100
BDD_coef$coefLMC_t24<-(BDD_coef$LMC_t24.y/BDD_coef$LMC_t24.x)*100
BDD_coef$coefTD<-(BDD_coef$TD.y/BDD_coef$TD.x)*100
BDD_coef$coefLT<-(BDD_coef$LT.y/BDD_coef$LT.x)*100

write.csv2(BDD_coef,"Data/BDD_coef.csv")
