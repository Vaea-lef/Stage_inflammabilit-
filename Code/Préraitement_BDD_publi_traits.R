
####################################################################
#################  PRETRAITEMENT BDD PUBLI TRAITS  #################
####################################################################




############# Traitement des données feuille ###########################

###  importation BDD_feuille en format CSV
setwd("C:/IRD/Stage_inflammabilit-") #définition du répertoire de travail
BDD_leaf<-read.csv2("Data/BDD_Feuille.csv", header = TRUE) #importation de la base


### Manipulation des données : calcul de Leaf Moisture Content (LMC) à t0 et t24,Leaf Dry Matter Content (LDMC), Specific leaf area (SLA), Leaf Thikness (LT) et Perte d'eau (PEF)

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


### création d'une nouvelle base de données "feuille" calculée avec ajout des colonnes calculées
BDD_leaf_calcule_publi<-data.frame(BDD_leaf,LMC_t0,LMC_t24,PEF,LDMC,SLA,LT)
head(BDD_leaf_calcule_publi) #pour voir si les colonnes sont bien présentes

###### Export de la BDD_Leaf_F CSV
write.csv2(BDD_leaf_calcule_publi,"Data/Publi/BDD_leaf_calcule_publi.csv")






###################  traitement des données de Tige  #############################

###  importation BDD_traits CSV
BDD_traits<-read.csv2("Data/BDD_Traits.csv")

### Manipulation des données : calcul du Volume et de la densité de la tige (TD) ainsi que Twig Moisture Content (TMC) et Twig Dry Matter Content (TDMC)

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

### création d'une nouvelle base de données calculée avec ajout des colonnes
BDD_traits_calcule_publi<-data.frame(BDD_traits,TV,TD,TMC_t0,TMC_t24,TDMC,TDIA)
head(BDD_traits_calcule_publi) #pour voir si les colonnes sont bien présentes

###### Export de la BDD traits calculée CSV

# Exporter les données avec les colonnes calculées CSV
write.csv2(BDD_traits_calcule_publi,"Data/Publi/BDD_traits_calcule_publi.csv")





############# Traitement des données Inflammabilité ###########################

###  importation BDD_traits CSV
BDD_infla<-read.csv2("Data/BDD_Inflammabilite.csv")

### Manipulation des données :calcul de la bulk density (SV), burning time, 

#calcul de la BD density (kg/m3)
BD<-((BDD_infla$masse/1000)/(((2/3)*pi*BDD_infla$longeur*((BDD_infla$largeur/2)*(BDD_infla$hauteur/2)))/1000000))
BD<-round(BD,2)
BD

#calcul de Burning Time (BT) (s) temps de combustion
BDD_infla$temps_total <- as.numeric(BDD_infla$temps_total)
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
BDD_infla_calcule_publi<-data.frame(BDD_infla,score_DI,BT, BT_test,FI,BD,VPD)
head(BDD_infla_calcule_publi) #pour voir la BDD finale

###### Export de la base infla calculée CSV
# Exporter les données avec les colonnes calculées CSV
write.csv2(BDD_infla_calcule_publi,"Data/Publi/BDD_infla_calcule_publi.csv")





######### Assembler les bases de données en une seule base avec TOUTES les infos  #########

##  importation base espèces
BDD_esp<-read.csv2("Data/BDD_Espece.csv")
#rennomer colonne
names(BDD_esp)[which(names(BDD_esp) == "Nom.scientifique")] <- "Nom_scientifique"
#afficher BDD esp
BDD_esp

#importation de la base échantillonnage
BDD_echantillonnage<-read.csv2("Data/BDD_Echantillonnage.csv")
BDD_echantillonnage

#assemblage : BDD esp et echantillonnage
data1 <- merge(BDD_echantillonnage, BDD_esp, by.x = "ID_espece", by.y = "ID.espece", all.x = TRUE)
data1

#assemblage data1 avec BDD inflammabilité calculée
data2<-merge(data1,BDD_infla_calcule_publi,"ID_echantillon","ID_echantillon",all.x=T)
data2

#assemblage data2 avec BDD traits
data3<-merge(data2,BDD_traits_calcule_publi,"ID_echantillon","ID_echantillon",all.x=T)
data3

#assemblage data2 avec BDD Feuille
BDD_finale_publi<-merge(data3,BDD_leaf_calcule_publi,"ID_echantillon","ID_echantillon",all.x=T)
BDD_finale_publi

# Couper le nom scientifique en morceaux au niveau de l'espace
morceaux <- strsplit(BDD_finale_publi$Nom_scientifique, " ")

# Le 1er morceau = genre
BDD_finale_publi$Genre <- sapply(morceaux, function(x) x[1])

# Le 2ème morceau = espèce
BDD_finale_publi$Espece <- sapply(morceaux, function(x) x[2])

# Vérifier
head(BDD_finale_publi[, c("Nom_scientifique", "Genre", "Espece")])

#export de la BDD finale 
write.csv2(BDD_finale_publi,"Data/Publi/BDD_finale_publi.csv")










############# création d'une BDD avec seulement les infos pour les analyses ############

##### D'abord on nettoie (on enlève les 3 variétés de Bourao)
# Vecteur des ID à supprimer (Bourao R, G, P)
var_a_supprimer <- c("Hibiscus tiliaceus R", "Hibiscus tiliaceus G", "Hibiscus tiliaceus P")

# Trouver leslignes des ids à supprimer
lignes_a_supprimer <- c()
for (id in var_a_supprimer) {lignes_a_supprimer <- c(lignes_a_supprimer, which(BDD_finale_publi$Nom_scientifique == id))}

# Supprimer ces lignes dans l'objet d'origine
BDD_finale_publi <- BDD_finale_publi[-lignes_a_supprimer, ]
dim(BDD_finale_publi)

# Nom espèce Bourao V changement 
BDD_finale_publi$Nom_scientifique[
  BDD_finale_publi$Nom_scientifique == "Hibiscus tiliaceus V"
] <- "Hibiscus tiliaceus"



### on garde que les variables intéréssantes pour l'analyse
BDD_final_ana<-subset(BDD_finale_publi, select=c(Nom_scientifique,Genre,Espece, Milieu_recolte,ID_espece,ID_echantillon,ID_Feuille,T_ambiante,Vent,Humidite,DI,DI_test,score_DI,BT,BT_test,MT,BB,BB_test,FI,Nb_ramifications,BD,TMC_t0,TMC_t24,TDMC,TD,TDIA,Gmin,LMC_t0,LMC_t24,PEF,LDMC,Surface_F,SLA,LT,VPD))
head(BDD_final_ana)
dim(BDD_final_ana)
#export de la BDD 
write.csv2(BDD_final_ana,"Data/BDD_final_ana.csv")










############# Base à l'échelle de l'échantillon (moyenne des données feuille pour chaque échantillon) #######################

#création de table avec moyenne et sd pour chaque variable en fonction du nom de l'espèce
temp<-BDD_final_ana[,8:35] ###sélection des colonnes comprenant les variables pour les intégrer dans la boucle
temp

#création d'une bdd d'origine pour moyenne (sert pour merge)
BDD_moy_ech <- aggregate(temp[,1] ~ ID_echantillon + Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_final_ana, FUN = mean, na.rm = TRUE)
BDD_moy_ech[,7] <- round(BDD_moy_ech[,7], 2)
colnames(BDD_moy_ech)[7] <- colnames(temp)[1]

#création d'un bdd d'origine pour sd (sert pour merge)
BDD_sd_ech <- aggregate(temp[,1] ~ ID_echantillon + Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_final_ana, FUN = sd, na.rm = TRUE)
BDD_sd_ech[,7] <- round(BDD_sd_ech[,7], 2)
colnames(BDD_sd_ech)[7] <- colnames(temp)[1]


# Boucle pour les calcul des moyennes et écart-types
for (i in 2:ncol(temp)) {
  # Moyenne
  temp_moy <- aggregate(temp[, i] ~ ID_echantillon + Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_final_ana, FUN = mean, na.rm = TRUE)
  temp_moy[,7] <- round(temp_moy[,7], 2)
  colnames(temp_moy)[7] <- colnames(temp)[i]
  BDD_moy_ech <- merge(BDD_moy_ech, temp_moy, by = c("ID_echantillon", "Nom_scientifique", "Genre", "Espece", "ID_espece", "Milieu_recolte"), all = TRUE)
  
  # Écart-type
  temp_sd <- aggregate(temp[, i] ~ ID_echantillon + Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_final_ana, FUN = sd, na.rm = TRUE)
  temp_sd[,7] <- round(temp_sd[,7], 2)
  colnames(temp_sd)[7] <- colnames(temp)[i]
  BDD_sd_ech <- merge(BDD_sd_ech, temp_sd, by = c("ID_echantillon", "Nom_scientifique", "Genre", "Espece", "ID_espece", "Milieu_recolte"), all = TRUE)
}

# Vecteur des ID à supprimer (échantillon qui n'étaient pas bien placés au dessus du chalumeau)
ids_a_supprimer <- c("02_5", "13_6", "30_5", "39_6", "42_5")

# Trouver leslignes des ids à supprimer
lignes_a_supprimer <- c()
for (id in ids_a_supprimer) {lignes_a_supprimer <- c(lignes_a_supprimer, which(BDD_moy_ech$ID_echantillon == id))}

# Supprimer ces lignes dans l'objet d'origine
BDD_moy_ech <- BDD_moy_ech[-lignes_a_supprimer, ]


## complete les NA de BD
head(BDD_moy_ech)
BDD_moy_ech$BD_mean <- BDD_moy_ech$BD           # copier la colonne originale
BDD_moy_ech$BD_mean[is.na(BDD_moy_ech$BD_mean)] <- mean(BDD_moy_ech$BD, na.rm = TRUE)
head(BDD_moy_ech)

# Export des données
write.csv2(BDD_moy_ech, "Data/Publi/BDD_moy_ech_publi.csv", row.names = FALSE)
write.csv2(BDD_sd_ech, "Data/Publi/BDD_sd_ech_publi.csv", row.names = FALSE)




####### suppression des échantillons qui ont une valeur par défaut pour MT (pour analyser MT)
BDD_moy_echMT <- BDD_moy_ech[BDD_moy_ech$MT != 150, ]

####### suppression des échantillons qui ont une valeur par défaut pour reste composates (pour analyser DI, BB et BT)
BDD_moy_ech3 <- BDD_moy_ech[BDD_moy_ech$DI_test != 10, ]

# Export des données
write.csv2(BDD_moy_echMT, "Data/Publi/BDD_moy_echMT.csv", row.names = FALSE)
write.csv2(BDD_moy_ech3, "Data/Publi/BDD_moy_ech3.csv", row.names = FALSE)
dim(BDD_moy_ech)










############# Base à l'échelle de l'espèce (complète : sans suppression des valeurs par défaut) ######################################
#création de table avec moyenne et sd pour chaque variable en fonction du nom de l'espèce
tem1<-BDD_moy_ech[,7:35] ###sélection des colonnes comprenant les variables pour les intégrer dans la boucle
tem1
#création d'un bdd d'origine pour moyenne (sert pour merge)
BDD_moy_esp <- aggregate(tem1[,1] ~ Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_moy_ech, FUN = mean, na.rm = TRUE)
BDD_moy_esp[,6] <- round(BDD_moy_esp[,6], 2)
colnames(BDD_moy_esp)[6] <- colnames(tem1)[1]

#création d'un bdd d'origine pour sd (sert pour merge)
BDD_sd_esp <- aggregate(tem1[,1] ~ Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_moy_ech, FUN = sd, na.rm = TRUE)
BDD_sd_esp[,6] <- round(BDD_sd_esp[,6], 2)
colnames(BDD_sd_esp)[6] <- colnames(tem1)[1]

# Boucle
for (i in 2:ncol(tem1)) {
  
  # Moyenne
  tem1_moy_esp <- aggregate(tem1[, i] ~ Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_moy_ech, FUN = mean, na.rm = TRUE)
  tem1_moy_esp[,6] <- round(tem1_moy_esp[,6], 2)
  colnames(tem1_moy_esp)[6] <- colnames(tem1)[i]
  BDD_moy_esp <- merge(BDD_moy_esp, tem1_moy_esp, by = c("Nom_scientifique", "Genre", "Espece", "ID_espece", "Milieu_recolte"), all = TRUE)
  
  # SD
  tem1_sd_esp <- aggregate(tem1[, i] ~ Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_moy_ech, FUN = sd, na.rm = TRUE)
  tem1_sd_esp[,6] <- round(tem1_sd_esp[,6], 2)
  colnames(tem1_sd_esp)[6] <- colnames(tem1)[i]
  BDD_sd_esp <- merge(BDD_sd_esp, tem1_sd_esp, by = c("Nom_scientifique", "Genre", "Espece", "ID_espece", "Milieu_recolte"), all = TRUE)
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



#export de la BDD 
write.csv2(BDD_moy_esp,"Data/Publi/BDD_moy_esp.csv")
write.csv2(BDD_sd_esp, "Data/Publi/BDD_sd_esp.csv")




############# Base à l'échelle de l'espèce (nettoyée MT) ######################################
#création de table avec moyenne et sd pour chaque variable en fonction du nom de l'espèce
tem2<-BDD_moy_echMT[,7:35] ###sélection des colonnes comprenant les variables pour les intégrer dans la boucle
tem2
#création d'un bdd d'origine pour moyenne (sert pour merge)
BDD_moy_espMT <- aggregate(tem2[,1] ~ Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_moy_echMT, FUN = mean, na.rm = TRUE)
BDD_moy_espMT[,6] <- round(BDD_moy_espMT[,6], 2)
colnames(BDD_moy_espMT)[6] <- colnames(tem2)[1]

#création d'un bdd d'origine pour sd (sert pour merge)
BDD_sd_espMT <- aggregate(tem2[,1] ~ Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_moy_echMT, FUN = sd, na.rm = TRUE)
BDD_sd_espMT[,6] <- round(BDD_sd_espMT[,6], 2)
colnames(BDD_sd_espMT)[6] <- colnames(tem2)[1]


#Boucle pour les calcul des moyennes et écart-types
for (i in 2:ncol(tem2)) {
  tem2_moy_esp <- aggregate(tem2[, i] ~ Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_moy_echMT, FUN = mean, na.rm = TRUE)
  tem2_moy_esp[,6] <- round(tem2_moy_esp[,6], 2)
  colnames(tem2_moy_esp)[6] <- colnames(tem2)[i]
  BDD_moy_espMT <- merge(BDD_moy_espMT, tem2_moy_esp, by = c("Nom_scientifique", "Genre", "Espece", "ID_espece", "Milieu_recolte"), all = TRUE)
  
  tem2_sd_esp <- aggregate(tem2[, i] ~ Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_moy_echMT, FUN = sd, na.rm = TRUE)
  tem2_sd_esp[,6] <- round(tem2_sd_esp[,6], 2)
  colnames(tem2_sd_esp)[6] <- colnames(tem2)[i]
  BDD_sd_espMT <- merge(BDD_sd_espMT, tem2_sd_esp, by = c("Nom_scientifique", "Genre", "Espece", "ID_espece", "Milieu_recolte"), all = TRUE)
}



BDD_FI_esp <- aggregate(FI ~ Nom_scientifique + ID_espece + Milieu_recolte, 
                        data = BDD_moy_echMT, 
                        FUN = sum, 
                        na.rm = TRUE)
colnames(BDD_FI_esp)[4] <- "Nb_FI"

BDD_FI_esp


BDD_moy_espMT <- merge(BDD_moy_espMT, BDD_FI_esp, 
                      by = c("Nom_scientifique", "ID_espece", "Milieu_recolte"), 
                      all.x = TRUE)

#Résultat
head(BDD_moy_espMT)
View(BDD_moy_espMT)

#export de la BDD 
write.csv2(BDD_moy_espMT,"Data/Publi/BDD_moy_espMT.csv")



############# Base à l'échelle de l'espèce (nettoyée reste) ######################################
tem2 <- BDD_moy_ech3[, 7:35]

# Moyenne
BDD_moy_esp3 <- aggregate(tem2[,1] ~ Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_moy_ech3, FUN = mean, na.rm = TRUE)
BDD_moy_esp3[,6] <- round(BDD_moy_esp3[,6], 2)
colnames(BDD_moy_esp3)[6] <- colnames(tem2)[1]

# SD
BDD_sd_esp3 <- aggregate(tem2[,1] ~ Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_moy_ech3, FUN = sd, na.rm = TRUE)
BDD_sd_esp3[,6] <- round(BDD_sd_esp3[,6], 2)
colnames(BDD_sd_esp3)[6] <- colnames(tem2)[1]

# Boucle
for (i in 2:ncol(tem2)) {
  
  tem2_moy_esp <- aggregate(tem2[, i] ~ Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_moy_ech3, FUN = mean, na.rm = TRUE)
  tem2_moy_esp[,6] <- round(tem2_moy_esp[,6], 2)
  colnames(tem2_moy_esp)[6] <- colnames(tem2)[i]
  BDD_moy_esp3 <- merge(BDD_moy_esp3, tem2_moy_esp, by = c("Nom_scientifique", "Genre", "Espece", "ID_espece", "Milieu_recolte"), all = TRUE)
  
  tem2_sd_esp <- aggregate(tem2[, i] ~ Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_moy_ech3, FUN = sd, na.rm = TRUE)
  tem2_sd_esp[,6] <- round(tem2_sd_esp[,6], 2)
  colnames(tem2_sd_esp)[6] <- colnames(tem2)[i]
  BDD_sd_esp3 <- merge(BDD_sd_esp3, tem2_sd_esp, by = c("Nom_scientifique", "Genre", "Espece", "ID_espece", "Milieu_recolte"), all = TRUE)
}

BDD_FI_esp <- aggregate(FI ~ Nom_scientifique + ID_espece + Milieu_recolte, 
                        data = BDD_moy_ech3, 
                        FUN = sum, 
                        na.rm = TRUE)
colnames(BDD_FI_esp)[4] <- "Nb_FI"

BDD_FI_esp


BDD_moy_esp3 <- merge(BDD_moy_esp3, BDD_FI_esp, 
                      by = c("Nom_scientifique", "ID_espece", "Milieu_recolte"), 
                      all.x = TRUE)

#Résultat
head(BDD_moy_esp3)
View(BDD_moy_esp3)

#export de la BDD 
write.csv2(BDD_moy_esp3,"Data/Publi/BDD_moy_esp3.csv")
















