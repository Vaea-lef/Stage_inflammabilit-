
#### TRAITEMENT DONNEE STAGE INFLAMMABILITE ####


# Installation des packages si besoin (pour l'import et l'export en format excel directement)
install.packages("readxl")
install.packages("writexl")

# Chargement des packages
library(readxl)     # pour lire les fichiers Excel
library(writexl)    # pour écrire un fichier Excel





############# Traitement des données feuille ###########################


##  importation BDD_feuille en format CSV
setwd("C:/Users/vaeal/OneDrive/VAEA_IRD/Data/BDD_brutes/CSV") #définition du répertoire de travail
BDD_leaf<-read.csv2("BDD_Feuille.csv", header = TRUE) #importation de la base

##  importation BDD_feuille excel
setwd("C:/Users/vaeal/OneDrive/VAEA_IRD/Data/BDD_brutes")
BDD_leaf <- read_excel("BDD_Feuille.xlsx")

# pour éventuellement visualiser la base importée
View(BDD_leaf) 

## Manipulation des données : calcul de Leaf Moisture Content (LMC) à t0 et t24,Leaf Dry Matter Content (LDMC) et Specific leaf area (SLA) 

#calcul contenu en eau "avant brûlage" (LMC_t24) 
LMC_t24<-((BDD_leaf$Masse_F_T24-BDD_leaf$Masse_F_Tsec)/BDD_leaf$Masse_F_Tsec)*100
LMC_t24<-round(LMC_t24,2)
LMC_t24

#calcul contenu en eau "frais" (LMC_t0)
LMC_t0<-((BDD_leaf$Masse_F_T0-BDD_leaf$Masse_F_Tsec)/BDD_leaf$Masse_F_Tsec)*100
LMC_t0<-round(LMC_t0,2)
LMC_t0

#calcul matière sèche (LDMC) (mg/g)
LDMC<-((BDD_leaf$Masse_F_Tsec * 1000)/BDD_leaf$Masse_F_T0)
LDMC<-round(LDMC,2)
LDMC

#calcul surface foliaire spécifique (SLA) (mm²/g)
SLA<-((BDD_leaf$Surface_F*100)/BDD_leaf$Masse_F_Tsec)
SLA<-round(SLA,2)
SLA

#Leaf Thickness = épaisseur de la feuille (mm)
Leaf_ep<-((1/SLA)*LDMC)
Leaf_ep<-round(Leaf_ep,4)
Leaf_ep

# création d'une nouvelle base de données "feuille" calculée avec ajout des colonnes
BDD_leaf_calcule<-data.frame(BDD_leaf,LMC_t24,LMC_t0,LDMC,SLA,Leaf_ep)
View(BDD_leaf_calcule) #pour voir la BDD finale

###### Export de la BDD excel et CSV

setwd("C:/Users/vaeal/OneDrive/VAEA_IRD/Data/BDD_ traitees")  ## changement du répertoire de travail :choix de la "zone" souhaitée d'enregistrement des fichiers
# Exporter les données avec les colonnes calculées EXCEL
write_xlsx(BDD_leaf_calcule, "BDD_leaf_calcule.xlsx")

setwd("C:/Users/vaeal/OneDrive/VAEA_IRD/Data/BDD_ traitees/CSV_calcule")
# Exporter les données avec les colonnes calculées CSV
write.csv2(BDD_leaf_calcule,"BDD_leaf_calcule.csv")






###################  traitement des données de Tige  #############################

##  importation BDD_traits CSV
setwd("C:/Users/vaeal/OneDrive/VAEA_IRD/Data/BDD_brutes/CSV")
BDD_traits<-read.csv2("BDD_Traits.csv")

# pour éventuellement visualiser la base importée
View(BDD_traits) 

## Manipulation des données : calcul du Volume et de la densité de la tige ainsi que Twig Moisture Content (TMC) et Twig Dry Matter Content (TDMC)

#calcul volume de la tige
Volume_tige<-(pi*((((BDD_traits$diametre_T_1_mm + BDD_traits$diametre_T_2_mm)/2)/2)/10)* BDD_traits$Longeur_T_cm)
Volume_tige<-round(Volume_tige,2)
Volume_tige

#calcul densitée de la tige (g/cm3)
Densite_tige<-(BDD_traits$Masse_T_T0/Volume)
Densite_tige<-round(Densite_tige,2)
Densite_tige

#calcul contenu en eau tige "avant brûlage" (TMC)
TMC_t24<-((BDD_leaf$Masse_F_T24-BDD_leaf$Masse_F_Tsec)/BDD_leaf$Masse_F_Tsec)*100
TMC_t24<-round(TMC_t24,2)
TMC_t24

#calcul TMC frais
TMC_t0<-((BDD_leaf$Masse_F_T0-BDD_leaf$Masse_F_Tsec)/BDD_leaf$Masse_F_Tsec)*100
TMC_t0<-round(TMC_t0,2)
TMC_t0

#calcul matière sèche tige (TDMC) (mg/g)
TDMC<-((BDD_leaf$Masse_F_Tsec * 1000)/BDD_leaf$Masse_F_T0)
TDMC<-round(TDMC,2)
TDMC

# création d'une nouvelle base de données calculée avec ajout des colonnes
BDD_traits_calcule<-data.frame(BDD_traits,Volume_tige,Densite_tige)
View(BDD_traits_calcule) #pour voir la BDD finale

###### Export de la BDD excel et CSV

setwd("C:/Users/vaeal/OneDrive/VAEA_IRD/Data/BDD_ traitees")
# Exporter les données avec les colonnes calculées EXCEL
write_xlsx(BDD_traits_calcule, "BDD_traits_calcule.xlsx")

setwd("C:/Users/vaeal/OneDrive/VAEA_IRD/Data/BDD_ traitees/CSV_calcule")
# Exporter les données avec les colonnes calculées CSV
write.csv2(BDD_traits_calcule,"BDD_traits_calcule.csv")




############# Traitement des données Inflammabilité ###########################

setwd("C:/Users/vaeal/OneDrive/VAEA_IRD/Data/BDD_brutes/CSV")
BDD_infla<-read.csv2("BDD_Inflammabilite.csv")

#calcul du volume des échantillons (m3)
Volume_Echantillon<-((2/3)*pi*BDD_infla$longeur*((BDD_infla$largeur/2)*(BDD_infla$hauteur/2)))/1000000
Volume_Echantillon<-round(Volume_Echantillon,3)
Volume_Echantillon

#calcul de la Bulk density (kg/m3)
Bulk<-((BDD_infla$masse/1000)/Volume_Echantillon)
Bulk<-round(Bulk,2)
Bulk

#calcul de Burning Time (BT) (s) temps de combustion
BT<-(BDD_infla$temps_total - (120 + BDD_infla$DI))
BT<-round(BT)
BT

# création d'une nouvelle base de données calculée avec ajout des colonnes
BDD_infla_calcule<-data.frame(BDD_infla,Volume_Echantillon,Bulk,BT)
View(BDD_infla_calcule) #pour voir la BDD finale

#Export de la base calculée

setwd("C:/Users/vaeal/OneDrive/VAEA_IRD/Data/BDD_ traitees")
# Exporter les données avec les colonnes calculées EXCEL
write_xlsx(BDD_infla_calcule, "BDD_infla_calcule.xlsx")

setwd("C:/Users/vaeal/OneDrive/VAEA_IRD/Data/BDD_ traitees/CSV_calcule")
# Exporter les données avec les colonnes calculées CSV
write.csv2(BDD_infla_calcule,"BDD_infla_calcule.csv")



######### assembler les bases de données calculée en une seule base #########

#assemblage : BDD traits et leaf
data1<-merge(BDD_leaf_calcule,BDD_traits_calcule)
View(data1)
data1<-data.frame(data1)

#assemblage data1 avec BDD inflammabilité calculée
data2<-merge(data1,BDD_infla_calcule)
View(data2)


#importation de la base échantillonnage

setwd("C:/Users/vaeal/OneDrive/VAEA_IRD/Data/BDD_brutes/CSV")
BDD_ech<-read.csv2("BDD_Echantillonnage.csv")

#assemblage data2 avec BDD Echantillonnage
BDD_finale<-merge(data2,BDD_ech)
View(BDD_finale)

#export de la BDD finale 

setwd("C:/Users/vaeal/OneDrive/VAEA_IRD/Data/BDD_ traitees")
# Exporter les données avec les colonnes calculées EXCEL
write_xlsx(BDD_finale, "BDD_finale.xlsx")

setwd("C:/Users/vaeal/OneDrive/VAEA_IRD/Data/BDD_ traitees/CSV_calcule")
# Exporter les données avec les colonnes calculées CSV
write.csv2(BDD_finale,"BDD_finale.csv")
