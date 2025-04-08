
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

# création d'une nouvelle base de données "feuille" calculée avec ajout des colonnes calculées
BDD_leaf_calcule<-data.frame(BDD_leaf,LMC_t0,LMC_t24,LDMC,SLA,LT)
View(BDD_leaf_calcule) #pour voir la BDD finale

###### Export de la BDD_Leaf_F CSV
write.csv2(BDD_leaf_F,"Data/BDD_leaf_calcule.csv")






###################  traitement des données de Tige  #############################

##  importation BDD_traits CSV
BDD_traits<-read.csv2("Data/BDD_Traits.csv")

# pour éventuellement visualiser la base importée
View(BDD_traits) 

## Manipulation des données : calcul du Volume et de la densité de la tige ainsi que Twig Moisture Content (TMC) et Twig Dry Matter Content (TDMC)

#calcul volume de la tige
V_tige<-(pi*((((BDD_traits$diametre_T_1_mm + BDD_traits$diametre_T_2_mm)/2)/2)/10)* BDD_traits$Longeur_T_cm)
V_tige<-round(V_tige,2)
V_tige

#calcul densitée de la tige (g/cm3)
D_tige<-(BDD_traits$Masse_T_T0/V_tige)
D_tige<-round(D_tige,2)
D_tige

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

# création d'une nouvelle base de données calculée avec ajout des colonnes
BDD_traits_calcule<-data.frame(BDD_traits,V_tige,D_tige,TMC_t0,TMC_t24,TDMC)
View(BDD_traits_calcule) #pour voir la BDD finale

###### Export de la BDD traits calculée CSV

# Exporter les données avec les colonnes calculées CSV
write.csv2(BDD_traits_calcule,"Data/BDD_traits_calcule.csv")




############# Traitement des données Inflammabilité ###########################


BDD_infla<-read.csv2("Data/BDD_Inflammabilite.csv")
View(BDD_infla)

#calcul du volume des échantillons (m3)
V_Ech<-((2/3)*pi*BDD_infla$longeur*((BDD_infla$largeur/2)*(BDD_infla$hauteur/2)))/1000000
V_Ech<-round(V_Ech,3)
V_Ech

#calcul de la Bulk density (kg/m3)
Bulk<-((BDD_infla$masse/1000)/V_Ech)
Bulk<-round(Bulk,2)
Bulk

#calcul de Burning Time (BT) (s) temps de combustion
BT<-(BDD_infla$temps_total - (120 + BDD_infla$DI))
BT<-round(BT)
BT

# création d'une nouvelle base de données calculée avec ajout des colonnes
BDD_infla_calcule<-data.frame(BDD_infla,V_Ech,Bulk,BT)
View(BDD_infla_calcule) #pour voir la BDD finale

###### Export de la base infla calculée CSV
# Exporter les données avec les colonnes calculées CSV
write.csv2(BDD_infla_calcule,"Data/BDD_infla_calcule.csv")





######### Assembler les bases de données calculée en une seule base #########

#assemblage : BDD traits et leaf
data1<-merge(BDD_leaf_calcule,BDD_traits_calcule)
View(data1)
data1<-data.frame(data1)

#assemblage data1 avec BDD inflammabilité calculée
data2<-merge(data1,BDD_infla_calcule)
View(data2)


#importation de la base échantillonnage
BDD_ech<-read.csv2("Data/BDD_Echantillonnage.csv")

#assemblage data2 avec BDD Echantillonnage
BDD_finale<-merge(data2,BDD_ech)
View(BDD_finale)

#export de la BDD finale 

write.csv2(BDD_finale,"Data/BDD_finale.csv")
