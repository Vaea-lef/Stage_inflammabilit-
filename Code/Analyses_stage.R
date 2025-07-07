
#### ANALYSE DE DONNEE STAGE INFLAMMABILITE ####


################# IMPORTATION ET CHARGEMENT PACKAGE ###########################

##  importation BDD_ana_ech en format CSV
setwd("C:/IRD/Stage_inflammabilit-") #définition du répertoire de travail
BDD_ech<-read.csv2("Data/BDD_moy_ech.csv", header = TRUE) #importation de l
BDD_ech

##  importation BDD_sd_esp en format CSV
BDD_sd_esp<-read.csv2("Data/BDD_sd_esp.csv", header = TRUE)
BDD_sd_esp

##  importation BDD_moy_esp en format CSV
BDD_esp<-read.csv2("Data/BDD_moy_esp.csv", header = TRUE) #importation de la base
BDD_esp
names(BDD_esp)[which(names(BDD_esp) == "Nb_ramifications")] <- "Nb_rami"
dim(BDD_esp)

##  Charger les packages nécessaires pour ACP
library(FactoMineR)
library(factoextra)





################  ACP INFLAMMABILITE ####################

# Sélection des colonnes des composantes de l'inflammabilité
colonnes_infla <- BDD_esp[, c(5:8)]

# Vérification des données
head(colonnes_infla)

# Application de l'ACP
res.pca <- PCA(colonnes_infla, scale.unit = TRUE, graph = FALSE)

# Résumé des résultats
summary(res.pca)

# Graphique des contributions des variables aux composantes principales
fviz_pca_var(res.pca, col.var = "contrib", gradient.cols = c("#6fec00", "#ff9e00", "#Ff0000"),repel = TRUE)

# Graphique des individus 
fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)

# Graphique combiné des variables et des individus
fviz_pca_biplot(res.pca,col.var = "contrib", gradient.cols =c("#00AFBB", "#E7B800", "#FC4E07") , repel = TRUE)

# Afficher l'ébouli
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50), main="Graphique de l'ébouli")








# Ajout d'un score d'inflammabilité basé sur la coordonnée de l'axe 1 et moyenne pondérée axe 1 et 2
coord <- res.pca$ind$coord  # coordonnées des individus
BDD_esp$score <- coord[, 1]  # Dim 1 = axe 1
BDD_esp$score2 <- (0.8*coord[,1] + 0.2*coord[,2]) / 2

# normalisation du score d'inflammabilité entre -1 et 1
min_score <- min(BDD_esp$score2, na.rm = TRUE)
max_score <- max(BDD_esp$score2, na.rm = TRUE)
BDD_esp$score_normalise <- -1 + (BDD_esp$score2 - min_score) * 2 / (max_score - min_score)

# visualiser les scores
View(BDD_esp)





#################### ACP TRAITS ############################

# Sélection des colonnes des traits fonctionnels
colonnes_traits <- na.omit(BDD_esp[,c(9:24)])

# Vérifier les données
colonnes_traits

# Application de l'ACP
res.pca <- PCA(colonnes_traits, scale.unit = TRUE, graph = FALSE)

# Résumé des résultats
summary(res.pca)

# Graphique des contributions des variables aux composantes principales
fviz_pca_var(res.pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

# Graphique des individus 
fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)

# Graphique combiné des variables et des individus
fviz_pca_biplot(res.pca, col.var = "contrib", gradient.cols =c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)

# Afficher l'ébouli
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50), main="Graphique de l'ébouli")






######## ACP INFLA avec projection des axes TRAITS ####################

# Sélection des lignes complètes (pas de NA dans les colonnes 3 à 22)
colonnes_complet <- na.omit(BDD_esp[, 5:24])

## ACP ##
# colonnes 1 à 4 : inflammabilité (axes actifs)
# colonnes 5 à 20 : traits (axes supplémentaires)
res.pca <- PCA(colonnes_complet, scale.unit = TRUE, 
               quanti.sup = 5:20, graph = FALSE)

# Visualisation ACP
fviz_pca_var(res.pca, col.var = "red", repel = TRUE) +
  ggtitle("Variables actives (infla) et supplémentaires (traits)")

# ACP avec les individus
fviz_pca_biplot(res.pca, col.var = "red",repel = TRUE)





################ BOXPLOT INFLA ESPECE ########################

######## DI ###########
# Calcul des médianes par espèce
med_DI <- aggregate(DI ~ Nom_scientifique, data = BDD_ech, median)

# Création ordre décroissant des médianes
ordre <- med_DI$Nom_scientifique[order(med_DI$DI, decreasing = TRUE)]

# facteur avec ce nouvel ordre (pour que l'ordre soit capté par le boxplot)
BDD_ech$Nom_scientifique <- factor(BDD_ech$Nom_scientifique, levels = ordre)

# Traçage du boxplot avec l'ordre décroissant 
boxplot(DI ~ Nom_scientifique, data = BDD_ech,
        las = 2, cex.axis = 0.7,
        main = "DI par espèce (ordre décroissant de la médiane)",
        xlab = "Espèce", ylab = "DI")



######## BT ###########
# Calcul des médianes par espèce
med_BT <- aggregate(BT ~ Nom_scientifique, data = BDD_ech, median)

# Création ordre décroissant des médianes
ordre <- med_BT$Nom_scientifique[order(med_BT$BT, decreasing = TRUE)]

# facteur avec ce nouvel ordre (pour que l'ordre soit capté par le boxplot)
BDD_ech$Nom_scientifique <- factor(BDD_ech$Nom_scientifique, levels = ordre)

# Traçage du boxplot avec l'ordre décroissant 
boxplot(BT ~ Nom_scientifique, data = BDD_ech,
        las = 2, cex.axis = 0.7,
        main = "BT par espèce (ordre décroissant de la médiane)",
        xlab = "Espèce", ylab = "BT")

######## MT ###########
# Calcul des médianes par espèce
med_MT <- aggregate(MT ~ Nom_scientifique, data = BDD_ech, median)

# Création ordre décroissant des médianes
ordre <- med_MT$Nom_scientifique[order(med_MT$MT, decreasing = TRUE)]

# facteur avec ce nouvel ordre (pour que l'ordre soit capté par le boxplot)
BDD_ech$Nom_scientifique <- factor(BDD_ech$Nom_scientifique, levels = ordre)

# Traçage du boxplot avec l'ordre décroissant 
boxplot(MT ~ Nom_scientifique, data = BDD_ech,
        las = 2, cex.axis = 0.7,
        main = "MT par espèce (ordre décroissant de la médiane)",
        xlab = "Espèce", ylab = "MT")


######## BB ###########
# Calcul des médianes par espèce
med_BB <- aggregate(BB ~ Nom_scientifique, data = BDD_ech, median)

# Création ordre décroissant des médianes
ordre <- med_BB$Nom_scientifique[order(med_BB$BB, decreasing = TRUE)]

# facteur avec ce nouvel ordre (pour que l'ordre soit capté par le boxplot)
BDD_ech$Nom_scientifique <- factor(BDD_ech$Nom_scientifique, levels = ordre)

# Traçage du boxplot avec l'ordre décroissant 
boxplot(BB ~ Nom_scientifique, data = BDD_ech,
        las = 2, cex.axis = 0.7,
        main = "BB par espèce (ordre décroissant de la médiane)",
        xlab = "Espèce", ylab = "BB")






################# Heatmap (couleur en fonction de chaque composante d'inflammabilité) ##################

# packages nécessaires
library(ggplot2)
library(tidyr)  # pour changer le format de la BDD en format "long"

# Copier la BDD (pour ne pas la modifier)
df_prep <- BDD_esp

# Normalisser les valeur entre 0 et 1 : fonction
normalize <- function(x) {  return((x - min(x)) / (max(x) - min(x)))}

df_norm <- df_prep
df_norm$MT <- normalize(df_prep$MT)
df_norm$DI <- normalize(df_prep$DI)
df_norm$BB <- normalize(df_prep$BB)
df_norm$BT <- normalize(df_prep$BT)

# Transformer en format long (une ligne par composante)
df_long <- pivot_longer(df_norm,
                        cols = c("MT", "DI", "BB", "BT"),
                        names_to = "Variable",
                        values_to = "Valeur")

# Créer l’ordre des espèces en fonction du score d'inflammabilité
ordre_esp <- df_prep[order(df_prep$score_normalise, decreasing = FALSE), "Nom_scientifique"]

# enlever les doublons et les NA
ordre_esp <- unique(na.omit(ordre_esp))

# Appliquer le facteur pour ordonner les espèces
df_long$Nom_scientifique <- factor(df_long$Nom_scientifique, levels = rev(ordre_esp))

# Générer la heatmap
ggplot(df_long, aes(x = Variable, y = Nom_scientifique, fill = Valeur)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(
    colors = c("darkgreen", "yellow", "red"),
    breaks = c(0.1, 0.5, 0.9),
    labels = c("Faible", "Moyenne", "Élevée"),
    name = "Inflammabilité"
  ) +
  labs(title = "Heatmap des composantes d’inflammabilité",
       x = "Composantes",
       y = "Espèces") +
  scale_y_discrete(drop = FALSE) +  # Affiche toutes les espèces, même si certaines valeurs sont NA
  theme(axis.text.y = element_text(size = 6))  # Réduction de la taille de texte si besoin





####################### matrice de corrélation (pour le choix des traits) ######################

### des traits foncitonnels
library (corrplot)
mat_cor_trait<-cor(colonnes_traits)
corrplot(mat_cor_trait)
round(mat_cor_trait, 2)   ## afficher les valeurs
# Visualiser la corrélation
corrplot(mat_cor_trait, method = "color", type = "upper", tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")
corrplot(mat_cor_trait, method = "color", tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")






#################### TEST ANOVA POUR VERIFICATION VARIABILITE ESPECES ######################
# inflammabilité
  #DI
variabilité_esp <- aov(BDD_ech$DI~BDD_ech$Nom_scientifique)
summary(variabilité_esp)
  #BT
variabilité_esp <- aov(BDD_ech$BT~BDD_ech$Nom_scientifique)
summary(variabilité_esp)
  #MT
variabilité_esp <- aov(BDD_ech$MT~BDD_ech$Nom_scientifique)
summary(variabilité_esp)
  #BB
variabilité_esp <- aov(BDD_ech$MT~BDD_ech$Nom_scientifique)
summary(variabilité_esp)

# traits
  #Nb_rami
variabilité_esp <- aov(BDD_ech$Nb_rami~BDD_ech$Nom_scientifique)
summary(variabilité_esp)
  #SV
variabilité_esp <- aov(BDD_ech$SV~BDD_ech$Nom_scientifique)
summary(variabilité_esp)
  #SD
variabilité_esp <- aov(BDD_ech$SD~BDD_ech$Nom_scientifique)
summary(variabilité_esp)
  #TDIA
variabilité_esp <- aov(BDD_ech$TDIA~BDD_ech$Nom_scientifique)
summary(variabilité_esp)
  #TD
variabilité_esp <- aov(BDD_ech$TD~BDD_ech$Nom_scientifique)
summary(variabilité_esp)
  #Gmin
variabilité_esp <- aov(BDD_ech$Gmin~BDD_ech$Nom_scientifique)
summary(variabilité_esp)
  #LMC_t0
variabilité_esp <- aov(BDD_ech$LMC_t0~BDD_ech$Nom_scientifique)
summary(variabilité_esp)
  #Surface_F
variabilité_esp <- aov(BDD_ech$Surface_F~BDD_ech$Nom_scientifique)
summary(variabilité_esp)
  #SLA
variabilité_esp <- aov(BDD_ech$SLA~BDD_ech$Nom_scientifique)
summary(variabilité_esp)









############# Proportions of variance explained (Alam) #############
library(lme4) #Package pour modèle linéaire mixte (random intercept model)

# BT
  # Modèle linéaire mixte
m_BT <- lmer(BT ~ 1 + (1 | Nom_scientifique), data = BDD_ech)
  #Analyse des résidus 
plot(m_BT)
hist(resid(m_BT)) 
  # Extraire les composantes de variance
var_BT <- as.data.frame(VarCorr(m_BT))
var_BT_inter <- var_BT$vcov[1]
var_BT_intra <- attr(VarCorr(m_BT), "sc")^2
  # Pourcentage de variance
total_BT <- var_BT_inter + var_BT_intra
100 * var_BT_inter / total_BT  # variance inter (espèce)
100 * var_BT_intra / total_BT  # variance intra (résiduelle)

# MT
m_MT <- lmer(MT ~ 1 + (1 | Nom_scientifique), data = BDD_ech)
plot(m_MT)            
hist(resid(m_MT))
     
var_MT <- as.data.frame(VarCorr(m_MT))
var_MT_inter <- var_MT$vcov[1]
var_MT_intra <- attr(VarCorr(m_MT), "sc")^2

total_MT <- var_MT_inter + var_MT_intra
100 * var_MT_inter / total_MT  
100 * var_MT_intra / total_MT  

# DI
m_DI <- lmer(DI ~ 1 + (1 | Nom_scientifique), data = BDD_ech)
plot(m_DI)            
hist(resid(m_DI))

var_DI <- as.data.frame(VarCorr(m_DI))
var_DI_inter <- var_DI$vcov[1]
var_DI_intra <- attr(VarCorr(m_DI), "sc")^2

total_DI <- var_DI_inter + var_DI_intra
100 * var_DI_inter / total_DI  
100 * var_DI_intra / total_DI 

# BB
m_BB <- lmer(BB ~ 1 + (1 | Nom_scientifique), data = BDD_ech)
plot(m_BB)            
hist(resid(m_BB))

var_BB <- as.data.frame(VarCorr(m_BB))
var_BB_inter <- var_BB$vcov[1]
var_BB_intra <- attr(VarCorr(m_BB), "sc")^2

total_BB <- var_BB_inter + var_BB_intra
100 * var_BB_inter / total_BB  
100 * var_BB_intra / total_BB 

# Nb_ramifications
m_Nb_ramifications <- lmer(Nb_ramifications ~ 1 + (1 | Nom_scientifique), data = BDD_ech)

var_Nb_ramifications <- as.data.frame(VarCorr(m_Nb_ramifications))
var_Nb_ramifications_inter <- var_Nb_ramifications$vcov[1]
var_Nb_ramifications_intra <- attr(VarCorr(m_Nb_ramifications), "sc")^2

total_Nb_ramifications <- var_Nb_ramifications_inter + var_Nb_ramifications_intra
100 * var_Nb_ramifications_inter / total_Nb_ramifications  
100 * var_Nb_ramifications_intra / total_Nb_ramifications 

# SD
m_SD <- lmer(SD ~ 1 + (1 | Nom_scientifique), data = BDD_ech)

var_SD <- as.data.frame(VarCorr(m_SD))
var_SD_inter <- var_SD$vcov[1]
var_SD_intra <- attr(VarCorr(m_SD), "sc")^2

total_SD <- var_SD_inter + var_SD_intra
100 * var_SD_inter / total_SD  
100 * var_SD_intra / total_SD 

# TDIA
m_TDIA <- lmer(TDIA ~ 1 + (1 | Nom_scientifique), data = BDD_ech)

var_TDIA <- as.data.frame(VarCorr(m_TDIA))
var_TDIA_inter <- var_TDIA$vcov[1]
var_TDIA_intra <- attr(VarCorr(m_TDIA), "sc")^2

total_TDIA <- var_TDIA_inter + var_TDIA_intra
100 * var_TDIA_inter / total_TDIA  
100 * var_TDIA_intra / total_TDIA 

# LT
m_LT <- lmer(LT ~ 1 + (1 | Nom_scientifique), data = BDD_ech)

var_LT <- as.data.frame(VarCorr(m_LT))
var_LT_inter <- var_LT$vcov[1]
var_LT_intra <- attr(VarCorr(m_LT), "sc")^2

total_LT <- var_LT_inter + var_LT_intra
100 * var_LT_inter / total_LT  
100 * var_LT_intra / total_LT 

# LMC_t0
m_LMC_t0 <- lmer(LMC_t0 ~ 1 + (1 | Nom_scientifique), data = BDD_ech)

var_LMC_t0 <- as.data.frame(VarCorr(m_LMC_t0))
var_LMC_t0_inter <- var_LMC_t0$vcov[1]
var_LMC_t0_intra <- attr(VarCorr(m_LMC_t0), "sc")^2

total_LMC_t0 <- var_LMC_t0_inter + var_LMC_t0_intra
100 * var_LMC_t0_inter / total_LMC_t0  
100 * var_LMC_t0_intra / total_LMC_t0 

# LMC_t24
m_LMC_t24 <- lmer(LMC_t24 ~ 1 + (1 | Nom_scientifique), data = BDD_ech)

var_LMC_t24 <- as.data.frame(VarCorr(m_LMC_t24))
var_LMC_t24_inter <- var_LMC_t24$vcov[1]
var_LMC_t24_intra <- attr(VarCorr(m_LMC_t24), "sc")^2

total_LMC_t24 <- var_LMC_t24_inter + var_LMC_t24_intra
100 * var_LMC_t24_inter / total_LMC_t24  
100 * var_LMC_t24_intra / total_LMC_t24 

# PEF
m_PEF <- lmer(PEF ~ 1 + (1 | Nom_scientifique), data = BDD_ech)

var_PEF <- as.data.frame(VarCorr(m_PEF))
var_PEF_inter <- var_PEF$vcov[1]
var_PEF_intra <- attr(VarCorr(m_PEF), "sc")^2

total_PEF <- var_PEF_inter + var_PEF_intra
100 * var_PEF_inter / total_PEF  
100 * var_PEF_intra / total_PEF 

# Gmin
m_Gmin <- lmer(Gmin ~ 1 + (1 | Nom_scientifique), data = BDD_ech)

var_Gmin <- as.data.frame(VarCorr(m_Gmin))
var_Gmin_inter <- var_Gmin$vcov[1]
var_Gmin_intra <- attr(VarCorr(m_Gmin), "sc")^2

total_Gmin <- var_Gmin_inter + var_Gmin_intra
100 * var_Gmin_inter / total_Gmin  
100 * var_Gmin_intra / total_Gmin 

# SLA
m_SLA <- lmer(SLA ~ 1 + (1 | Nom_scientifique), data = BDD_ech)

var_SLA <- as.data.frame(VarCorr(m_SLA))
var_SLA_inter <- var_SLA$vcov[1]
var_SLA_intra <- attr(VarCorr(m_SLA), "sc")^2

total_SLA <- var_SLA_inter + var_SLA_intra
100 * var_SLA_inter / total_SLA  
100 * var_SLA_intra / total_SLA 

# Surface_F
m_Surface_F <- lmer(Surface_F ~ 1 + (1 | Nom_scientifique), data = BDD_ech)

var_Surface_F <- as.data.frame(VarCorr(m_Surface_F))
var_Surface_F_inter <- var_Surface_F$vcov[1]
var_Surface_F_intra <- attr(VarCorr(m_Surface_F), "sc")^2

total_Surface_F <- var_Surface_F_inter + var_Surface_F_intra
100 * var_Surface_F_inter / total_Surface_F  
100 * var_Surface_F_intra / total_Surface_F 







################### matrice de corrélation pour observer l'effet des traits sur infla ################
### de tout 
library (corrplot)
mat_cor_all<-cor(colonnes_complet)
corrplot(mat_cor_all)
round(mat_cor_all, 2)   ## afficher les valeurs
corrplot(mat_cor_all, method = "color", tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")












#################### PLOT ##################################
# histogramme pour vérifier la normalité des variables réponse (infla)
hist(BDD_ech$MT)#normal
hist(log(BDD_ech$MT))
hist(BDD_ech$BT)#necéssite peut être tranfo log 
hist(1/BDD_ech$)
hist(log(BDD_ech$BT)+1)
log(0)
hist(BDD_ech$DI)#voir 
hist(BDD_ech$BB) #voir

###### MT #########
#plot avec traits (VE)
plot(BDD_ech$LT, BDD_ech$MT, xlab = "LT", ylab = "MT")
plot(BDD_ech$LMC_t0, BDD_ech$MT, xlab = "LMC_t0", ylab = "MT")
plot(BDD_ech$LMC_t24, BDD_ech$MT, xlab = "LMC_t24", ylab = "MT")
plot(BDD_ech$Nb_ramifications, BDD_ech$MT, xlab = "Nombre de ramifications", ylab = "MT")
plot(BDD_ech$PEF, BDD_ech$MT, xlab = "PEF", ylab = "MT")
plot(BDD_ech$Gmin, BDD_ech$MT, xlab = "Gmin", ylab = "MT")
plot(BDD_ech$SLA, BDD_ech$MT, xlab = "SLA", ylab = "MT")
plot(BDD_ech$Surface_F, BDD_ech$MT, xlab = "Surface foliaire", ylab = "MT")
plot(BDD_ech$TDIA, BDD_ech$MT, xlab = "Diamètre de tige", ylab = "MT")
plot(BDD_ech$SD, BDD_ech$MT, xlab = "SD", ylab = "MT")

#modeles
m_nul=lm(MT~LT*LMC_t0*LMC_t24*Nb_ramifications*PEF*Gmin*SLA*Surface_F*TDIA*SD,data=BDD_ech)
summary(m_nul)
m1 <- lm(MT ~  LT+ LMC_t0 + Nb_ramifications  + Gmin + SLA + Surface_F  + SD, data = BDD_ech)
summary(m1)
m2 <- lm(MT ~  LT+ LMC_t0 + Nb_ramifications + SLA + Gmin  + SD, data = BDD_ech)
summary(m2)
m3 <- lm(MT ~  LT+ LMC_t0  + SLA + Gmin  + SD, data = BDD_ech)
summary(m3)
m4 <- lm(MT ~  LT+ LMC_t0  + Gmin  + SD, data = BDD_ech)
summary(m4)
m5 <- lm(MT ~   LMC_t0 + SD, data = BDD_ech)
summary(m5)

m6 <- lm(MT ~   LMC_t0 * SD, data = BDD_ech)
summary(m6)

m7 <- lm(MT ~   LMC_t0 , data = BDD_ech)
summary(m7)
m8 <- lm(MT ~   SD, data = BDD_ech)
summary(m8)

AIC(m_nul,m1,m2,m3,m4,m5,m6) 
plot(m6)
hist(resid(m6))


###### BT #########
#plot avec traits (VE)
plot(BDD_ech$LT, BDD_ech$BT, xlab = "LT", ylab = "BT")
plot(BDD_ech$LMC_t0, BDD_ech$BT, xlab = "LMC_t0", ylab = "BT")
plot(BDD_ech$LMC_t24, BDD_ech$BT, xlab = "LMC_t24", ylab = "BT")
plot(BDD_ech$Nb_ramifications, BDD_ech$BT, xlab = "Nombre de ramifications", ylab = "BT")
plot(BDD_ech$PEF, BDD_ech$BT, xlab = "PEF", ylab = "BT")
plot(BDD_ech$Gmin, BDD_ech$BT, xlab = "Gmin", ylab = "BT")
plot(BDD_ech$SLA, BDD_ech$BT, xlab = "SLA", ylab = "BT")
plot(BDD_ech$Surface_F, BDD_ech$BT, xlab = "Surface foliaire", ylab = "BT")
plot(BDD_ech$TDIA, BDD_ech$BT, xlab = "Diamètre de tige", ylab = "BT")
plot(BDD_ech$SD, BDD_ech$BT, xlab = "SD", ylab = "BT")

#modeles
m1 <- lm(BT ~  LT+ LMC_t0 + Nb_ramifications  + Gmin + SLA + Surface_F  + SD, data = BDD_ech)
summary(m1)
m2 <- lm(BT ~  LT+ LMC_t0 + Nb_ramifications + SLA + Surface_F + SD, data = BDD_ech)
summary(m2)
m3 <- lm(BT ~  LT+ LMC_t0  + SLA + Surface_F  + SD, data = BDD_ech)
summary(m3)
m4 <- lm(BT ~  LT+ LMC_t0  + SLA  + SD, data = BDD_ech)
summary(m4)
m5 <- lm(BT ~   LT+ LMC_t0 + SD, data = BDD_ech)
summary(m5)

m6 <- lm(BT ~   LT* LMC_t0 * SD, data = BDD_ech)
summary(m6)

m7 <- lm(BT ~   LT + LMC_t0 + SD + LT*SD + LT*LMC_t0, data = BDD_ech)
summary(m7)
m8 <- lm(BT ~   LT + LMC_t0 + SD + LT*SD, data = BDD_ech)
summary(m8)
m9 <- lm(BT ~   LT*SD, data = BDD_ech)
summary(m9)


AIC(m_nul,m1,m2,m3,m4,m5,m6,m7,m8,m9) 
plot(m6)
hist(resid(m6))



###### BB #########
#plot avec traits (VE)
plot(BDD_ech$LT, BDD_ech$BB, xlab = "LT", ylab = "BB")
plot(BDD_ech$LMC_t0, BDD_ech$BB, xlab = "LMC_t0", ylab = "BB")
plot(BDD_ech$LMC_t24, BDD_ech$BB, xlab = "LMC_t24", ylab = "BB")
plot(BDD_ech$Nb_ramifications, BDD_ech$BB, xlab = "Nombre de ramifications", ylab = "BB")
plot(BDD_ech$PEF, BDD_ech$BB, xlab = "PEF", ylab = "BB")
plot(BDD_ech$Gmin, BDD_ech$BB, xlab = "Gmin", ylab = "BB")
plot(BDD_ech$SLA, BDD_ech$BB, xlab = "SLA", ylab = "BB")
plot(BDD_ech$Surface_F, BDD_ech$BB, xlab = "Surface foliaire", ylab = "BB")
plot(BDD_ech$TDIA, BDD_ech$BB, xlab = "Diamètre de tige", ylab = "BB")
plot(BDD_ech$SD, BDD_ech$BB, xlab = "SD", ylab = "BB")

#modeles
m1 <- lm(BB ~  LT+ LMC_t0 + Nb_ramifications  + Gmin + SLA + Surface_F  + SD, data = BDD_ech)
summary(m1)
m2 <- lm(BB ~  LT+ LMC_t0 + Nb_ramifications + Gmin + SLA + Surface_F , data = BDD_ech)
summary(m2)
m3 <- lm(BB ~  LT+ LMC_t0 + Nb_ramifications + SLA + Gmin   , data = BDD_ech)
summary(m3)
m4 <- lm(BB ~  LT+ LMC_t0 + Nb_ramifications + SLA  , data = BDD_ech)
summary(m4)
m5 <- lm(BB ~   LT* LMC_t0 * Nb_ramifications * SLA  , data = BDD_ech)
summary(m5)
m6 <- lm(BB ~   LT* LMC_t0 * Nb_ramifications + LT* LMC_t0 * SLA + LT* Nb_ramifications * SLA + LMC_t0 * Nb_ramifications * SLA, data = BDD_ech)
summary(m6)
m7 <- lm(BB ~   LT* LMC_t0 + LT*Nb_ramifications + LT*SLA + LMC_t0*Nb_ramifications + LMC_t0*SLA + Nb_ramifications*SLA  , data = BDD_ech)
summary(m7)
m8 <- lm(BB ~   LT* LMC_t0 + Nb_ramifications + SLA, data = BDD_ech)
summary(m8)

m9 <- lm(BB ~   LT* LMC_t0 + Nb_ramifications , data = BDD_ech)
summary(m9)

m10 <- lm(BB ~   LT+ LMC_t0 + Nb_ramifications, data = BDD_ech)
summary(m10)

AIC(m2,m3,m4,m5,m6,m7,m8,m9,m10) 
plot(m9)
hist(resid(m9))












############### plot avec ecart-type #######################

plot(BDD_esp$LMC_t24,BDD_esp$MT,xlim=c(0,1000),ylim=c(100,900))
segments(x0=BDD_esp$LMC_t24-BDD_sd_esp$LMC_t24,
         x1=BDD_esp$LMC_t24+BDD_sd_esp$LMC_t24,
         y0=BDD_esp$MT,
         y1=BDD_esp$MT)
segments(x0=BDD_esp$LMC_t24,
         x1=BDD_esp$LMC_t24,
         y0=BDD_esp$MT-BDD_sd_esp$MT,
         y1=BDD_esp$MT+BDD_sd_esp$MT)
points(BDD_esp$LMC_t24,BDD_esp$MT,pch=21,bg="lightblue",cex=2)










################ graphique pour visualiser les espèces et leur inflammabilité ####################
library(ggplot2)

# définition du min et du max pour le graph 
scoremin <- min(BDD_esp$score_normalise)
scoremax <- max(BDD_esp$score_normalise)


# graphique
ggplot(BDD_esp, aes(x = reorder(Nom_scientifique, score_normalise), y = score_normalise, fill = score_normalise)) +
  geom_bar(stat = "identity", color = NA) +  
  coord_flip() +
  scale_fill_gradient2(
    low = "darkgreen",
    mid = "yellow",
    high = "red",
    midpoint = 0,
    limits = c(scoremin, scoremax),) +
  labs(x = "Espèces",
       y = "Score d'inflammabilité",
       fill = "Score",
       title = "Score d'inflammabilité par espèce"
      ) 



