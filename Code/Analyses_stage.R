
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
colonnes_infla <- BDD_esp[, c(4:7)]

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
colonnes_traits <- na.omit(BDD_esp[,c(8:22)])

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
colonnes_complet <- na.omit(BDD_esp[, 4:22])

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
corrplot(mat_cor_trait, method = "color", type = "upper",
         tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")
corrplot(mat_cor_trait, method = "color", 
         tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")
# sélection automatique des traits
library(caret)
traits_non_corrélés <- colonnes_traits[, -findCorrelation(mat_cor_trait, cutoff = 0.6)]
traits_non_corrélés







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
  #PET
variabilité_esp <- aov(BDD_ech$PET~BDD_ech$Nom_scientifique)
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
library(lme4)

# Modèle linéaire mixte
m_BT <- lmer(LMC_t24 ~ 1 + (1 | Nom_scientifique), data = BDD_ech)
summary(m_BT)

# Extraire les composantes de variance
var_BT <- as.data.frame(VarCorr(m_BT))
var_BT_inter <- var_BT$vcov[1]
var_BT_intra <- attr(VarCorr(m_BT), "sc")^2

# Pourcentage de variance
total_BT <- var_BT_inter + var_BT_intra
100 * var_BT_inter / total_BT  # variance inter (espèce)
100 * var_BT_intra / total_BT  # variance intra (résiduelle)




################### matrice de corrélation pour observer l'effet des traits sur infla ################



### de tout 
library (corrplot)
mat_cor_all<-cor(colonnes_complet)
corrplot(mat_cor_all)
round(mat_cor_all, 2)   ## afficher les valeurs
corrplot(mat_cor_all, method = "color", 
         tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")








#################### PLOT ##################################

#histogrammes pour vérifier la normalité des données 
hist(BDD_esp$score2)

hist(BDD_ech$MT)

hist(BDD_esp$BT)
logBT<-log(BDD_esp$BT)
hist(logBT)


hist(BDD_esp$BB)

hist(BDD_esp$DI)
BDD_esp$sqrtDI=sqrt(BDD_esp$DI)
hist(BDD_esp$sqrtDI)




### LMC-t24 avec les 4 composantes infla ####
plot(BDD_esp$LMC_t24,BDD_esp$BT)
text(BDD_esp$LMC_t24,BDD_esp$BT,BDD_esp$Nom_scientifique)

plot(BDD_esp$LMC_t24,BDD_esp$BB)
text(BDD_esp$LMC_t24,BDD_esp$BB,BDD_esp$Nom_scientifique)

plot(BDD_ech$LMC_t24,BDD_ech$MT)
text(BDD_esp$LMC_t24,BDD_esp$MT,BDD_esp$Nom_scientifique)

plot(BDD_esp$LMC_t24,BDD_esp$DI)
text(BDD_esp$LMC_t24,BDD_esp$DI,BDD_esp$Nom_scientifique)

#### SLA ####
plot(BDD_esp$SLA,BDD_esp$BT)
text(BDD_esp$SLA,BDD_esp$BT,BDD_esp$Nom_scientifique)

#### Nb rami #####
plot(BDD_esp$Nb_ramifications,BDD_esp$BB)
text(BDD_esp$Nb_ramifications,BDD_esp$BB,BDD_esp$Nom_scientifique)

#### surface #####
plot(BDD_esp$Surface_F,BDD_esp$BT)
text(BDD_esp$Surface_F,BDD_esp$BT,BDD_esp$Nom_scientifique)
plot(BDD_esp$Surface_F,BDD_esp$score2)
text(BDD_esp$Surface_F,BDD_esp$score2,BDD_esp$Nom_scientifique,)

#### score ###
plot(BDD_esp$score2,BDD_esp$score2)
text(BDD_esp$score2,BDD_esp$score2,BDD_esp$Nom_scientifique,)

plot(BDD_esp$score2, BDD_esp$score2, type = "n")  # Ne trace que l'espace vide
text(BDD_esp$score2, BDD_esp$score2, BDD_esp$Nom_scientifique, cex = 0.5, srt = -45, pos = 4)

library(ggplot2)
library(ggrepel)

ggplot(BDD_esp, aes(x = score2, y = score2)) + 
  geom_point() +
  geom_text_repel(aes(Nom_scientifique), size = 3) +
  theme_minimal()

m<-lm(score2~LMC_t24,data=BDD_esp)
summary(m)
plot(m)

anova(m)

hist(BDD_ech$MT)
hist(BDD_ech$BT)
hist(BDD_ech$DI)
hist(BDD_ech$BB)






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


################# Heatmap (couleur en fonction de chaque composante d'inflammabilité) ##################

# packages nécessaires
library(ggplot2)
library(tidyr)  # pour changer le format de la BDD en format "long"

# Copier la BDD
df_prep <- BDD_esp

# Normalisser les valeur entre 0 et 1 : fonction
normalize <- function(x) {  return((x - min(x)) / (max(x) - min(x)))}

df_norm <- df_prep
df_norm$MT <- normalize(df_prep$MT)
df_norm$DI <- normalize(df_prep$DI)
df_norm$BB <- normalize(df_prep$BB)
df_norm$BT <- normalize(df_prep$BT)

# Transformer en format long
df_long <- pivot_longer(df_norm,
                        cols = c("MT", "DI", "BB", "BT"),
                        names_to = "Variable",
                        values_to = "Valeur")
#définir lordre des espèces en fonction du score d'inflammabilité calculé plus haut
ordre_esp <- BDD_esp[order(BDD_esp$score_normalise), "Nom_scientifique"]


# Facteur pour trier les espèces de haut en bas (du plus inflammable au moins inflammable)
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
       y = "Espèces") 

