
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
dim(BDD_esp)
names(BDD_esp)[which(names(BDD_esp) == "Nb_ramifications")] <- "Nb_rami"
dim(BDD_esp)
View(BDD_esp)

################## DISTRIBUTION DES DONNEES #####################
#Infla
hist(BDD_esp$DI)      # plus de grande valeurs que petites
hist(BDD_esp$DI_test) # same DI
hist(scale(log(BDD_esp$BT)))      # plus de petites valeurs que grande
hist(BDD_esp$BT_test) # same BT
hist(scale(BDD_esp$BB))      # normale
hist(qlogis(BDD_esp$BB_test)) # same BB
hist(scale(BDD_esp$MT))      # normale 

#traits

hist(log(BDD_esp$Nb_rami))  #normale
hist(log(BDD_esp$SD))       #normale
hist(log(BDD_esp$TMC_t0))         #normale
hist(log(BDD_esp$TMC_t24))   #plutôt normale
hist(log(BDD_esp$TDMC))           #normale
hist(log(BDD_esp$TD))             #normale
hist(log(BDD_esp$TDIA))           #normale
hist(BDD_esp$LMC_t0)         #plutôt normale
hist(log(BDD_esp$LMC_t24))   #normale
hist(BDD_esp$LDMC)           #normale
hist(log(BDD_esp$Surface_F)) #normale
hist(BDD_esp$SLA)            #normale
hist(log(BDD_esp$LT))        #normale





############################################################## ECHELLE ESPECE #############################################################

################## DISTRIBUTION DES DONNEES #####################
#Infla
hist(BDD_esp$DI)      # asymétrique gauche
hist(BDD_esp$DI_test) # asymétrique droite
hist(BDD_esp$BT)      # asymétrique droite
hist(BDD_esp$BT_test) # asymétrique droite
hist(BDD_esp$BB)      # normale proportion
hist(BDD_esp$BB_test) # normale proportion
hist(BDD_esp$MT)      # normale 

hist(BDD_esp$T_ambiante)
plot(BDD_esp$T_ambiante,BDD_esp$DI_test)

#traits

hist(BDD_esp$Nb_rami)   #asymétrique droite
hist(BDD_esp$SD)        #asymétrique droite
hist(BDD_esp$TMC_t0)    #asymétrique droite ou normale (hésitation)
hist(BDD_esp$TMC_t24)   #asymétrique droite ou normale (hésitation)
hist(BDD_esp$TDMC)      #normale
hist(BDD_esp$TD)        #normale
hist(BDD_esp$TDIA)      #normale
hist(BDD_esp$LMC_t0)    #asymétrique droite ou normale (hésitation) 
hist(BDD_esp$LMC_t24)   #asymétrique droite ou normale (hésitation) 
hist(BDD_esp$LDMC)      #normale
hist(BDD_esp$Surface_F) #asymétrique droite
hist(BDD_esp$SLA)       #normale
hist(BDD_esp$LT)        #asymétrique droite









################ EFFET DES CONDITIONS METEO ##################
#MT
options(na.action = "na.omit")
m<-glm(MT~T_ambiante+Vent+Humidite,data=BDD_esp,family="gaussian")
summary(m)

pred<-predict(m,type = "response")
plot(BDD_esp$MT,pred)
abline(a=0,b=1)


#BB 
BDD_esp$BB_prop <- BDD_esp$BB_test/100
m<-glm(BB_prop~T_ambiante+Vent+Humidite,data=BDD_esp,family="binomial")
summary(m)

pred<-predict(m,type = "response")
plot(BDD_esp$BB_prop,pred)
abline(a=0,b=1)


#DI
m<-glm(DI_test~T_ambiante+Vent,data=BDD_esp,family="Gamma")
summary(m)

pred<-predict(m)
plot(BDD_esp$DI_test,pred)
abline(a=0,b=1)

plot(BDD_esp$T_ambiante,1/BDD_esp$DI_test)


#BT
m<-glm(BT_test~T_ambiante,data=BDD_esp,family="Gamma")
summary(m)

pred<-predict(m,type = "response")
plot(BDD_esp$BT,pred)
abline(a=0,b=1)

plot(BDD_esp$Vent,BDD_esp$BT_test)










############################### ACP ######################################

##  Charger les packages nécessaires pour ACP
library(FactoMineR)
library(factoextra)


################  ACP INFLAMMABILITE ####################

# Sélection des colonnes des composantes de l'inflammabilité
colonnes_infla <- BDD_esp[, c(10,12,13,15)]

# Vérification des données
head(colonnes_infla)

# Centrage-réduction des données
colonnes_infla_cr <- scale(colonnes_infla)

# Vérification des données standardisées
head(colonnes_infla_cr)

# Application de l'ACP
res.pca <- PCA(colonnes_infla_cr, scale.unit = FALSE, graph = FALSE)

# Résumé des résultats
summary(res.pca)


# Graphique des contributions des variables aux composantes principales
fviz_pca_var(res.pca, col.var = "contrib", gradient.cols = c("#6fec00", "#ff9e00", "#Ff0000"),repel = TRUE)


# Graphique combiné des variables et des individus
fviz_pca_biplot(res.pca,col.var = "contrib", gradient.cols =c("#00AFBB", "#E7B800", "#FC4E07") , repel = TRUE)

# Afficher l'ébouli
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50), main="Graphique de l'ébouli")

# Ajout d'un score d'inflammabilité basé sur la coordonnée de l'axe 1 
coord <- res.pca$ind$coord  # coordonnées des individus
BDD_esp$score <- coord[,1]


# normalisation du score d'inflammabilité entre -1 et 1
min_score <- min(BDD_esp$score, na.rm = TRUE)
max_score <- max(BDD_esp$score, na.rm = TRUE)
BDD_esp$score_normalise <- -1 + (BDD_esp$score - min_score) * 2 / (max_score - min_score)


HCPC(res.pca,method="ward")


distcoord<-dist(coord[,1:2])
clustcoord<-hclust(distcoord,method="ward.D2")
plot(clustcoord)
# visualiser les scores
View(BDD_esp)






####################### CLASSEMENT ESPECES ############################

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






#################### ACP TRAITS ############################

# Sélection des colonnes des traits fonctionnels
colonnes_traits <- na.omit(BDD_esp[, setdiff(16:30, c(23, 26))])

# Vérifier les données
colonnes_traits


colonnes_traits_cr <- scale(log(colonnes_traits))
colonnes_traits_cr
# Application de l'ACP
res.pca <- PCA(colonnes_traits_cr, scale.unit = FALSE, graph = FALSE)

# Résumé des résultats
summary(res.pca)

# Graphique des contributions des variables aux composantes principales
fviz_pca_var(res.pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

# Graphique combiné des variables et des individus
fviz_pca_biplot(res.pca, col.var = "contrib", gradient.cols =c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)

# Afficher l'ébouli
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50), main="Graphique de l'ébouli")


####################### matrice de corrélation (pour le choix des traits) ######################

### des traits foncitonnels
library (corrplot)
mat_cor_trait<-cor(colonnes_traits,method="spearman")
corrplot(mat_cor_trait)
round(mat_cor_trait, 2)   ## afficher les valeurs
# Visualiser la corrélation
corrplot(mat_cor_trait, method = "color", type = "upper", tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")
corrplot(mat_cor_trait, method = "color", tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")

library(caret)
traits_non_corrélés <- colonnes_traits[, -findCorrelation(mat_cor_trait, cutoff = 0.7)]
traits_non_corrélés










######## ACP INFLA avec projection des axes TRAITS (TEST) ####################

# Sélection des lignes complètes (pas de NA dans les colonnes 3 à 22)
colonnes_complet_test <- na.omit(BDD_esp[, setdiff(8:30, c(8,9,11,14,23, 26))])
dim (colonnes_complet_test)

colonnes_complet_test_cr <- scale(colonnes_complet_test)
head(colonnes_complet_test_cr)
## ACP ##
# colonnes 1 à 4 : inflammabilité (axes actifs)
# colonnes 5 à 20 : traits (axes supplémentaires)
res.pca <- PCA(colonnes_complet_test_cr, scale.unit = FALSE, 
               quanti.sup = 5:17, graph = FALSE)

# Visualisation ACP
fviz_pca_var(res.pca, col.var = "red", repel = TRUE) +
  ggtitle("Variables actives (infla) et supplémentaires (traits)")

# ACP avec les individus
fviz_pca_biplot(res.pca, col.var = "red",repel = TRUE)









############# PLOT INFLA ESPECE ######################

library(ggplot2)



# Calcul des moyennes par espèce
moyennes_MT <- aggregate(MT ~ Nom_scientifique, data = BDD_ech, FUN = mean, na.rm = TRUE)

# Ordonner les espèces croissante
moyennes_MT <- moyennes_MT[order(moyennes_MT$MT), ]
BDD_ech$Nom_scientifique <- factor(BDD_ech$Nom_scientifique, levels = moyennes_MT$Nom_scientifique)

# Créer un dataframe pour les moyennes (pour ggplot)
moyennes_MT$Type <- "Moyenne"
BDD_ech$Type <- "Individuel"

# Plot avec légende
MT<-ggplot() +
  geom_point(data = BDD_ech, aes(x = MT, y = Nom_scientifique, color = Type), size = 2) +
  geom_point(data = moyennes_MT, aes(x = MT, y = Nom_scientifique, color = Type), size = 3) +
  scale_color_manual(values = c("Individuel" = "black", "Moyenne" = "#f14900")) +
  labs(x = "MT", y = "Espèces", color = "") +
  theme_bw()
MT


# Calcul des moyennes par espèce
moyennes_DI_test <- aggregate(DI_test ~ Nom_scientifique, data = BDD_ech, FUN = mean, na.rm = TRUE)

# Ordonner les espèces croissante (en fonction de MT)
moyennes_DI_test <- moyennes_DI_test[order(moyennes_MT$MT),  ]
BDD_ech$Nom_scientifique <- factor(BDD_ech$Nom_scientifique, levels = moyennes_DI_test$Nom_scientifique)

# Créer un dataframe pour les moyennes (pour ggplot)
moyennes_DI_test$Type <- "Moyenne"
BDD_ech$Type <- "Individuel"

# Plot avec légende
DI<- ggplot() +
  geom_point(data = BDD_ech, aes(x = DI_test, y = Nom_scientifique, color = Type), size = 2) +
  geom_point(data = moyennes_DI_test, aes(x = DI_test, y = Nom_scientifique, color = Type), size = 3) +
  scale_color_manual(values = c("Individuel" = "black", "Moyenne" = "#f14900")) +
  labs(x = "DI", y = "Espèce", color = "") +
  theme_bw()
DI


# Calcul des moyennes par espèce
moyennes_BT_test <- aggregate(BT_test ~ Nom_scientifique, data = BDD_ech, FUN = mean, na.rm = TRUE)

# Ordonner les espèces croissante
moyennes_BT_test <- moyennes_BT_test[order(moyennes_MT$MT), ]
BDD_ech$Nom_scientifique <- factor(BDD_ech$Nom_scientifique, levels = moyennes_BT_test$Nom_scientifique)

# Créer un dataframe pour les moyennes (pour ggplot)
moyennes_BT_test$Type <- "Moyenne"
BDD_ech$Type <- "Individuel"

# Plot avec légende
BT<-ggplot() +
  geom_point(data = BDD_ech, aes(x = BT_test, y = Nom_scientifique, color = Type), size = 2) +
  geom_point(data = moyennes_BT_test, aes(x = BT_test, y = Nom_scientifique, color = Type), size = 3) +
  scale_color_manual(values = c("Individuel" = "black", "Moyenne" = "#f14900")) +
  labs(x = "BT", y = "Espèce", color = "") +
  theme_bw()
BT


# Calcul des moyennes par espèce
moyennes_BB_test <- aggregate(BB_test ~ Nom_scientifique, data = BDD_ech, FUN = mean, na.rm = TRUE)

# Ordonner les espèces croissante
moyennes_BB_test <- moyennes_BB_test[order(moyennes_MT$MT),  ]
BDD_ech$Nom_scientifique <- factor(BDD_ech$Nom_scientifique, levels = moyennes_BB_test$Nom_scientifique)

# Créer un dataframe pour les moyennes (pour ggplot)
moyennes_BB_test$Type <- "Moyenne"
BDD_ech$Type <- "Individuel"

# Plot avec légende
BB<-ggplot() +
  geom_point(data = BDD_ech, aes(x = BB_test, y = Nom_scientifique, color = Type), size = 2) +  
  geom_point(data = moyennes_BB_test, aes(x = BB_test, y = Nom_scientifique, color = Type), size = 3)+
  scale_color_manual(values = c("Individuel" = "black", "Moyenne" = "#f14900")) +
  labs(x = "BB", y = "Espèce", color = "") +
  theme_bw()
BB















################ BOXPLOT INFLA ESPECE ########################

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


######## DI_TEST ###########
# Calcul des médianes par espèce
med_DI_test <- aggregate(DI_test ~ Nom_scientifique, data = BDD_ech, median)

# Création ordre décroissant des médianes
ordre <- med_DI_test$Nom_scientifique[order(med_DI_test$DI_test, decreasing = TRUE)]

# facteur avec ce nouvel ordre (pour que l'ordre soit capté par le boxplot)
BDD_ech$Nom_scientifique <- factor(BDD_ech$Nom_scientifique, levels = ordre)

# Traçage du boxplot avec l'ordre décroissant 
boxplot(DI_test ~ Nom_scientifique, data = BDD_ech,
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


######## BT_TEST ###########
# Calcul des médianes par espèce
med_BT_test <- aggregate(BT_test ~ Nom_scientifique, data = BDD_ech, median)

# Création ordre décroissant des médianes
ordre <- med_BT_test$Nom_scientifique[order(med_BT_test$BT_test, decreasing = TRUE)]

# facteur avec ce nouvel ordre (pour que l'ordre soit capté par le boxplot)
BDD_ech$Nom_scientifique <- factor(BDD_ech$Nom_scientifique, levels = ordre)

# Traçage du boxplot avec l'ordre décroissant 
boxplot(BT_test ~ Nom_scientifique, data = BDD_ech,
        las = 2, cex.axis = 0.7,
        main = "BT par espèce (ordre décroissant de la médiane)",
        xlab = "Espèce", ylab = "BT")



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


######## BB_TEST ###########
# Calcul des médianes par espèce
med_BB_test <- aggregate(BB_test ~ Nom_scientifique, data = BDD_ech, median)

# Création ordre décroissant des médianes
ordre <- med_BB_test$Nom_scientifique[order(med_BB_test$BB_test, decreasing = TRUE)]

# facteur avec ce nouvel ordre (pour que l'ordre soit capté par le boxplot)
BDD_ech$Nom_scientifique <- factor(BDD_ech$Nom_scientifique, levels = ordre)

# Traçage du boxplot avec l'ordre décroissant 
boxplot(BB_test ~ Nom_scientifique, data = BDD_ech,
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
df_norm$score_DI <- normalize(df_prep$score_DI)
df_norm$BB_test <- normalize(df_prep$BB_test)
df_norm$BT_test <- normalize(df_prep$BT_test)
View(df_norm)
# Transformer en format long (une ligne par composante)
df_long <- pivot_longer(df_norm,
                        cols = c("MT", "score_DI", "BB_test", "BT_test"),
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























#################### TEST ANOVA POUR VARIABILITE INTRA et INTER ESPECES ######################
#  Liste des variables à analyser 
vars <- c("BB_test", "BT_test", "MT", "DI_test", "TD", "Surface_F", "Nb_ramifications",
          "LMC_t24", "LDMC", "SD", "LT", "SLA")

#  Transformation log(var + 0.01)
BDD_ech_log <- BDD_ech
BDD_ech_log[vars] <- log(BDD_ech[vars] + 0.01)

# 3. Initialisation d'une table pour stocker les résultats
var_partition <- data.frame(Variable = character(),
                            p_value = numeric(),
                            Prop_intra = numeric(),
                            Prop_inter = numeric(),
                            stringsAsFactors = FALSE)

# 4. Boucle sur chaque variable pour faire l'ANOVA et extraire les variances
for (var in vars) {
  # Formule de l'ANOVA
  collage <- as.formula(paste(var, "~ Nom_scientifique"))
  
  # ANOVA
  anova_var <- aov(collage, data = BDD_ech_log)
  aov_summary <- summary(anova_var)[[1]]
  
  # Extraction des sommes de carrés
  SS_among <- aov_summary["Nom_scientifique", "Sum Sq"]
  SS_within <- aov_summary["Residuals", "Sum Sq"]
  SS_total <- SS_among + SS_within
  
  # Stockage des résultats
  var_partition <- rbind(var_partition, data.frame(
    Variable = var,
    p_value = aov_summary["Nom_scientifique", "Pr(>F)"],
    Prop_inter = SS_among / SS_total,
    Prop_intra = SS_within / SS_total
  ))
}

# 5. Affichage des résultats
print(var_partition)




#################### GLM POUR VARIABILITE INTRA et INTER ESPECES ######################
#chargement pckage 
library(lme4)

# BT_test : durée de combustion (Gamma)
mod_BT <- glmer(BT_test ~ 1 + (1 | Nom_scientifique),family = "Gamma", data = BDD_ech)
summary(mod_BT)


# DI_test : délai d’inflammation (Gamma)
mod_DI <- glmer(DI_test ~ 1 + (1 | Nom_scientifique),family = "Gamma", data = BDD_ech)
summary(mod_DI)


# MT : température max (Gaussian)
mod_MT <- lmer(MT ~ 1 + (1 | Nom_scientifique), data = BDD_ech)
summary(mod_MT)


# BB_prop : proportion brûlée (Binomiale)
BDD_ech$BB_prop <- BDD_ech$BB_test/100
mod_BB <- glmer(BB_prop ~ 1 + (1 | Nom_scientifique),family = "binomial", data = BDD_ech)
summary(mod_BB)






















################### matrice de corrélation pour observer l'effet des traits sur infla ################
### de tout 
library (corrplot)
mat_cor_all<-cor(colonnes_complet,method="spearman")
corrplot(mat_cor_all)
round(mat_cor_all, 2)   ## afficher les valeurs
corrplot(mat_cor_all, method = "color", tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")











############################ MODELES #####################################
#standardisation des données 
SD_cr<-scale(BDD_esp$SD)
TD_cr<-scale(BDD_esp$TD)
LA_cr<-scale(BDD_esp$Surface_F)
LDMC_cr<-scale(BDD_esp$LDMC)
LT_cr<-scale(BDD_esp$LT)



###### MT #########
#plot avec traits (VE)
plot(BDD_esp$TD, BDD_esp$MT)
plot(BDD_esp$LT, BDD_esp$MT)
plot(BDD_esp$LDMC, BDD_esp$MT)
plot(BDD_esp$Surface_F, BDD_esp$MT)
plot(BDD_esp$SD, BDD_esp$MT)

#modèles
m_MT0<-glm(MT~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")
summary(m_MT0)

m_MT1<-glm(MT~SD_cr+TD_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")        # meilleur modèle
summary(m_MT1)

m_MT2<-glm(MT~SD_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")
summary(m_MT2)

m_MT3<-glm(MT~SD_cr+LDMC_cr,data=BDD_esp,family="gaussian")
summary(m_MT3)

AIC(m_MT0,m_MT1,m_MT2,m_MT3) 


pred<-predict(m,type = "response")
plot(BDD_esp$MT,pred)
abline(a=0,b=1)



###### BB #########
#plot avec traits (VE)
plot(BDD_esp$TD, BDD_esp$BB_test)
plot(BDD_esp$LT, BDD_esp$BB_test)
plot(BDD_esp$LDMC, BDD_esp$BB_test)
plot(BDD_esp$Surface_F, BDD_esp$BB_test)
plot(BDD_esp$SD, BDD_esp$BB_test)

#modèles
BDD_esp$BB_prop <- BDD_esp$BB_test/100

m_BB0<-glm(BB_prop~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp,family="binomial")
summary(m_BB0)

m_BB1<-glm(BB_prop~LDMC_cr,data=BDD_esp,family="binomial")        # meilleur modèle
summary(m_BB1)

AIC(m_BB0,m_BB1) #!!! pas le même nombre d'observation car SD a deux esp de moins 

pred<-predict(m,type = "response")
plot(BDD_esp$BB_prop,pred)



###### BT #########
#plot avec traits (VE)
plot(BDD_esp$TD, BDD_esp$BT_test)
plot(BDD_esp$LT, BDD_esp$BT_test)
plot(BDD_esp$LDMC, BDD_esp$BT_test)
plot(BDD_esp$Surface_F, BDD_esp$BT_test)
plot(BDD_esp$SD, BDD_esp$BT_test)

#modèles
m_BT0<-glm(BT_test~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp,family="Gamma")
summary(m_BT0)

m_BT1<-glm(BT_test~SD_cr+TD_cr+LDMC_cr,data=BDD_esp,family="Gamma")       
summary(m_BT1)

m_BT2<-glm(BT_test~TD_cr+LDMC_cr,data=BDD_esp,family="Gamma")
summary(m_BT2)

AIC(m_BT0,m_BT1) #!!! pas le même nombre d'observation car SD a deux esp de moins 


pred<-predict(m,type = "response")
plot(BDD_esp$MT,pred)
abline(a=0,b=1)



###### DI #########
#plot avec traits (VE)
plot(BDD_esp$TD, BDD_esp$DI_test)
plot(BDD_esp$LT, BDD_esp$DI_test)
plot(BDD_esp$LDMC, BDD_esp$DI_test)
plot(BDD_esp$Surface_F, BDD_esp$DI_test)
plot(BDD_esp$SD, BDD_esp$DI_test)

#modèles
m_DI0<-glm(DI_test~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp,family="Gamma")
summary(m_DI0)

m_DI1<-glm(DI_test~TD_cr+LDMC_cr+LT_cr,data=BDD_esp,family="Gamma")       
summary(m_DI1)

AIC(m_DI0,m_DI1) #!!! pas le même nombre d'observation car SD a deux esp de moins 



###### Score ########
hist(BDD_esp$score_normalise)
#plot avec traits (VE)
plot(BDD_esp$TD, BDD_esp$score_normalise)
plot(BDD_esp$LT, BDD_esp$score_normalise)
plot(BDD_esp$LDMC, BDD_esp$score_normalise)
plot(BDD_esp$Surface_F, BDD_esp$score_normalise)
plot(BDD_esp$SD, BDD_esp$score_normalise)

#modèles
m_score0<-glm(score_normalise~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")
summary(m_score0)

m_score1<-glm(score_normalise~SD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")       
summary(m_score1)

m_score2<-glm(score_normalise~SD_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")       
summary(m_score2)

m_score3<-glm(score_normalise~LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")  ##meilleur modèle     
summary(m_score3)

m_score4<-glm(score_normalise~LDMC_cr,data=BDD_esp,family="gaussian")       
summary(m_score4)

AIC(m_score4,m_score3)



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














BDD_esp$BB_prop <- BDD_esp$BB_test/100

m<-glm(BB_prop~Surface_F+TD+LMC_t24+LT,data=BDD_esp,family="binomial")
summary(m)
pred<-predict(m,type = "response")
plot(BDD_esp$BB_prop,pred)
length(BDD_esp$BB_prop)
length(pred)
abline(a=0,b=1)


options(na.action = "na.omit")
m<-glm(MT~LMC_t24+LDMC+LT+Nb_rami+SLA+Surface_F+TD,data=BDD_esp,family="gaussian")
summary(m)
pred<-predict(m,type = "response")
plot(BDD_esp$MT,pred)
abline(a=0,b=1)

library(MuMIn)

dredge(m)
test<-na.exclude(BDD_esp)
plot(BDD_esp$LDMC,BDD_esp$LMC_t24)

m<-glm(LMC_t24~Gmin+LDMC,data=BDD_esp,family="Gamma")
summary(m)
plot(BDD_esp$Gmin,BDD_esp$LMC_t24)


foo<-seq(min(BDD_esp$Gmin),max(BDD_esp$Gmin),1)
pred_data1<-data.frame(Gmin=foo,
                       LDMC=rep(quantile(BDD_esp$LDMC,p=0.25),length(foo)))

pred_data2<-data.frame(Gmin=foo,
                       LDMC=rep(quantile(BDD_esp$LDMC,p=0.75),length(foo)))

pred_data3<-data.frame(Gmin=foo,
                       LDMC=rep(quantile(BDD_esp$LDMC,p=0.05),length(foo)))


pred1<-predict(m,type="response",newdata=pred_data1,se.fit=TRUE)
pred2<-predict(m,type="response",newdata=pred_data2,se.fit=TRUE)
pred3<-predict(m,type="response",newdata=pred_data3,se.fit=TRUE)


pred

plot(BDD_esp$LMC_t24~BDD_esp$Gmin)
lines(pred1$fit~pred_data$Gmin,lwd=2)
lines(pred2$fit~pred_data$Gmin,lwd=2)
lines(pred3$fit~pred_data$Gmin,lwd=2)



lines(pred$fit~pred_data$Gmin,lwd=2)
lines(pred$fit+(1.96*pred$se.fit)~pred_data$Gmin,lty=3)
lines(pred$fit-(1.96*pred$se.fit)~pred_data$Gmin,lty=3)


hist(BDD_esp$LMC_t24)
