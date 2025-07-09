
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




















##  Charger les packages nécessaires pour ACP
library(FactoMineR)
library(factoextra)


################  ACP INFLAMMABILITE ####################

# Sélection des colonnes des composantes de l'inflammabilité
colonnes_infla <- BDD_esp[, c(5,7,9,10)]

# Vérification des données
head(colonnes_infla)

# Centrage-réduction des données
colonnes_infla_cr <- scale(log(colonnes_infla))

# Vérification des données standardisées
head(colonnes_infla_cr)

# Application de l'ACP
res.pca <- PCA(colonnes_infla_cr, scale.unit = FALSE, graph = FALSE)

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

# Ajout d'un score d'inflammabilité basé sur la coordonnée de l'axe 1 
coord <- res.pca$ind$coord  # coordonnées des individus
BDD_esp$score <- coord[,1]


# normalisation du score d'inflammabilité entre -1 et 1
min_score <- min(BDD_esp$score, na.rm = TRUE)
max_score <- max(BDD_esp$score, na.rm = TRUE)
BDD_esp$score_normalise <- -1 + (BDD_esp$score - min_score) * 2 / (max_score - min_score)







################  ACP INFLAMMABILITE TEST ####################

# Sélection des colonnes des composantes de l'inflammabilité
colonnes_infla_test <- BDD_esp[, c(7,9,10,12)]

# Vérification des données
head(colonnes_infla_test)

# Centrage-réduction des données
colonnes_infla_test_cr <- scale(colonnes_infla_test)

# Vérification des données standardisées
head(colonnes_infla_test_cr)

# Application de l'ACP
res.pca <- PCA(colonnes_infla_test_cr, scale.unit = FALSE, graph = FALSE)

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
BDD_esp$score_test <- coord[,1]

# normalisation du score d'inflammabilité entre -1 et 1
min_score <- min(BDD_esp$score_test, na.rm = TRUE)
max_score <- max(BDD_esp$score_test, na.rm = TRUE)
BDD_esp$score_test_normalise <- -1 + (BDD_esp$score_test - min_score) * 2 / (max_score - min_score)

HCPC(res.pca)


distcoord<-dist(coord[,1:2])
clustcoord<-hclust(distcoord,method="ward.D2")
plot(clustcoord)
# visualiser les scores
View(BDD_esp)





#################### ACP TRAITS ############################

# Sélection des colonnes des traits fonctionnels
colonnes_traits <- na.omit(BDD_esp[, setdiff(13:27, c(20, 23))])

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

# Graphique des individus 
fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)

# Graphique combiné des variables et des individus
fviz_pca_biplot(res.pca, col.var = "contrib", gradient.cols =c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)

# Afficher l'ébouli
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50), main="Graphique de l'ébouli")






######## ACP INFLA avec projection des axes TRAITS ####################

# Sélection des lignes complètes (pas de NA dans les colonnes 3 à 22)
colonnes_complet <- na.omit(BDD_esp[, setdiff(5:29, c(6,8,11,19, 22,27,28,29,30))])
dim (colonnes_complet)
## ACP ##
# colonnes 1 à 4 : inflammabilité (axes actifs)
# colonnes 5 à 20 : traits (axes supplémentaires)
res.pca <- PCA(colonnes_complet, scale.unit = TRUE, 
               quanti.sup = 5:18, graph = FALSE)

# Visualisation ACP
fviz_pca_var(res.pca, col.var = "red", repel = TRUE) +
  ggtitle("Variables actives (infla) et supplémentaires (traits)")

# ACP avec les individus
fviz_pca_biplot(res.pca, col.var = "red",repel = TRUE)






######## ACP INFLA avec projection des axes TRAITS (TEST) ####################

# Sélection des lignes complètes (pas de NA dans les colonnes 3 à 22)
colonnes_complet_test <- na.omit(BDD_esp[, setdiff(5:26, c(5,6,8,11,20, 23))])
dim (colonnes_complet_test)

colonnes_complet_test_cr <- scale(colonnes_complet_test)

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


################# Heatmap TEST ##################

# packages nécessaires
library(ggplot2)
library(tidyr)  # pour changer le format de la BDD en format "long"

# Copier la BDD (pour ne pas la modifier)
df_prep <- BDD_esp

# Normalisser les valeur entre 0 et 1 : fonction
normalize <- function(x) {  return((x - min(x)) / (max(x) - min(x)))}

df_norm <- df_prep
df_norm$MT <- normalize(df_prep$MT)
df_norm$DI_test <- normalize(df_prep$DI_test)
df_norm$BB_test <- normalize(df_prep$BB_test)
df_norm$BT_test <- normalize(df_prep$BT_test)

# Transformer en format long (une ligne par composante)
df_long <- pivot_longer(df_norm,
                        cols = c("MT", "DI_test", "BB_test", "BT_test"),
                        names_to = "Variable",
                        values_to = "Valeur")

# Créer l’ordre des espèces en fonction du score d'inflammabilité
ordre_esp <- df_prep[order(df_prep$score_test_normalise, decreasing = FALSE), "Nom_scientifique"]

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
mat_cor_trait<-cor(colonnes_traits,method="spearman")
corrplot(mat_cor_trait)
round(mat_cor_trait, 2)   ## afficher les valeurs
# Visualiser la corrélation
corrplot(mat_cor_trait, method = "color", type = "upper", tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")
corrplot(mat_cor_trait, method = "color", tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")

plot(log(BDD_esp$Surface_F)~BDD_esp$Nb_rami)

hist(BDD_esp$Surface_F)











#################### TEST ANOVA POUR VARIABILITE INTRA et INTER ESPECES ######################
# Variables d'intérêt
variables <- c("BB_test", "BT_test", "MT", "DI_test", "TD", "Surface_F", "LMC_t24", "SD", "LT", "SLA")

# Initialisation du tableau résultat simplifié
results_df <- data.frame(
  variable = character(),
  prop_inter = numeric(),
  prop_intra = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Boucle sur chaque variable
for (i in variables) {
  
# Log-transformation + 0.01
y <- log(BDD_ech[[i]] + 0.01)
  
# Formule pour ANOVA
formula <- as.formula(paste("y ~ Nom_scientifique"))
  
# modèle anova
fit <- aov(formula, data = BDD_ech)
  
# Résumé ANOVA
anova_res <- summary(fit)[[1]]
  
# Somme des carrés pour l'espèce et le résiduel
ss_species <- anova_res["Nom_scientifique", "Sum Sq"]
ss_residual <- anova_res["Residuals", "Sum Sq"]
ss_total <- ss_species + ss_residual
  
# Proportions de variance
prop_inter <- ss_species / ss_total
prop_intra <- ss_residual / ss_total
  
# p-value
p_value <- anova_res["Nom_scientifique", "Pr(>F)"]
  
# Ajouter au tableau résultat simplifié
results_df <- rbind(results_df, data.frame(
  variable = v,
  prop_inter = prop_inter,
  prop_intra = prop_intra,
  p_value = p_value,
  stringsAsFactors = FALSE
  ))
}

# Afficher les résultats simplifiés
print(results_df)

summary(fit)






















################### matrice de corrélation pour observer l'effet des traits sur infla ################
### de tout 
library (corrplot)
mat_cor_all<-cor(colonnes_complet,method="spearman")
corrplot(mat_cor_all)
round(mat_cor_all, 2)   ## afficher les valeurs
corrplot(mat_cor_all, method = "color", tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")












#################### PLOT DISTRIBUTION ##################################
# histogramme pour vérifier la normalité des variables réponse (infla)
hist(BDD_ech$MT)       #normale
hist(sqrt(BDD_ech$MT))
shapiro.test(BDD_ech$MT)
hist(BDD_ech$BT_test)       #necéssite peut être tranfo log 
hist(log(BDD_ech$BT_test + 1))  #normale

hist(BDD_ech$DI_test)#voir je trouve pas la bonne fonction

hist(BDD_ech$BB) #voir
hist(log(BDD_ech$BB+1))
shapiro.test(log(BDD_ech$BB+1))




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
