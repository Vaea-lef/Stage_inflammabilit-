
#### ANALYSE DE DONNEE STAGE INFLAMMABILITE ####
#packages
library(FactoMineR)
library(factoextra)
library (corrplot)




################# IMPORTATION ET CHARGEMENT PACKAGE ###########################

##  importation BDD_ana_ech en format CSV
setwd("C:/IRD/Stage_inflammabilit-") #définition du répertoire de travail
BDD_ech<-read.csv2("Data/BDD_moy_ech.csv", header = TRUE) #importation de l
BDD_ech
names(BDD_ech)[which(names(BDD_ech) == "Nb_ramifications")] <- "Nb_rami"
View(BDD_ech)

BDD_echMT<-read.csv2("Data/BDD_moy_ech1.csv", header = TRUE) #importation de l
BDD_echMT
names(BDD_echMT)[which(names(BDD_echMT) == "Nb_ramifications")] <- "Nb_rami"
View(BDD_echMT)

BDD_ech3<-read.csv2("Data/BDD_moy_ech2.csv", header = TRUE) #importation de l
BDD_ech3
names(BDD_ech3)[which(names(BDD_ech3) == "Nb_ramifications")] <- "Nb_rami"
View(BDD_ech3)

##  importation BDD_sd_esp en format CSV
BDD_sd_esp<-read.csv2("Data/BDD_sd_esp.csv", header = TRUE)
BDD_sd_esp

##  importation BDD_moy_esp en format CSV
BDD_esp<-read.csv2("Data/BDD_moy_esp.csv", header = TRUE) #importation de la base
BDD_esp
dim(BDD_esp)
names(BDD_esp)[which(names(BDD_esp) == "Nb_ramifications")] <- "Nb_rami"
View(BDD_esp)

BDD_esp_netMT<-read.csv2("Data/BDD_moy_esp1.csv", header = TRUE) #importation de la base
BDD_esp_netMT
dim(BDD_esp_netMT)
names(BDD_esp_netMT)[which(names(BDD_esp_netMT) == "Nb_ramifications")] <- "Nb_rami"
View(BDD_esp_netMT)

BDD_esp_net3<-read.csv2("Data/BDD_moy_esp2.csv", header = TRUE) #importation de la base
BDD_esp_net3
dim(BDD_esp_net3)
names(BDD_esp_net3)[which(names(BDD_esp_net3) == "Nb_ramifications")] <- "Nb_rami"
View(BDD_esp_net3)

################## DISTRIBUTION DES DONNEES #####################
#Infla
hist(BDD_esp$DI)      # plus de grande valeurs que petites
hist(BDD_esp$DI_test) # same DI
hist(scale(log(BDD_esp$BT)))      # plus de petites valeurs que grande
hist(BDD_esp$BT_test) # same BT
hist(scale(BDD_esp$BB))      # normale
hist(qlogis(BDD_esp$BB_test)) # same BB
hist(scale(BDD_esp$MT))      # normale 

#Infla
hist(BDD_ech3$DI_test)      # plus de grande valeurs que petites
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
hist((BDD_esp$LT))        #normale





######################################################## ECHELLE ESPECE #############################################################

################## DISTRIBUTION DES DONNEES #####################
#Infla
hist(BDD_esp$DI)      # asymétrique gauche
hist(BDD_esp$DI_test,xlab="DI",main="DI distribution",xlim=c(0,10),breaks=seq(0,10,0.5)) # asymétrique droite
hist(BDD_esp$BT)      # asymétrique droite
hist(BDD_esp$BT_test,xlab="BT",main="BT distribution",ylim=c(0,25)) # asymétrique droite
hist(BDD_esp$BB)      # normale proportion
hist(BDD_esp_net$BB_test,,xlab="BB",main="BB distribution") # normale proportion
hist(BDD_esp_net$MT,xlab="MT",main="MT distribution",xlim=c(100,1000),ylim=c(0,25),breaks=seq(100,1000,100))      # normale 
hist(BDD_esp$BB_prop)
hist(BDD_esp$T_ambiante,xlab="Temperature (°C)",main="T distribution")
hist(BDD_esp$Vent,xlab="Wind (km/h)",main="Wind distribution")
hist(BDD_esp$Humidite,xlab="Humidity (%)",main="Humidity distribution")
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







############################### ACP ######################################

##  Charger les packages nécessaires pour ACP
library(FactoMineR)
library(factoextra)


################  ACP INFLAMMABILITE ####################

# Sélection des colonnes des composantes de l'inflammabilité
colonnes_infla <- BDD_esp[, c(10,12,13,15,32)]

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

n <- BDD_esp[BDD_esp$Nom_scientifique != "Scaevola taccada", ]
BDD_esp_net$score_normalise <- n$score_normalise

HCPC(res.pca,method="complete")


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
  geom_bar() +  
  coord_flip() +
  scale_fill_gradient2(
    low = "darkgreen",
    mid = "yellow",
    high = "red",
    midpoint = 0.48,
    limits = c(scoremin, scoremax),) +
  labs(x = "Espèces",
       y = "Score d'inflammabilité",
       fill = "Score",
       title = "Score d'inflammabilité par espèce"
  ) 

################# sans ggplot ######################
# Palette de couleur
pal_vert_jaune <- colorRampPalette(c("darkgreen", "yellow"))
pal_jaune_rouge <- colorRampPalette(c("yellow", "red"))

# Définir les bornes du graph
mean_score <- mean(BDD_esp$score, na.rm = TRUE)
min_score <- min(BDD_esp$score, na.rm = TRUE)
max_score <- max(BDD_esp$score, na.rm = TRUE)

# Nombre de nuances de chaque côté
n_colors <- 100
n_left <- round((mean_score - min_score) / (max_score - min_score) * n_colors)
n_right <- n_colors - n_left

# Générer  couleurs
couleurs <- c(pal_vert_jaune(n_left), pal_jaune_rouge(n_right))

# Attribuer à chaque espèce une couleur selon sa position par rapport à la min–max
score_scaled <- round((BDD_esp$score - min_score) / (max_score - min_score) * (length(couleurs) - 1)) + 1

# Ordre des points
o <- order(BDD_esp$score)
couleur_points <- couleurs[score_scaled][o]

# Tracé
par(mar = c(4, 9, 0, 0))
plot(BDD_esp$score[o], 1:length(BDD_esp$Nom_scientifique), axes="n")
# Axes
axis(2, at = 1:length(BDD_esp$Nom_scientifique), labels = BDD_esp$Nom_scientifique[o], las = 1, cex.axis = 0.7)
mtext("Espèces", side = 2, line = 8.5, cex = 1)
axis(1)
mtext("Flammability score", side = 1, line = 3, cex = 1)

# Lignes de référence
abline(v = mean_score, lwd = 2)
abline(v = quantile(BDD_esp$score, na.rm = TRUE)[2], lwd = 2, lty = 3)
abline(v = quantile(BDD_esp$score, na.rm = TRUE)[4], lwd = 2, lty = 3)

# Points colorés
points(BDD_esp$score[o], 1:length(BDD_esp$Nom_scientifique), 
       pch = 21, cex = 1.5, bg = couleur_points)




#################### ACP TRAITS ############################

# Sélection des colonnes des traits fonctionnels
colonnes_traits <- na.omit(BDD_esp[, setdiff(17:31, c(23,24, 27))])

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
traits_non_corrélés <- colonnes_traits[, -findCorrelation(mat_cor_trait, cutoff = 0.68)]
traits_non_corrélés

plot(BDD_esp$Gmin~BDD_esp$PEF)








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

boxplot(log(BDD_ech$BT_test+0.01)~BDD_ech$Nom_scientifique)
boxplot(BDD_ech$BT_test~BDD_ech$Nom_scientifique)
log(100)

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


# Calcul des moyennes par espèce
moyennes_FI <- aggregate(FI ~ Nom_scientifique, data = BDD_ech, FUN = mean, na.rm = TRUE)

# Ordonner les espèces croissante
moyennes_FI <- moyennes_FI[order(moyennes_MT$MT),  ]
BDD_ech$Nom_scientifique <- factor(BDD_ech$Nom_scientifique, levels = moyennes_FI$Nom_scientifique)

# Créer un dataframe pour les moyennes (pour ggplot)
moyennes_FI$Type <- "Moyenne"
BDD_ech$Type <- "Individuel"

# Plot avec légende
FI<-ggplot() +
  geom_point(data = BDD_ech, aes(x = FI, y = Nom_scientifique, color = Type), size = 2) +  
  geom_point(data = moyennes_FI, aes(x = FI, y = Nom_scientifique, color = Type), size = 3)+
  scale_color_manual(values = c("Individuel" = "black", "Moyenne" = "#f14900")) +
  labs(x = "FI", y = "Espèce", color = "") +
  theme_bw()
FI













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
vars <- c("BB_test", "BT_test", "MT", "DI_test","FI", "TD", "Surface_F",
          "LDMC", "SD", "LT")

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




kruskal.test(BDD_echMT$MT~BDD_echMT$Nom_scientifique)
kruskal.test(BDD_ech3$BB_test~BDD_ech3$Nom_scientifique)
kruskal.test(BDD_ech3$BT_test~BDD_ech3$Nom_scientifique)
kruskal.test(BDD_ech3$DI_test~BDD_ech3$Nom_scientifique)

kruskal.test(BDD_ech$LDMC~BDD_ech$Nom_scientifique)
kruskal.test(BDD_ech$LMC_t24~BDD_ech$Nom_scientifique)
kruskal.test(BDD_ech$SD~BDD_ech$Nom_scientifique)
kruskal.test(BDD_ech$TD~BDD_ech$Nom_scientifique)
kruskal.test(BDD_ech$Surface_F~BDD_ech$Nom_scientifique)
kruskal.test(BDD_ech$LT~BDD_ech$Nom_scientifique)


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
mod_BB <- lmer(BB_prop ~ 1 + (1 | Nom_scientifique), data = BDD_ech)
summary(mod_BB)

r.squaredGLMM(mod_BT)
r.squaredGLMM(mod_DI)
r.squaredGLMM(mod_MT)
r.squaredGLMM(mod_BB)

# SD
mod_SD <- glmer(SD ~ 1 + (1 | Nom_scientifique),family = "Gamma", data = BDD_ech)
summary(mod_SD)
r.squaredGLMM(mod_SD)

# TD
mod_TD <- lmer(TD ~ 1 + (1 | Nom_scientifique), data = BDD_ech)
summary(mod_TD)
r.squaredGLMM(mod_TD)

# Surface_F
mod_Surface_F <- glmer(Surface_F ~ 1 + (1 | Nom_scientifique),family = "Gamma", data = BDD_ech)
summary(mod_Surface_F)
r.squaredGLMM(mod_Surface_F)

# LDMC : température max (Gaussian)
mod_LDMC <- lmer(LDMC ~ 1 + (1 | Nom_scientifique), data = BDD_ech)
summary(mod_LDMC)
r.squaredGLMM(mod_LDMC)

# SLA
mod_SLA <- lmer(SLA ~ 1 + (1 | Nom_scientifique), data = BDD_ech)
summary(mod_SLA)
r.squaredGLMM(mod_SLA)








################### matrice de corrélation pour observer l'effet des traits sur infla ################
### de tout 
library (corrplot)
mat_cor_all<-cor(colonnes_complet,method="spearman")
corrplot(mat_cor_all)
round(mat_cor_all, 2)   ## afficher les valeurs
corrplot(mat_cor_all, method = "color", tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")














############################ MODELES #####################################
#standardisation des données 
BDD_esp$SD_cr<-scale(BDD_esp$SD)
BDD_esp$TD_cr<-scale(BDD_esp$TD)
BDD_esp$LA_cr<-scale(BDD_esp$Surface_F)
BDD_esp$LDMC_cr<-scale(BDD_esp$LDMC)
BDD_esp$LT_cr<-scale(BDD_esp$LT)
BDD_esp$Vent_cr<-scale(BDD_esp$Vent)
BDD_esp$Temp_cr<-scale(BDD_esp$T_ambiante)
BDD_esp$LMC_t24_cr<-scale(BDD_esp$LMC_t24)
BDD_esp$Nb_rami_cr<-scale(BDD_esp$Nb_rami)

BDD_esp_net_MT$SD_cr<-scale(BDD_esp_net_MT$SD)
BDD_esp_net_MT$TD_cr<-scale(BDD_esp_net_MT$TD)
BDD_esp_net_MT$LA_cr<-scale(BDD_esp_net_MT$Surface_F)
BDD_esp_net_MT$LDMC_cr<-scale(BDD_esp_net_MT$LDMC)
BDD_esp_net_MT$LT_cr<-scale(BDD_esp_net_MT$LT)
BDD_esp_net_MT$Vent_cr<-scale(BDD_esp_net_MT$Vent)
BDD_esp_net_MT$Temp_cr<-scale(BDD_esp_net_MT$T_ambiante)
BDD_esp_net_MT$LMC_t24_cr<-scale(BDD_esp_net_MT$LMC_t24)
BDD_esp_net_MT$Nb_rami_cr<-scale(BDD_esp_net_MT$Nb_rami)

BDD_esp_net3$SD_cr<-scale(BDD_esp_net3$SD)
BDD_esp_net3$TD_cr<-scale(BDD_esp_net3$TD)
BDD_esp_net3$LA_cr<-scale(BDD_esp_net3$Surface_F)
BDD_esp_net3$LDMC_cr<-scale(BDD_esp_net3$LDMC)
BDD_esp_net3$LT_cr<-scale(BDD_esp_net3$LT)
BDD_esp_net3$Vent_cr<-scale(BDD_esp_net3$Vent)
BDD_esp_net3$Temp_cr<-scale(BDD_esp_net3$T_ambiante)
BDD_esp_net3$LMC_t24_cr<-scale(BDD_esp_net3$LMC_t24)
BDD_esp_net3$Nb_rami_cr<-scale(BDD_esp_net3$Nb_rami)

###### FI ########
# Comptage du nombre d'essais par espèce dans BDD_ech
essais_par_espece <- table(BDD_ech$Nom_scientifique)
essais_par_espece

# Création d'une nouvelle colonne Nb_essais dans BDD_esp en s'appuyant sur le nom scientifique
BDD_esp$Nb_essais <- essais_par_espece[BDD_esp$Nom_scientifique]

# Vérification
head(BDD_esp)

plot(BDD_esp$FI,BDD_esp$LT)
plot(BDD_esp$LT,BDD_esp$FI)




mFI<-glm(cbind(Nb_FI,Nb_essais-Nb_FI)~LMC_t24,data=BDD_esp,family="binomial")
summary(mFI)
pred<-predict(mFI,type="response")

ndata<-data.frame(LMC_t24=seq(min(BDD_esp$LMC_t24),max(BDD_esp$LMC_t24),1))
pred<-predict(mFI,type="response",newdata=ndata)

plot(pred~ndata$LMC_t24,type="l",lwd=2,xlab="LMC_t24",ylab="Ignition Frequency",ylim=c(0,1))
points(BDD_esp$LMC_t24,(BDD_esp$Nb_FI/nb_essais))

summary(pred)


mFI1<-glm(cbind(Nb_FI,Nb_essais-Nb_FI)~LDMC,data=BDD_esp,family="binomial")
summary(mFI1)
pred<-predict(mFI1,type="response")

ndata<-data.frame(LDMC=seq(min(BDD_esp$LDMC),max(BDD_esp$LDMC),1))
pred<-predict(mFI1,type="response",newdata=ndata)

plot(pred~ndata$LDMC,type="l",lwd=2,xlab="LDMC",ylab="Ignition Frequency",ylim=c(0,1))
points(BDD_esp$LDMC,(BDD_esp$Nb_FI/nb_essais))

summary(pred)

mFI<-glm(cbind(Nb_FI,Nb_essais-Nb_FI)~LT,data=BDD_esp,family="binomial")
summary(mFI)
pred<-predict(mFI,type="response")

ndata<-data.frame(LT=seq(min(BDD_esp$LT),max(BDD_esp$LT),1))
pred<-predict(mFI,type="response",newdata=ndata)

plot(pred~ndata$LT,type="l",lwd=2,xlab="LT",ylab="Ignition Frequency",ylim=c(0,1))
points(BDD_esp$LT,(BDD_esp$Nb_FI/nb_essais))




mFI1<-glm(cbind(Nb_FI,Nb_essais-Nb_FI)~LMC_t24+LT+Surface_F+TD+SD,data=BDD_esp,family="binomial")
summary(mFI1)


mFI2<-glm(cbind(Nb_FI,Nb_essais-Nb_FI)~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LDMC_cr+LT_cr,data=BDD_esp,family="binomial")
summary(mFI2)

mFI3<-glm(cbind(Nb_FI,Nb_essais-Nb_FI)~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp,family="binomial")
summary(mFI3)


ndata<-data.frame(LMC_t24=seq(min(BDD_esp$LMC_t24),max(BDD_esp$LMC_t24),1))
pred<-predict(mFI3,type="response",newdata=ndata)

plot(pred~ndata$LMC_t24,type="l",lwd=2,xlab="LMC_t24",ylab="Ignition Frequency",ylim=c(0,1))
points(BDD_esp$LMC_t24,(BDD_esp$Nb_FI/nb_essais))









mFI4<-glm(cbind(Nb_FI,Nb_essais-Nb_FI)~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LT_cr,data=BDD_esp,family="binomial")
summary(mFI4)
AIC(mFI2,mFI3,mFI4)

###################### graph bleu rouge ###############

# Extraire coefficients (sans l'intercept)
par(mar = c(5,7,5,5))
coefs <- summary(mFI4)$coefficients[-1, ]

# Variables utiles
estimates <- coefs[, "Estimate"]
stderr <- coefs[, "Std. Error"]
pval <- coefs[, "Pr(>|z|)"]
labels <- rownames(coefs)

# Calcul des intervalles de confiance plus ou moins SE
ci <- cbind(estimates - stderr, estimates + stderr)

# Étoiles de significativité
stars <- ifelse(pval < 0.001, "***",
                ifelse(pval < 0.01, "**",
                       ifelse(pval < 0.05, "*",
                              ifelse(pval < 0.1, ".", ""))))

# Couleurs selon signe du coefficient
cols <- ifelse(estimates < 0, "#e90000", "#117304")

# Ordre des variables (du bas vers le haut)
y_pos <- length(estimates):1

# Plot de base
plot(estimates, y_pos,
     xlim = c(-3,3),
     ylim = c(1,length(estimates) +0.25 ),
     pch = 16, col = cols,
     xlab = "Estimate",
     ylab = "",
     axes = FALSE,
     main = "Ignition Frequency",
     cex.main = 1.5,cex = 1.5)

# Ligne verticale à zéro
abline(v = 0, lty = 2)

# Barres d'erreur (IC)
segments(ci[,1], y_pos, ci[,2], y_pos, col = cols)

# Axe Y avec noms des variables
axis(2, at = y_pos, labels = c("SD", "TD", "LA", "LMC_t24","LT"), las = 1,cex.axis = 1.5)

# Axe X
axis(1,cex.axis = 1.5)

# Valeurs des coefficients + étoiles
text(estimates, y_pos + 0.15,
     labels = paste0(round(estimates, 2), stars),
     col = cols, font = 2, cex = 1.7)

# Ajouter un texte avec le pseudo R² (à modifier si besoin)
mtext(expression(R^2~"= 0.72"), side = 3, adj = 0, line = 0.5, cex = 1.5)









###### MT #########
#plot avec traits (VE)
plot(BDD_esp$TD, BDD_esp$MT)
plot(BDD_esp$LT, BDD_esp$MT)
plot(BDD_esp$LDMC, BDD_esp$MT)
plot(BDD_esp$Surface_F, BDD_esp$MT)
plot(BDD_esp$SD, BDD_esp$MT)

#modèles
m_MT0<-glm(MT~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp_net,family="gaussian")
summary(m_MT0)
m_MT1<-glm(MT~SD_cr+LMC_t24_cr+LDMC_cr,data=BDD_esp,family="gaussian")
summary(m_MT1)
m_MT1<-glm(MT~SD_cr+LMC_t24_cr+LDMC_cr,data=BDD_esp,family="gaussian")
summary(m_MT1)

m_MT2<-glm(MT~SD_cr+TD_cr+LA_cr+Nb_rami_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")
summary(m_MT2)

AIC(m_MT0,m_MT1)
anova(m_MT0,m_MT2)

###################### graph coefs ###############

# Extraire coefficients (sans l'intercept)
par(mar = c(5,7,5,5))
coefs <- summary(m_MT0)$coefficients[-1, ]

# Variables utiles
estimates <- coefs[, "Estimate"]
stderr <- coefs[, "Std. Error"]
pval <- coefs[, "Pr(>|t|)"]
labels <- rownames(coefs)

# Calcul des intervalles de confiance plus ou moins SE
ci <- cbind(estimates - stderr, estimates + stderr)

# Étoiles de significativité
stars <- ifelse(pval < 0.001, "***",
                ifelse(pval < 0.01, "**",
                       ifelse(pval < 0.05, "*",
                              ifelse(pval < 0.1, ".", ""))))

# Couleurs selon signe du coefficient
cols <- ifelse(estimates < 0, "#e90000", "#117304")

# Ordre des variables (du bas vers le haut)
y_pos <- length(estimates):1

# Plot de base
plot(estimates, y_pos,
     xlim = c(-200,200),
     ylim = c(1,length(estimates) +0.25 ),
     pch = 16, col = cols,
     xlab = "Estimate",
     ylab = "",
     axes = FALSE,
     main = "Maximum Temperature",
     cex.main = 1.5,cex = 1.5)

# Ligne verticale à zéro
abline(v = 0, lty = 2)

# Barres d'erreur (IC)
segments(ci[,1], y_pos, ci[,2], y_pos, col = cols)

# Axe Y avec noms des variables
axis(2, at = y_pos, labels = c("SD", "TD", "LA", "LMC_t24","LT"), las = 1,cex.axis = 1.5)

# Axe X
axis(1,cex.axis = 1.5)

# Valeurs des coefficients + étoiles
text(estimates, y_pos + 0.15,
     labels = paste0(round(estimates, 2), stars),
     col = cols, font = 2, cex =  1.7)

# Ajouter un texte avec le pseudo R² (à modifier si besoin)
mtext(expression(R^2~"= 0.58"), side = 3, adj = 0, line = 0.5, cex = 1.5)







m_MT01<-glm(MT~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")
summary(m_MT01)
plot(m_MT01)

m_MT1<-glm(MT~SD_cr+TD_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")        # meilleur modèle
summary(m_MT1)

m_MT2<-glm(MT~SD_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")
summary(m_MT2)

m_MT3<-glm(MT~SD_cr+LDMC_cr,data=BDD_esp,family="gaussian")
summary(m_MT3)

AIC(m_MT0,m_MT01,m_MT1,m_MT2,m_MT3) 




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

plot(BDD_esp$LMC_t24~BDD_esp$Gmin)
lines(pred1$fit~pred_data$Gmin,lwd=2)
lines(pred2$fit~pred_data$Gmin,lwd=2)
lines(pred3$fit~pred_data$Gmin,lwd=2)


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
hist(BDD_esp$BB_prop)
BDD_esp$BB_prop <- BDD_esp$BB_test/100
BDD_esp_net$BB_prop <- BDD_esp_net$BB_test/100


m_BB0<-glm(BB_prop~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp_net,family="gaussian")
summary(m_BB0)
m_BB1<-glm(BB_prop~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")
summary(m_BB1)
m_BB2<-glm(BB_prop~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")
summary(m_BB2)

AIC (m_BB0,m_BB1,m_BB2)

###################### graph coefs ###############

# Extraire coefficients (sans l'intercept)
par(mar = c(5,7,5,5))
coefs <- summary(m_BB0)$coefficients[-1, ]

# Variables utiles
estimates <- coefs[, "Estimate"]
stderr <- coefs[, "Std. Error"]
pval <- coefs[, "Pr(>|t|)"]
labels <- rownames(coefs)

# Calcul des intervalles de confiance plus ou moins SE
ci <- cbind(estimates - stderr, estimates + stderr)

# Étoiles de significativité
stars <- ifelse(pval < 0.001, "***",
                ifelse(pval < 0.01, "**",
                       ifelse(pval < 0.05, "*",
                              ifelse(pval < 0.1, ".", ""))))

# Couleurs selon signe du coefficient
cols <- ifelse(estimates < 0, "#e90000", "#117304")

# Ordre des variables (du bas vers le haut)
y_pos <- length(estimates):1

# Plot de base
plot(estimates, y_pos,
     xlim = c(-0.3,0.3),
     ylim = c(1,length(estimates) +0.25 ),
     pch = 16, col = cols,
     xlab = "Estimate",
     ylab = "",
     axes = FALSE,
     main = "Burnt Biomass",
     cex.main = 1.5,cex = 1.5)

# Ligne verticale à zéro
abline(v = 0, lty = 2)

# Barres d'erreur (IC)
segments(ci[,1], y_pos, ci[,2], y_pos, col = cols)

# Axe Y avec noms des variables
axis(2, at = y_pos, labels = c("SD", "TD", "LA", "LMC_t24","LT"), las = 1,cex.axis = 1.5)

# Axe X
axis(1,cex.axis = 1.5)

# Valeurs des coefficients + étoiles
text(estimates, y_pos + 0.15,
     labels = paste0(round(estimates, 2), stars),
     col = cols, font = 2, cex =  1.7)

# Ajouter un texte avec le pseudo R² (à modifier si besoin)
mtext(expression(R^2~"= 0.46"), side = 3, adj = 0, line = 0.5, cex = 1.5)




m_BB0<-glm(BB_prop~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr+Vent_cr+Temp_cr,data=BDD_esp,family="quasibinomial")
summary(m_BB0)

m_BB01<-glm(BB_prop~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp,family="binomial")
summary(m_BB0)

m_BB1<-glm(BB_prop~LDMC_cr,data=BDD_esp,family="binomial")        # meilleur modèle
summary(m_BB1)

library(betareg)
model <- betareg(BB_prop~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr+Vent_cr+Temp_cr,data=BDD_esp)
summary(model)


AIC(m_BB0,model) #!!! pas le même nombre d'observation car SD a deux esp de moins 



m_BB0<-glm(BB_prop~LDMC+LT,data=BDD_esp,family="gaussian")
summary(m_BB0)
m_BB1<-glm(BB_prop~LDMC+LT,data=BDD_esp,family="quasibinomial")
summary(m_BB1)
m_BB2<-glm(BB_prop~LDMC+LT,data=BDD_esp,family="binomial")
summary(m_BB2)
library(betareg)
model <- betareg(BB_prop~LDMC+LT,data=BDD_esp)
summary(model)

c<-seq(1,1000,0.1)

new<-data.frame(LDMC=c)

pred1<-predict(m_BB0,newdata=new,type = "response")
plot(BDD_esp$BB_prop,pred1)
pred2<-predict(m_BB1,type = "response",newdata=new)
plot(BDD_esp$BB_prop,pred2)
pred3<-predict(m_BB2,type = "response",newdata=new)

pred1
plot(pred1~c,type="l",lwd=3)
lines(pred2~c,lwd=3,lty=3)
lines(pred3~c,lwd=3,lty=2)
points(BDD_esp$BB_prop~BDD_esp$LDMC)





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
m_BT1<-glm(BT_test~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LT_cr,data=BDD_esp,family="Gamma")
summary(m_BT1)
m_BT2<-glm(BT_test~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LDMC_cr+LT_cr,data=BDD_esp,family="Gamma")
summary(m_BT2)

AIC(m_BT0,m_BT1,m_BT2)

###################### graph coefs ###############

# Extraire coefficients (sans l'intercept)
par(mar = c(5,7,5,5))
coefs <- summary(m_BT0)$coefficients[-1, ]

# Variables utiles
estimates <- coefs[, "Estimate"]
stderr <- coefs[, "Std. Error"]
pval <- coefs[, "Pr(>|t|)"]
labels <- rownames(coefs)

# Calcul des intervalles de confiance plus ou moins SE
ci <- cbind(estimates - stderr, estimates + stderr)

# Étoiles de significativité
stars <- ifelse(pval < 0.001, "***",
                ifelse(pval < 0.01, "**",
                       ifelse(pval < 0.05, "*",
                              ifelse(pval < 0.1, ".", ""))))

# Couleurs selon signe du coefficient
cols <- ifelse(estimates < 0, "#e90000", "#117304")

# Ordre des variables (du bas vers le haut)
y_pos <- length(estimates):1

# Plot de base
plot(estimates, y_pos,
     xlim = c(-0.015,0.015),
     ylim = c(1,length(estimates) +0.25 ),
     pch = 16, col = cols,
     xlab = "Estimate",
     ylab = "",
     axes = FALSE,
     main = "Burning Time",
     cex.main = 1.5,cex = 1.5)

# Ligne verticale à zéro
abline(v = 0, lty = 2)

# Barres d'erreur (IC)
segments(ci[,1], y_pos, ci[,2], y_pos, col = cols)

# Axe Y avec noms des variables
axis(2, at = y_pos, labels = c("SD", "TD", "LA", "LMC_t24","LT"), las = 1,cex.axis = 1.5)

# Axe X
axis(1,cex.axis = 1.5)

# Valeurs des coefficients + étoiles
text(estimates, y_pos + 0.15,
     labels = paste0(round(estimates, 4), stars),
     col = cols, font = 2, cex =  1.7)

# Ajouter un texte avec le pseudo R² (à modifier si besoin)
mtext(expression(R^2~"= 0.51"), side = 3, adj = 0, line = 0.5, cex = 1.5)




m_BT1<-glm(BT_test~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp,family="Gamma")       
summary(m_BT1)

m_BT2<-glm(BT_test~SD_cr+TD_cr+LA_cr+LDMC_cr,data=BDD_esp,family="Gamma")    ### meilleur modèle 
summary(m_BT2)

m_BT3<-glm(BT_test~SD_cr+TD_cr+LDMC_cr,data=BDD_esp,family="Gamma")
summary(m_BT3)


AIC(m_BT0,m_BT1,m_BT2,m_BT3,m_BT4) #!!! pas le même nombre d'observation car SD a deux esp de moins 


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
m_DI1<-glm(DI_test~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LT_cr,data=BDD_esp,family="Gamma")
summary(m_DI1)
m_DI2<-glm(DI_test~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LDMC_cr+LT_cr,data=BDD_esp,family="Gamma")
summary(m_DI2)

AIC(m_DI0,m_DI1,m_DI2)

###################### graph coefs ###############

# Extraire coefficients (sans l'intercept)
par(mar = c(5,7,5,5))
coefs <- summary(m_DI0)$coefficients[-1, ]

# Variables utiles
estimates <- coefs[, "Estimate"]
stderr <- coefs[, "Std. Error"]
pval <- coefs[, "Pr(>|t|)"]
labels <- rownames(coefs)

# Calcul des intervalles de confiance plus ou moins SE
ci <- cbind(estimates - stderr, estimates + stderr)

# Étoiles de significativité
stars <- ifelse(pval < 0.001, "***",
                ifelse(pval < 0.01, "**",
                       ifelse(pval < 0.05, "*",
                              ifelse(pval < 0.1, ".", ""))))

# Couleurs selon signe du coefficient
cols <- ifelse(estimates < 0, "#e90000", "#117304")

# Ordre des variables (du bas vers le haut)
y_pos <- length(estimates):1

# Plot de base
plot(estimates, y_pos,
     xlim = c(-0.6,0.6),
     ylim = c(1,length(estimates) +0.25 ),
     pch = 16, col = cols,
     xlab = "Estimate",
     ylab = "",
     axes = FALSE,
     main = "Ignition Delay",
     cex.main = 1.5,cex = 1.5)

# Ligne verticale à zéro
abline(v = 0, lty = 2)

# Barres d'erreur (IC)
segments(ci[,1], y_pos, ci[,2], y_pos, col = cols)

# Axe Y avec noms des variables
axis(2, at = y_pos, labels = c("SD", "TD", "LA", "LMC_t24","LT"), las = 1,cex.axis = 1.5)

# Axe X
axis(1,cex.axis = 1.5)

# Valeurs des coefficients + étoiles
text(estimates, y_pos + 0.15,
     labels = paste0(round(estimates, 4), stars),
     col = cols, font = 2, cex =  1.7)

# Ajouter un texte avec le pseudo R² (à modifier si besoin)
mtext(expression(R^2~"= 0.62"), side = 3, adj = 0, line = 0.5, cex = 1.5)



m_DI1<-glm(DI_test~TD_cr+LDMC_cr+LA_cr+LT_cr+Vent_cr+Temp_cr,data=BDD_esp,family="Gamma")       
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
m_score1<-glm(score_normalise~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LT_cr,data=BDD_esp,family="gaussian")
summary(m_score1)
m_score2<-glm(score_normalise~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")
summary(m_score2)


AIC(m_score0,m_score1,m_score2)


###################### graph coefs ###############

# Extraire coefficients (sans l'intercept)
par(mar = c(5,7,5,5))
coefs <- summary(m_score1)$coefficients[-1, ]

# Variables utiles
estimates <- coefs[, "Estimate"]
stderr <- coefs[, "Std. Error"]
pval <- coefs[, "Pr(>|t|)"]
labels <- rownames(coefs)

# Calcul des intervalles de confiance plus ou moins SE
ci <- cbind(estimates - stderr, estimates + stderr)

# Étoiles de significativité
stars <- ifelse(pval < 0.001, "***",
                ifelse(pval < 0.01, "**",
                       ifelse(pval < 0.05, "*",
                              ifelse(pval < 0.1, ".", ""))))

# Couleurs selon signe du coefficient
cols <- ifelse(estimates < 0, "#e90000", "#117304")

# Ordre des variables (du bas vers le haut)
y_pos <- length(estimates):1

# Plot de base
plot(estimates, y_pos,
     xlim = c(-0.6,0.6),
     ylim = c(1,length(estimates) +0.25 ),
     pch = 16, col = cols,
     xlab = "Estimate",
     ylab = "",
     axes = FALSE,
     main = "Flammability score",
     cex.main = 1.5,cex = 1.5)

# Ligne verticale à zéro
abline(v = 0, lty = 2)

# Barres d'erreur (IC)
segments(ci[,1], y_pos, ci[,2], y_pos, col = cols)

# Axe Y avec noms des variables
axis(2, at = y_pos, labels = c("SD", "TD", "LA", "LMC_t24","LT"), las = 1,cex.axis = 1.5)

# Axe X
axis(1,cex.axis = 1.5)

# Valeurs des coefficients + étoiles
text(estimates, y_pos + 0.15,
     labels = paste0(round(estimates, 4), stars),
     col = cols, font = 2, cex =  1.7)

# Ajouter un texte avec le pseudo R² (à modifier si besoin)
mtext(expression(R^2~"= 0.75"), side = 3, adj = 0, line = 0.5, cex = 1.5)




m_score1<-glm(score_normalise~SD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")       
summary(m_score1)

m_score2<-glm(score_normalise~SD_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")       
summary(m_score2)

m_score3<-glm(score_normalise~LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")  ##meilleur modèle     
summary(m_score3)

m_score4<-glm(score_normalise~LDMC_cr,data=BDD_esp,family="gaussian")       
summary(m_score4)

AIC(m_score4,m_score3)































############################################################## ECHELLE ESPECE #############################################################
BDD_ech <- BDD_ech[!is.na(BDD_ech$MT), ]
################## DISTRIBUTION DES DONNEES #####################
#Infla
hist(BDD_ech$DI)      # asymétrique gauche
hist(BDD_ech$DI_test) # asymétrique droite
hist(BDD_ech$BT)      # asymétrique droite
hist(BDD_ech$BT_test) # asymétrique droite
hist(BDD_ech$BB)      # normale proportion
hist(BDD_ech$BB_test) # normale proportion
hist(BDD_ech$MT)      # normale 

#traits
hist(BDD_ech$Nb_ramifications)   #asymétrique droite
hist(BDD_ech$SD)                 #asymétrique droite
hist(BDD_ech$TMC_t0)             #asymétrique droite ou normale (hésitation)
hist(BDD_ech$TMC_t24)            #asymétrique droite ou normale (hésitation)
hist(BDD_ech$TDMC)               #normale
hist(BDD_ech$TD)                 #normale
hist(BDD_ech$TDIA)               #normale
hist(BDD_ech$LMC_t0)             #asymétrique droite ou normale (hésitation) 
hist(BDD_ech$LMC_t24)            #asymétrique droite ou normale (hésitation) 
hist(BDD_ech$LDMC)               #normale
hist(BDD_ech$Surface_F)          #asymétrique droite
hist(BDD_ech$SLA)                #normale
hist(BDD_ech$LT)                 #asymétrique droite




################ EFFET DES CONDITIONS METEO ##################
#MT
options(na.action = "na.omit")
m<-glm(MT~Vent,data=BDD_ech,family="gaussian")
summary(m)

#BB 
BDD_ech$BB_prop <- BDD_ech$BB_test/100
m<-glm(BB_prop~Vent,data=BDD_ech,family="binomial")
summary(m)

#DI
m<-glm(DI_test~T_ambiante+Vent,data=BDD_ech,family="Gamma")
summary(m)
plot(BDD_ech$T_ambiante,1/BDD_ech$DI_test)

#BT
m<-glm(BT_test~T_ambiante,data=BDD_ech,family="Gamma")
summary(m)
plot(BDD_ech$Vent,BDD_ech$BT_test)



View(BDD_ech)


############################### ACP ######################################

##  Charger les packages nécessaires pour ACP
library(FactoMineR)
library(factoextra)


################  ACP INFLAMMABILITE ####################

# Sélection des colonnes des composantes de l'inflammabilité
colonnes_infla_ech <- BDD_ech[, c(10,12,13,15)]

# Vérification des données
head(colonnes_infla_ech)

# Centrage-réduction des données
colonnes_infla_ech_cr <- scale(colonnes_infla_ech)

# Vérification des données standardisées
head(colonnes_infla_ech_cr)

# Application de l'ACP
res.pca <- PCA(colonnes_infla_ech_cr, scale.unit = FALSE, graph = FALSE)

# Résumé des résultats
summary(res.pca)


# Graphique des contributions des variables aux composantes principales
fviz_pca_var(res.pca, col.var = "contrib", gradient.cols = c("#6fec00", "#ff9e00", "#Ff0000"),repel = TRUE)


# Graphique combiné des variables et des individus
fviz_pca_biplot(res.pca,col.var = "contrib", gradient.cols =c("#00AFBB", "#E7B800", "#FC4E07") , repel = TRUE,label = "var")

# Afficher l'ébouli
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50), main="Graphique de l'ébouli")

# Ajout d'un score d'inflammabilité basé sur la coordonnée de l'axe 1 
coord <- res.pca$ind$coord  # coordonnées des individus
BDD_ech$score <- coord[,1]


# normalisation du score d'inflammabilité entre -1 et 1
min_score <- min(BDD_ech$score, na.rm = TRUE)
max_score <- max(BDD_ech$score, na.rm = TRUE)
BDD_ech$score_normalise <- -1 + (BDD_ech$score - min_score) * 2 / (max_score - min_score)


HCPC(res.pca,method="ward")




####################### CLASSEMENT ESPECES ############################

################ graphique pour visualiser les espèces et leur inflammabilité ####################
library(ggplot2)


# 2. Calcul des moyennes de score par espèce
BDD_moyennes <- aggregate(score_normalise ~ Nom_scientifique, data = BDD_ech, FUN = mean, na.rm = TRUE)
BDD_moyennes
# 3. Définir les bornes du gradient
scoremin <- min(BDD_moyennes$score_normalise)
scoremax <- max(BDD_moyennes$score_normalise)

# 4. Graphique
ggplot(BDD_moyennes, aes(x = reorder(Nom_scientifique, score_normalise), y = score_normalise, fill = score_normalise)) +
  geom_bar(stat = "identity", color = NA) +  
  coord_flip() +
  scale_fill_gradient2(
    low = "darkgreen",
    mid = "yellow",
    high = "red",
    midpoint = 0,
    limits = c(scoremin, scoremax)
  ) +
  labs(
    x = "Espèces",
    y = "Score d'inflammabilité (moyenne)",
    fill = "Score moyen",
    title = "Score d'inflammabilité moyen par espèce"
  ) +
  theme_minimal()


################# Heatmap (couleur en fonction de chaque composante d'inflammabilité) ##################

# packages nécessaires
library(ggplot2)
library(tidyr)  # pour changer le format de la BDD en format "long"

# Retirer lignes NA MT
df_prep <- BDD_ech[!is.na(BDD_ech$MT), ]

# Fonction normalisation
normalize <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

# Normalisation sur les données brutes (pas sur les données déjà normalisées)
df_norm <- df_prep
df_norm$MT <- normalize(df_prep$MT)
df_norm$score_DI <- normalize(df_prep$score_DI)
df_norm$BB_test <- normalize(df_prep$BB_test)
df_norm$BT_test <- normalize(df_prep$BT_test)

# Moyenne par espèce
MT_moy <- aggregate(MT ~ Nom_scientifique, data = df_norm, FUN = mean, na.rm = TRUE)
DI_moy <- aggregate(score_DI ~ Nom_scientifique, data = df_norm, FUN = mean, na.rm = TRUE)
BB_moy <- aggregate(BB_test ~ Nom_scientifique, data = df_norm, FUN = mean, na.rm = TRUE)
BT_moy <- aggregate(BT_test ~ Nom_scientifique, data = df_norm, FUN = mean, na.rm = TRUE)

df_moy <- merge(MT_moy, DI_moy, by = "Nom_scientifique")
df_moy <- merge(df_moy, BB_moy, by = "Nom_scientifique")
df_moy <- merge(df_moy, BT_moy, by = "Nom_scientifique")

# Normalisation par variable des moyennes pour que toutes les variables soient comparables en couleurs
vars <- c("MT", "score_DI", "BB_test", "BT_test")
for (v in vars) {
  rng <- range(df_moy[[v]], na.rm = TRUE)
  if (diff(rng) == 0) {
    df_moy[[v]] <- 0
  } else {
    df_moy[[v]] <- (df_moy[[v]] - rng[1]) / diff(rng)
  }
}

# Ordre des espèces par moyenne globale (score_normalise)
score_moy <- aggregate(score_normalise ~ Nom_scientifique, data = df_norm, FUN = mean, na.rm = TRUE)
ordre_esp <- score_moy[order(score_moy$score_normalise), "Nom_scientifique"]

library(tidyr)
df_long <- pivot_longer(df_moy,
                        cols = vars,
                        names_to = "Variable",
                        values_to = "Valeur")

df_long$Nom_scientifique <- factor(df_long$Nom_scientifique, levels = rev(ordre_esp))

library(ggplot2)
ggplot(df_long, aes(x = Variable, y = Nom_scientifique, fill = Valeur)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(
    colors = c("darkgreen", "yellow", "red"),
    breaks = c(0.1, 0.5, 0.9),
    labels = c("Faible", "Moyenne", "Élevée"),
    name = "Inflammabilité (moyenne normalisée)"
  ) +
  labs(title = "Heatmap des composantes d’inflammabilité (moyennes par espèce)",
       x = "Composantes",
       y = "Espèces") +
  scale_y_discrete(drop = FALSE) +
  theme(axis.text.y = element_text(size = 6))




#################### ACP TRAITS ############################

# Sélection des colonnes des traits fonctionnels
colonnes_traits_ech <- na.omit(BDD_ech[, setdiff(16:30, c(23, 26))])

# Vérifier les données
colonnes_traits_ech


colonnes_traits_ech_cr <- scale(log(colonnes_traits_ech+0.01))
colonnes_traits_ech_cr

summary(colonnes_traits_ech_cr)
any(is.na(colonnes_traits_ech_cr))  # TRUE si NA
any(is.infinite(as.matrix(colonnes_traits_ech_cr)))  # TRUE si Inf ou -Inf


# Application de l'ACP
res.pca <- PCA(colonnes_traits_ech_cr, scale.unit = FALSE, graph = FALSE)

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
mat_cor_trait_ech<-cor(colonnes_traits_ech,method="spearman")
corrplot(mat_cor_trait_ech)
round(mat_cor_trait_ech, 2)   ## afficher les valeurs
# Visualiser la corrélation
corrplot(mat_cor_trait_ech, method = "color", type = "upper", tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")
corrplot(mat_cor_trait_ech, method = "color", tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")

library(caret)
traits_non_corrélés_ech <- colonnes_traits_ech[, -findCorrelation(mat_cor_trait_ech, cutoff = 0.7)]
traits_non_corrélés_ech


plot(log(BDD_ech$LDMC),log(BDD_ech$LMC_t24))
plot(BDD_ech$MT~BDD_ech$LDMC)
plot(BDD_ech$MT~BDD_ech$LMC_t24,xlim=(c(0,500)))
plot(BDD_ech$MT,BDD_ech$LMC_t24)

summary(lm(BDD_ech$MT~BDD_ech$LDMC))
summary(lm(BDD_ech$MT~BDD_ech$LMC_t24))
plot(BDD_ech$LMC_t24,BDD_ech$LMC_t0)


######## ACP INFLA avec projection des axes TRAITS (TEST) ####################

# Sélection des lignes complètes (pas de NA dans les colonnes 3 à 22)
colonnes_complet_ech_test <- na.omit(BDD_ech[, setdiff(8:30, c(8,9,11,14,23, 26))])
dim (colonnes_complet_ech_test)

colonnes_complet_ech_test_cr <- scale(colonnes_complet_ech_test)
head(colonnes_complet_ech_test_cr)
## ACP ##
# colonnes 1 à 4 : inflammabilité (axes actifs)
# colonnes 5 à 20 : traits (axes supplémentaires)
res.pca <- PCA(colonnes_complet_ech_test_cr, scale.unit = FALSE, 
               quanti.sup = 5:17, graph = FALSE)

# Visualisation ACP
fviz_pca_var(res.pca, col.var = "red", repel = TRUE) +
  ggtitle("Variables actives (infla) et supplémentaires (traits)")

# ACP avec les individus
fviz_pca_biplot(res.pca, col.var = "red",repel = TRUE)








############################ MODELES #####################################
library(lmerTest)
#standardisation des données 
BDD_ech$SD_cr<-scale(BDD_ech$SD)
BDD_ech$TD_cr<-scale(BDD_ech$TD)
BDD_ech$LA_cr<-scale(BDD_ech$Surface_F)
BDD_ech$LDMC_cr<-scale(BDD_ech$LDMC)
BDD_ech$LT_cr<-scale(BDD_ech$LT)
BDD_ech$Vent_cr<-scale(BDD_ech$Vent)
BDD_ech$Temp_cr<-scale(BDD_ech$T_ambiante)

###### MT #########
#plot avec traits (VE)
plot(BDD_ech$TD, BDD_ech$MT)
plot(BDD_ech$LT, BDD_ech$MT)
plot(BDD_ech$LDMC, BDD_ech$MT)
plot(BDD_ech$Surface_F, BDD_ech$MT)
plot(BDD_ech$SD, BDD_ech$MT)

#modèles
modech_MT0 <- glm(MT ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr, family="gaussian",data = BDD_ech)
summary(modech_MT0)

modech_MT1 <- lmer(MT ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr + (1 | Nom_scientifique), data = BDD_ech)
summary(modech_MT1)
r.squaredGLMM(modech_MT1)

plot(modech_MT1)

AIC(modech_MT0,modech_MT1) 



###### BB #########
#plot avec traits (VE)
plot(BDD_ech$TD, BDD_ech$BB)
plot(BDD_ech$LT, BDD_ech$BB)
plot(BDD_ech$LDMC, BDD_ech$BB)
plot(BDD_ech$Surface_F, BDD_ech$BB)
plot(BDD_ech$SD, BDD_ech$BB)

#modèles
BDD_ech$BB_prop <- BDD_ech$BB_test/100
modech_BB0 <- glm(BB_prop ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,family = "gaussian", data = BDD_ech)
summary(modech_BB0)

modech_BB4 <- glm(BB_prop ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,family = "binomial", data = BDD_ech)
summary(modech_BB4)

modech_BB1 <- lmer(BB_prop ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr + (1 | Nom_scientifique), data = BDD_ech)    #### meilleur modèle 
summary(modech_BB1)

modech_BB2 <- glmer(BB_prop ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr + (1 | Nom_scientifique), family = "binomial", data = BDD_ech)
summary(modech_BB2)

AIC(modech_BB0, modech_BB1,modech_BB2,modech_BB4) 




###### BT #########
#plot avec traits (VE)
plot(BDD_ech$TD, BDD_ech$BT_test)
plot(BDD_ech$LT, BDD_ech$BT_test)
plot(BDD_ech$LDMC, BDD_ech$BT_test)
plot(BDD_ech$Surface_F, BDD_ech$BT_test)
plot(BDD_ech$SD, BDD_ech$BT_test)

#modèles

modech_BT0 <- glm(BT_test ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,family = "Gamma", data = BDD_ech)
summary(modech_BT0)

modech_BT2 <- glmer(BT_test ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr + (1 | Nom_scientifique), family = "Gamma", data = BDD_ech)  #meilleur modèle
summary(modech_BT2)

AIC(modech_BT0,modech_BT2) 



###### DI #########
#plot avec traits (VE)
plot(BDD_ech$TD, BDD_ech$DI_test)
plot(BDD_ech$LT, BDD_ech$DI_test)
plot(BDD_ech$LDMC, BDD_ech$DI_test)
plot(BDD_ech$Surface_F, BDD_ech$DI_test)
plot(BDD_ech$SD, BDD_ech$DI_test)

#modèles

modech_DI0 <- glm(DI_test ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,family = "Gamma", data = BDD_ech)
summary(modech_DI0)

modech_DI2 <- glmer(DI_test ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr + (1 | Nom_scientifique), family = "Gamma", data = BDD_ech)  #meilleur modèle
summary(modech_DI2)

AIC(modech_DI0,modech_DI2) 



###### Score ########
hist(BDD_ech$score_normalise)
#plot avec traits (VE)
plot(BDD_ech$TD, BDD_ech$score_normalise)
plot(BDD_ech$LT, BDD_ech$score_normalise)
plot(BDD_ech$LDMC, BDD_ech$score_normalise)
plot(BDD_ech$Surface_F, BDD_ech$score_normalise)
plot(BDD_ech$SD, BDD_ech$score_normalise)

#modèles
modech_score_normalise0 <- glm(score_normalise ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr, family="gaussian",data = BDD_ech)
summary(modech_score_normalise0)

modech_score_normalise1 <- lmer(score_normalise ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr + (1 | Nom_scientifique), data = BDD_ech)
summary(modech_score_normalise1)

AIC(modech_score_normalise0,modech_score_normalise1) 
























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

###################

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
                       LDMC=rep(mean(BDD_esp$LDMC),length(foo)))

pred_data2<-data.frame(Gmin=foo,
                       LDMC=rep(quantile(BDD_esp$LDMC,p=0.75),length(foo)))

pred_data3<-data.frame(Gmin=foo,
                       LDMC=rep(quantile(BDD_esp$LDMC,p=0.05),length(foo)))


pred1<-predict(m,type="response",newdata=pred_data1,se.fit=TRUE)
pred2<-predict(m,type="response",newdata=pred_data2,se.fit=TRUE)
pred3<-predict(m,type="response",newdata=pred_data3,se.fit=TRUE)

plot(BDD_esp$LMC_t24~BDD_esp$Gmin)
lines(pred1$fit~pred_data$Gmin,lwd=2)
lines(pred2$fit~pred_data$Gmin,lwd=2)
lines(pred3$fit~pred_data$Gmin,lwd=2)



lines(pred$fit~pred_data$Gmin,lwd=2)
lines(pred$fit+(1.96*pred$se.fit)~pred_data$Gmin,lty=3)
lines(pred$fit-(1.96*pred$se.fit)~pred_data$Gmin,lty=3)


hist(BDD_esp$LMC_t24)


mFI<-glm(cbind(Nb_FI,nb_essais-Nb_FI)~LMC_t24,data=BDD_esp,family="binomial")
summary(mFI)

plot(BDD_esp$Nb_FI~BDD_esp$LMC_t24)







# Modèle
m <- glm(LMC_t24 ~ Gmin + LDMC, data = BDD_esp, family = "Gamma")

# Valeurs de Gmin
foo <- seq(min(BDD_esp$Gmin), max(BDD_esp$Gmin), length.out = 100)

# Trois niveaux de LDMC : moyen, haut (Q75), bas (Q5)
pred_data1 <- data.frame(Gmin = foo,
                         LDMC = rep(mean(BDD_esp$LDMC), length(foo)))

pred_data2 <- data.frame(Gmin = foo,
                         LDMC = rep(quantile(BDD_esp$LDMC, 0.75), length(foo)))

pred_data3 <- data.frame(Gmin = foo,
                         LDMC = rep(quantile(BDD_esp$LDMC, 0.05), length(foo)))

# Prédictions
pred1 <- predict(m, type = "response", newdata = pred_data1, se.fit = TRUE)
pred2 <- predict(m, type = "response", newdata = pred_data2, se.fit = TRUE)
pred3 <- predict(m, type = "response", newdata = pred_data3, se.fit = TRUE)

# Plot des points observés
plot(BDD_esp$Gmin, BDD_esp$LMC_t24, 
     xlab = "Gmin", ylab = "LMC_t24", 
     main = "Effet de Gmin sur LMC_t24 selon LDMC")

# Courbes de prédiction
lines(foo, pred1$fit, col = "blue", lwd = 2)  # LDMC moyen
lines(foo, pred2$fit, col = "darkgreen", lwd = 2)  # LDMC haut
lines(foo, pred3$fit, col = "orange", lwd = 2)  # LDMC bas

# Intervalle de confiance pour LDMC moyen
lines(foo, pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines(foo, pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)

# Légende
legend("topright", legend = c("LDMC moyen", "LDMC Q75", "LDMC Q5"),
       col = c("blue", "darkgreen", "orange"), lwd = 2, bty = "n")

























################################ graphs prédiction ##############################


vars <- c("SD_cr", "TD_cr", "LA_cr", "LMC_t24_cr", "LT_cr")
BDD_esp[vars] <- lapply(BDD_esp[vars], function(x) as.numeric(x))

############ FI ######################
# Modèle GLM binomial
mFI3 <- glm(cbind(Nb_FI, Nb_essais - Nb_FI) ~ SD_cr + TD_cr + LA_cr + LMC_t24_cr + LT_cr,data = BDD_esp, family = "binomial")

summary(mFI3)

moy <- mean(BDD_esp$LMC_t24, na.rm = TRUE)
ecart <- sd(BDD_esp$LMC_t24, na.rm = TRUE)

# Valeurs de LMC_t24
foo <- seq(min(BDD_esp$LMC_t24_cr), max(BDD_esp$LMC_t24_cr),length.out = 100)

# on fait varier LMDC et on fixe les autres variables
pred_data1 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_esp$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_esp$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp$LA_cr), length(foo)),
  LT_cr = rep(mean(BDD_esp$LT_cr), length(foo))
)


# Prédictions
pred1 <- predict(mFI3, type = "response", newdata = pred_data1, se.fit = TRUE)


# Plot des points observés
plot(BDD_esp$LMC_t24, BDD_esp$Nb_FI / BDD_esp$Nb_essais, 
     xlab = "LMC_t24", ylab = "Ignition Frequency",ylim = c(0,1), 
     main = "Effet de LMC_t24 sur LMC_t24 selon SD")

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "blue", lwd = 2) 

axis(2)

# Intervalle de confiance pour SD moyen
lines((foo*ecart+moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart+moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)



############ MT ######################
# Modèle GLM 
m_MT0<-glm(MT~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp_net,family="gaussian")
summary(m_MT0)
moy <- mean(BDD_esp$LDMC_cr, na.rm = TRUE)
ecart <- sd(BDD_esp$LDMC_cr, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_esp_net$LDMC_cr, na.rm = TRUE), max(BDD_esp_net$LDMC_cr,na.rm = TRUE))

# on fait varier 1 variable et on fixe les autres variables
pred_data1 <- data.frame(
  LDMC = foo,
  SD = rep(mean(BDD_esp_net$SD,na.rm = TRUE), length(foo)),
  TD = rep(mean(BDD_esp_net$TD), length(foo)),
  Surface_F = rep(mean(BDD_esp_net$Surface_F), length(foo)),
  LT = rep(mean(BDD_esp_net$LT), length(foo))
)


# Prédictions
pred1 <- predict(m_MT0, type = "response", newdata = pred_data1, se.fit = TRUE)


# Plot des points observés
plot(BDD_esp_net$LDMC, BDD_esp_net$MT, 
     xlab = "LDMC", ylab = "Maximum Temperature", 
     main = "MT selon LDMC")

# Courbes de prédiction
lines(foo, pred1$fit, col = "blue", lwd = 2) 

# Intervalle de confiance pour SD moyen
lines(foo, pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines(foo, pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)


#################pour les autres variables ###################
# Modèle GLM 
m_MT0<-glm(MT~SD+TD+Surface_F+LDMC+LT,data=BDD_esp_net,family="gaussian")
summary(m_MT0)

# Valeurs de LDMC
foo <- seq(min(BDD_esp_net$Surface_F,na.rm = TRUE), max(BDD_esp_net$Surface_F,na.rm = TRUE))

# on fait varier 1 variable et on fixe les autres variables
pred_data1 <- data.frame(
  Surface_F = foo,
  LDMC = rep(mean(BDD_esp_net$LDMC,na.rm = TRUE), length(foo)),
  TD = rep(mean(BDD_esp_net$TD,na.rm = TRUE), length(foo)),
  LT = rep(mean(BDD_esp_net$LT,na.rm = TRUE), length(foo)),
  SD = rep(mean(BDD_esp_net$SD,na.rm = TRUE), length(foo))
)


# Prédictions
pred1 <- predict(m_MT0, type = "response", newdata = pred_data1, se.fit = TRUE)


# Plot des points observés
plot(BDD_esp_net$Surface_F, BDD_esp_net$MT, 
     xlab = "Surface_F", ylab = "Maximum Temperature", 
     main = "MT selon Surface_F")

# Courbes de prédiction
lines(foo, pred1$fit, col = "#329c2f", lwd = 2) 

# Intervalle de confiance pour Surface_F moyen
lines(foo, pred1$fit + 1.96 * pred1$se.fit, col = "#329c2f", lty = 3)
lines(foo, pred1$fit - 1.96 * pred1$se.fit, col = "#329c2f", lty = 3)





##################### BB #######################
############ BB ######################
# Modèle GLM 
m_BB0<-glm(BB_prop~SD+TD+Surface_F+LDMC+LT,data=BDD_esp_net,family="gaussian")
summary(m_BB0)

# Valeurs de LDMC
foo <- seq(min(BDD_esp_net$LDMC, na.rm = TRUE), max(BDD_esp_net$LDMC,na.rm = TRUE))

# on fait varier 1 variable et on fixe les autres variables
pred_data1 <- data.frame(
  LDMC = foo,
  SD = rep(mean(BDD_esp_net$SD,na.rm = TRUE), length(foo)),
  TD = rep(mean(BDD_esp_net$TD), length(foo)),
  Surface_F = rep(mean(BDD_esp_net$Surface_F), length(foo)),
  LT = rep(mean(BDD_esp_net$LT), length(foo))
)


# Prédictions
pred1 <- predict(m_BB0, type = "response", newdata = pred_data1, se.fit = TRUE)


# Plot des points observés
plot(BDD_esp_net$LDMC, BDD_esp_net$BB_prop, 
     xlab = "LDMC", ylab = "Burnt Biomass", 
     main = "BB selon LDMC")

# Courbes de prédiction
lines(foo, pred1$fit, col = "blue", lwd = 2) 

# Intervalle de confiance pour SD moyen
lines(foo, pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines(foo, pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)


#################pour les autres variables ###################
# Modèle GLM 
m_BB0<-glm(BB_prop~SD+TD+Surface_F+LDMC+LT,data=BDD_esp_net,family="gaussian")
summary(m_BB0)

# Valeurs de LDMC
foo <- seq(min(BDD_esp_net$LT,na.rm = TRUE), max(BDD_esp_net$LT,na.rm = TRUE))

# on fait varier 1 variable et on fixe les autres variables
pred_data1 <- data.frame(
  LT = foo,
  LDMC = rep(mean(BDD_esp_net$LDMC,na.rm = TRUE), length(foo)),
  TD = rep(mean(BDD_esp_net$TD,na.rm = TRUE), length(foo)),
  Surface_F = rep(mean(BDD_esp_net$Surface_F,na.rm = TRUE), length(foo)),
  SD = rep(mean(BDD_esp_net$SD,na.rm = TRUE), length(foo))
)


# Prédictions
pred1 <- predict(m_BB0, type = "response", newdata = pred_data1, se.fit = TRUE)


# Plot des points observés
plot(BDD_esp_net$LT, BDD_esp_net$BB_prop, 
     xlab = "LT", ylab = "Burnt Biomass", 
     main = "BB selon LT")

# Courbes de prédiction
lines(foo, pred1$fit, col = "#329c2f", lwd = 2) 

# Intervalle de confiance pour LT moyen
lines(foo, pred1$fit + 1.96 * pred1$se.fit, col = "#329c2f", lty = 3)
lines(foo, pred1$fit - 1.96 * pred1$se.fit, col = "#329c2f", lty = 3)









##################### BT #######################
############ BT ######################
# Modèle GLM 
m_BT0<-glm(BT~SD+TD+Surface_F+LDMC+LT,data=BDD_esp_net,family="Gamma")
summary(m_BT0)

# Valeurs de LDMC
foo <- seq(min(BDD_esp_net$LDMC, na.rm = TRUE), max(BDD_esp_net$LDMC,na.rm = TRUE))

# on fait varier 1 variable et on fixe les autres variables
pred_data1 <- data.frame(
  LDMC = foo,
  SD = rep(mean(BDD_esp_net$SD,na.rm = TRUE), length(foo)),
  TD = rep(mean(BDD_esp_net$TD), length(foo)),
  Surface_F = rep(mean(BDD_esp_net$Surface_F), length(foo)),
  LT = rep(mean(BDD_esp_net$LT), length(foo))
)


# Prédictions
pred1 <- predict(m_BT0, type = "response", newdata = pred_data1, se.fit = TRUE)


# Plot des points observés
plot(BDD_esp_net$LDMC, BDD_esp_net$BT, 
     xlab = "LDMC", ylab = "Burning Time", 
     main = "BT selon LDMC")

# Courbes de prédiction
lines(foo, pred1$fit, col = "blue", lwd = 2) 

# Intervalle de confiance pour SD moyen
lines(foo, pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines(foo, pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)


#################pour les autres variables ###################
# Modèle GLM 
m_BT0<-glm(BT~SD+TD+Surface_F+LDMC+LT,data=BDD_esp_net,family="Gamma")
summary(m_BT0)

# Valeurs de LDMC
foo <- seq(min(BDD_esp_net$SD,na.rm = TRUE), max(BDD_esp_net$SD,na.rm = TRUE),length.out = 100)

# on fait varier 1 variable et on fixe les autres variables
pred_data1 <- data.frame(
  SD = foo,
  LDMC = rep(mean(BDD_esp_net$LDMC,na.rm = TRUE), length(foo)),
  TD = rep(mean(BDD_esp_net$TD,na.rm = TRUE), length(foo)),
  Surface_F = rep(mean(BDD_esp_net$Surface_F,na.rm = TRUE), length(foo)),
  LT = rep(mean(BDD_esp_net$LT,na.rm = TRUE), length(foo))
)


# Prédictions
pred1 <- predict(m_BT0, type = "response", newdata = pred_data1, se.fit = TRUE)


# Plot des points observés
plot(BDD_esp_net$SD, BDD_esp_net$BT, 
     xlab = "SD", ylab = "Burning Time", 
     main = "BT selon SD")

# Courbes de prédiction
lines(foo, pred1$fit, col = "#329c2f", lwd = 2) 

# Intervalle de confiance pour SD moyen
lines(foo, pred1$fit + 1.96 * pred1$se.fit, col = "#329c2f", lty = 3)
lines(foo, pred1$fit - 1.96 * pred1$se.fit, col = "#329c2f", lty = 3)





##################### DI #######################
############ DI ######################
# Modèle GLM 
m_DI0<-glm(DI~SD+TD+Surface_F+LMC_t24+LT,data=BDD_esp_net,family="Gamma")
summary(m_DI0)

# Valeurs de LMC_t24
foo <- seq(min(BDD_esp_net$LMC_t24, na.rm = TRUE), max(BDD_esp_net$LMC_t24,na.rm = TRUE))

# on fait varier 1 variable et on fixe les autres variables
pred_data1 <- data.frame(
  LMC_t24 = foo,
  SD = rep(mean(BDD_esp_net$SD,na.rm = TRUE), length(foo)),
  TD = rep(mean(BDD_esp_net$TD), length(foo)),
  Surface_F = rep(mean(BDD_esp_net$Surface_F), length(foo)),
  LT = rep(mean(BDD_esp_net$LT), length(foo))
)


# Prédictions
pred1 <- predict(m_DI0, type = "response", newdata = pred_data1, se.fit = TRUE)


# Plot des points observés
plot(BDD_esp_net$LMC_t24, BDD_esp_net$DI, 
     xlab = "LMC_t24", ylab = "Ignition Delay", 
     main = "DI selon LMC_t24")

# Courbes de prédiction
lines(foo, pred1$fit, col = "blue", lwd = 2) 

# Intervalle de confiance pour SD moyen
lines(foo, pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines(foo, pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)


#################pour les autres variables ###################
# Modèle GLM 
m_DI0<-glm(DI~SD+TD+Surface_F+LDMC+LT,data=BDD_esp_net,family="Gamma")
summary(m_DI0)

# Valeurs de LDMC
foo <- seq(min(BDD_esp_net$TD,na.rm = TRUE), max(BDD_esp_net$TD,na.rm = TRUE),length.out = 100)

# on fait varier 1 variable et on fixe les autres variables
pred_data1 <- data.frame(
  TD = foo,
  LDMC = rep(mean(BDD_esp_net$LDMC,na.rm = TRUE), length(foo)),
  SD = rep(mean(BDD_esp_net$SD,na.rm = TRUE), length(foo)),
  Surface_F = rep(mean(BDD_esp_net$Surface_F,na.rm = TRUE), length(foo)),
  LT = rep(mean(BDD_esp_net$LT,na.rm = TRUE), length(foo))
)


# Prédictions
pred1 <- predict(m_DI0, type = "response", newdata = pred_data1, se.fit = TRUE)


# Plot des points observés
plot(BDD_esp_net$TD, BDD_esp_net$DI, 
     xlab = "TD", ylab = "Ignition Delay", 
     main = "DI selon TD")

# Courbes de prédiction
lines(foo, pred1$fit, col = "#329c2f", lwd = 2) 

# Intervalle de confiance pour TD moyen
lines(foo, pred1$fit + 1.96 * pred1$se.fit, col = "#329c2f", lty = 3)
lines(foo, pred1$fit - 1.96 * pred1$se.fit, col = "#329c2f", lty = 3)














############################### variabilité ################################
kruskal.test(MT ~ Nom_scientifique, data = BDD_esp_netMT)

