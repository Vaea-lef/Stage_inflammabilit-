#### ANALYSE DE DONNEE STAGE INFLAMMABILITE ####


#packages
library(FactoMineR)
library(factoextra)
library (corrplot)




################# IMPORTATION ET CHARGEMENT PACKAGE ###########################

setwd("C:/IRD/Stage_inflammabilit-") #définition du répertoire de travail

##  importation BDD_ana_ech en format CSV
BDD_ech<-read.csv2("Data/BDD_moy_ech.csv", header = TRUE) #importation de l
names(BDD_ech)[which(names(BDD_ech) == "Nb_ramifications")] <- "Nb_rami"

BDD_echMT<-read.csv2("Data/BDD_moy_ech1.csv", header = TRUE) #importation de l
names(BDD_echMT)[which(names(BDD_echMT) == "Nb_ramifications")] <- "Nb_rami"

BDD_ech3<-read.csv2("Data/BDD_moy_ech2.csv", header = TRUE) #importation de l
names(BDD_ech3)[which(names(BDD_ech3) == "Nb_ramifications")] <- "Nb_rami"

##  importation BDD_sd_esp en format CSV
BDD_sd_esp<-read.csv2("Data/BDD_sd_esp.csv", header = TRUE)

##  importation BDD_moy_esp en format CSV
BDD_esp<-read.csv2("Data/BDD_moy_esp.csv", header = TRUE) #importation de la base
names(BDD_esp)[which(names(BDD_esp) == "Nb_ramifications")] <- "Nb_rami"

##  importation BDD_moy_esp spéciale MT en format CSV
BDD_esp_netMT<-read.csv2("Data/BDD_moy_esp1.csv", header = TRUE) #importation de la base
names(BDD_esp_netMT)[which(names(BDD_esp_netMT) == "Nb_ramifications")] <- "Nb_rami"

##  importation BDD_moy_esp spéciale BT, BB, DI en format CSV
BDD_esp_net3<-read.csv2("Data/BDD_moy_esp2.csv", header = TRUE) #importation de la base
names(BDD_esp_net3)[which(names(BDD_esp_net3) == "Nb_ramifications")] <- "Nb_rami"




################## DISTRIBUTION DES DONNEES #####################
#Infla
hist(BDD_esp_net3$DI_test,xlab="DI",main="DI distribution",xlim=c(0,6),breaks=seq(0,6,0.5)) # asymétrique droite
hist(BDD_esp_net3$BT_test,xlab="BT",main="BT distribution",ylim=c(0,25)) # asymétrique droite
hist(BDD_esp_net3$BB_test,,xlab="BB",main="BB distribution") # normale proportion
hist(BDD_esp_netMT$MT,xlab="MT",main="MT distribution",xlim=c(400,1000),ylim=c(0,15),breaks=seq(400,1000,50))      # normale 


#traits
hist(BDD_esp$Nb_rami)  
hist(BDD_esp$SD)       
hist(BDD_esp$TMC_t0)         
hist(BDD_esp$TMC_t24)  
hist(BDD_esp$TDMC)         
hist(BDD_esp$TD)            
hist(BDD_esp$TDIA)           
hist(BDD_esp$LMC_t0)         
hist(BDD_esp$LMC_t24)  
hist(BDD_esp$LDMC)          
hist(BDD_esp$Surface_F) 
hist(BDD_esp$SLA)            
hist(BDD_esp$LT)        





############################### ACP INFLA ######################################

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

HCPC(res.pca,method="ward")

?HCPC




####################### CLASSEMENT ESPECES ############################
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
round(mat_cor_trait, 2)   ## afficher les valeurs
corrplot(mat_cor_trait, method = "color", tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")

library(caret)
traits_non_corrélés <- colonnes_traits[, -findCorrelation(mat_cor_trait, cutoff = 0.68)]
traits_non_corrélés






############# PLOT INFLA ESPECE ######################

library(ggplot2)
############### MT ####################
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

############### DI ####################
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



############### BT ####################
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


############### BT ####################
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


############### FI ####################
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












############################ MODELES #####################################

#standardisation des données 
BDD_esp$SD_cr<-as.numeric(scale(BDD_esp$SD))
BDD_esp$TD_cr<-as.numeric(scale(BDD_esp$TD))
BDD_esp$LA_cr<-as.numeric(scale(BDD_esp$Surface_F))
BDD_esp$LDMC_cr<-as.numeric(scale(BDD_esp$LDMC))
BDD_esp$LT_cr<-as.numeric(scale(BDD_esp$LT))
BDD_esp$Vent_cr<-as.numeric(scale(BDD_esp$Vent))
BDD_esp$Temp_cr<-as.numeric(scale(BDD_esp$T_ambiante))
BDD_esp$LMC_t24_cr<-as.numeric(scale(BDD_esp$LMC_t24))
BDD_esp$Nb_rami_cr<-as.numeric(scale(BDD_esp$Nb_rami))

BDD_esp_netMT$SD_cr<-as.numeric(scale(BDD_esp_netMT$SD))
BDD_esp_netMT$TD_cr<-as.numeric(scale(BDD_esp_netMT$TD))
BDD_esp_netMT$LA_cr<-as.numeric(scale(BDD_esp_netMT$Surface_F))
BDD_esp_netMT$LDMC_cr<-as.numeric(scale(BDD_esp_netMT$LDMC))
BDD_esp_netMT$LT_cr<-as.numeric(scale(BDD_esp_netMT$LT))
BDD_esp_netMT$Vent_cr<-as.numeric(scale(BDD_esp_netMT$Vent))
BDD_esp_netMT$Temp_cr<-as.numeric(scale(BDD_esp_netMT$T_ambiante))
BDD_esp_netMT$LMC_t24_cr<-as.numeric(scale(BDD_esp_netMT$LMC_t24))
BDD_esp_netMT$Nb_rami_cr<-as.numeric(scale(BDD_esp_netMT$Nb_rami))

BDD_esp_net3$SD_cr<-as.numeric(scale(BDD_esp_net3$SD))
BDD_esp_net3$TD_cr<-as.numeric(scale(BDD_esp_net3$TD))
BDD_esp_net3$LA_cr<-as.numeric(scale(BDD_esp_net3$Surface_F))
BDD_esp_net3$LDMC_cr<-as.numeric(scale(BDD_esp_net3$LDMC))
BDD_esp_net3$LT_cr<-as.numeric(scale(BDD_esp_net3$LT))
BDD_esp_net3$Vent_cr<-as.numeric(scale(BDD_esp_net3$Vent))
BDD_esp_net3$Temp_cr<-as.numeric(scale(BDD_esp_net3$T_ambiante))
BDD_esp_net3$LMC_t24_cr<-as.numeric(scale(BDD_esp_net3$LMC_t24))
BDD_esp_net3$Nb_rami_cr<-as.numeric(scale(BDD_esp_net3$Nb_rami))





###### FI ########
# Comptage du nombre d'essais par espèce dans BDD_ech
essais_par_espece <- table(BDD_ech$Nom_scientifique)
essais_par_espece

# Création d'une nouvelle colonne Nb_essais dans BDD_esp en s'appuyant sur le nom scientifique
BDD_esp$Nb_essais <- essais_par_espece[BDD_esp$Nom_scientifique]

# Vérification
head(BDD_esp)

#modèle
mFI3<-glm(cbind(Nb_FI,Nb_essais-Nb_FI)~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp,family="binomial")
summary(mFI3)

mFI4<-glm(cbind(Nb_FI,Nb_essais-Nb_FI)~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LT_cr,data=BDD_esp,family="binomial")
summary(mFI4)

AIC(mFI3,mFI4)

###################### graph vert rouge ###############

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

###################### prédicion ############### mettre LMC_t24 ou LDMC
# Modèle GLM binomial (mFI3 ou mFI4)
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
pred1 <- predict(mFI4, type = "response", newdata = pred_data1, se.fit = TRUE)


# Plot des points observés
plot(BDD_esp$LMC_t24, BDD_esp$Nb_FI / BDD_esp$Nb_essais, 
     xlab = "LMC_t24", ylab = "Ignition Frequency",ylim = c(0,1), 
     main = "Effet de LMC_t24 sur LMC_t24 selon SD")

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "blue", lwd = 2) 

axis(2)

# Intervalle de confiance pour SD moyen
# Calcul borne supérieure (limitée à 1)
lines((foo * ecart + moy), pmin(pred1$fit + 1.96 * pred1$se.fit, 1), col = "blue", lty = 3)

# Calcul borne inférieure (limitée à 0)
lines((foo * ecart + moy), pmax(pred1$fit - 1.96 * pred1$se.fit, 0), col = "blue", lty = 3)






###### MT #########
#modèles
m_MT0<-glm(MT~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp_netMT,family="gaussian")
summary(m_MT0)
m_MT1<-glm(MT~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LT_cr,data=BDD_esp_netMT,family="gaussian")
summary(m_MT1)


AIC(m_MT0,m_MT1)

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


###################### prédicion ############### mettre LMC_t24 ou LDMC
# Modèle GLM binomial (m_MT0 ou m_MT1)
moy <- mean(BDD_esp_netMT$LDMC, na.rm = TRUE)
ecart <- sd(BDD_esp_netMT$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_esp_netMT$LDMC_cr), max(BDD_esp_netMT$LDMC_cr),length.out = 100)

# on fait varier LMDC et on fixe les autres variables
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp_netMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_esp_netMT$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp_netMT$LA_cr), length(foo)),
  LT_cr = rep(mean(BDD_esp_netMT$LT_cr), length(foo))
)

pred_data2 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp_netMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_esp_netMT$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp_netMT$LA_cr), length(foo)),
  LT_cr = rep(quantile(BDD_esp_netMT$LT_cr,0.75), length(foo))
)

pred_data3 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp_netMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_esp_netMT$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp_netMT$LA_cr), length(foo)),
  LT_cr = rep(quantile(BDD_esp_netMT$LT_cr,0.25), length(foo))
)

# Prédictions
pred1 <- predict(m_MT0, type = "response", newdata = pred_data1, se.fit = TRUE)
pred2 <- predict(m_MT0, type = "response", newdata = pred_data2, se.fit = TRUE)
pred3 <- predict(m_MT0, type = "response", newdata = pred_data3, se.fit = TRUE)

# Plot des points observés
plot(BDD_esp_netMT$LDMC, BDD_esp_netMT$MT, 
     xlab = "LDMC", ylab = "Maximum temperature", 
     main = "Effet de LDMC sur LDMC selon SD")

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 
lines((foo*ecart+moy), pred2$fit, col = "red" ,lty=3) 
lines((foo*ecart+moy), pred3$fit, col = "red",lty=3) 

axis(2)

# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)








###### BB #########
#modèles
BDD_esp_net3$BB_prop <- BDD_esp_net3$BB_test/100


m_BB0<-glm(BB_prop~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp_net3,family="gaussian")
summary(m_BB0)
m_BB1<-glm(BB_prop~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LT_cr,data=BDD_esp_net3,family="gaussian")
summary(m_BB1)



###################### graph coefs ###############

# Extraire coefficients (sans l'intercept)
par(mar = c(5,7,5,5))
coefs <- summary(m_BB1)$coefficients[-1, ]

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
mtext(expression(R^2~"= 0.43"), side = 3, adj = 0, line = 0.5, cex = 1.5)

################# Prédiction ########################
###################### prédicion ############### mettre LMC_t24 ou LDMC
# Modèle GLM binomial (m_BB0 ou m_BB1)
moy <- mean(BDD_esp_net3$LDMC, na.rm = TRUE)
ecart <- sd(BDD_esp_net3$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_esp_net3$LDMC_cr), max(BDD_esp_net3$LDMC_cr),length.out = 100)

# on fait varier LMDC et on fixe les autres variables
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp_net3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_esp_net3$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp_net3$LA_cr), length(foo)),
  LT_cr = rep(mean(BDD_esp_net3$LT_cr), length(foo))
)

pred_data2 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp_net3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_esp_net3$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp_net3$LA_cr), length(foo)),
  LT_cr = rep(quantile(BDD_esp_net3$LT_cr,0.75), length(foo))
)

pred_data3 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp_net3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_esp_net3$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp_net3$LA_cr), length(foo)),
  LT_cr = rep(quantile(BDD_esp_net3$LT_cr,0.25), length(foo))
)

# Prédictions
pred1 <- predict(m_BB0, type = "response", newdata = pred_data1, se.fit = TRUE)
pred2 <- predict(m_BB0, type = "response", newdata = pred_data2, se.fit = TRUE)
pred3 <- predict(m_BB0, type = "response", newdata = pred_data3, se.fit = TRUE)

# Plot des points observés
plot(BDD_esp_net3$LDMC, BDD_esp_net3$BB_prop, 
     xlab = "LDMC", ylab = "Burnt Biomass", 
     main = "Effet de LDMC sur BB selon SD")

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 
lines((foo*ecart+moy), pred2$fit, col = "red" ,lty=3) 
lines((foo*ecart+moy), pred3$fit, col = "red",lty=3) 

axis(2)

# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)




###### BT #########
#modèles
m_BT0<-glm(BT_test~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp_net3,family=Gamma(link="log"))
summary(m_BT0)
m_BT1<-glm(BT_test~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LT_cr,data=BDD_esp_net3,family=Gamma(link="log"))
summary(m_BT1)

AIC(m_BT0,m_BT1)

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
     xlim = c(-0.5,0.5),
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
axis(2, at = y_pos, labels = c("SD", "TD", "LA", "LDMC","LT"), las = 1,cex.axis = 1.5)

# Axe X
axis(1,cex.axis = 1.5)

# Valeurs des coefficients + étoiles
text(estimates, y_pos + 0.15,
     labels = paste0(round(estimates, 4), stars),
     col = cols, font = 2, cex =  1.7)

# Ajouter un texte avec le pseudo R² (à modifier si besoin)
mtext(expression(R^2~"= 0.62"), side = 3, adj = 0, line = 0.5, cex = 1.5)


############### prédiction ###############
###################### prédicion ############### mettre LMC_t24 ou LDMC
# Modèle GLM binomial (m_BT0 ou m_BT1)
moy <- mean(BDD_esp_net3$LDMC, na.rm = TRUE)
ecart <- sd(BDD_esp_net3$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_esp_net3$LDMC_cr), max(BDD_esp_net3$LDMC_cr),length.out = 100)

# on fait varier LMDC et on fixe les autres variables
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp_net3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_esp_net3$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp_net3$LA_cr), length(foo)),
  LT_cr = rep(mean(BDD_esp_net3$LT_cr), length(foo))
)

pred_data2 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp_net3$SD_cr,na.rm = TRUE), length(foo)),
  LT_cr = rep(mean(BDD_esp_net3$LT_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp_net3$LA_cr), length(foo)),
  TD_cr = rep(quantile(BDD_esp_net3$TD_cr,0.75), length(foo))
)

pred_data3 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp_net3$SD_cr,na.rm = TRUE), length(foo)),
  LT_cr = rep(mean(BDD_esp_net3$LT_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp_net3$LA_cr), length(foo)),
  TD_cr = rep(quantile(BDD_esp_net3$TD_cr,0.25), length(foo))
)

# Prédictions
pred1 <- predict(m_BT0, type = "response", newdata = pred_data1, se.fit = TRUE)
pred2 <- predict(m_BT0, type = "response", newdata = pred_data2, se.fit = TRUE)
pred3 <- predict(m_BT0, type = "response", newdata = pred_data3, se.fit = TRUE)

# Plot des points observés
plot(BDD_esp_net3$LDMC, BDD_esp_net3$BT, 
     xlab = "LDMC", ylab = "Burning Time", 
     main = "Effet de LDMC sur BT selon SD")

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 
lines((foo*ecart+moy), pred2$fit, col = "red" ,lty=3) 
lines((foo*ecart+moy), pred3$fit, col = "red",lty=3) 

axis(2)

# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)


###################
###### DI ######################################################################
###################

#modèles
m_DI0<-glm(DI_test~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp_net3,family=Gamma(link="log"))
summary(m_DI0)
m_DI1<-glm(DI_test~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LT_cr,data=BDD_esp_net3,family=Gamma(link="log"))
summary(m_DI1)

AIC(m_DI0,m_DI1)

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
axis(2, at = y_pos, labels = c("SD", "TD", "LA", "LDMC","LT"), las = 1,cex.axis = 1.5)

# Axe X
axis(1,cex.axis = 1.5)

# Valeurs des coefficients + étoiles
text(estimates, y_pos + 0.15,
     labels = paste0(round(estimates, 4), stars),
     col = cols, font = 2, cex =  1.7)

# Ajouter un texte avec le pseudo R² (à modifier si besoin)
mtext(expression(R^2~"= 0.52"), side = 3, adj = 0, line = 0.5, cex = 1.5)

############ prédiction ##########
###################### prédicion ############### mettre LMC_t24 ou LDMC
# Modèle GLM binomial (m_DI0 ou m_DI1)
moy <- mean(BDD_esp_net3$LDMC, na.rm = TRUE)
ecart <- sd(BDD_esp_net3$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_esp_net3$LDMC_cr), max(BDD_esp_net3$LDMC_cr),length.out = 100)

# on fait varier LMDC et on fixe les autres variables
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp_net3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_esp_net3$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp_net3$LA_cr), length(foo)),
  LT_cr = rep(mean(BDD_esp_net3$LT_cr), length(foo))
)


pred_data2 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp_net3$SD_cr,na.rm = TRUE), length(foo)),
  LT_cr = rep(mean(BDD_esp_net3$LT_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp_net3$LA_cr), length(foo)),
  TD_cr = rep(quantile(BDD_esp_net3$TD_cr,0.75), length(foo))
)

pred_data3 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp_net3$SD_cr,na.rm = TRUE), length(foo)),
  LT_cr = rep(mean(BDD_esp_net3$LT_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp_net3$LA_cr), length(foo)),
  TD_cr = rep(quantile(BDD_esp_net3$TD_cr,0.25), length(foo))
)

# Prédictions
pred1 <- predict(m_DI0, type = "response", newdata = pred_data1, se.fit = TRUE)
pred2 <- predict(m_DI0, type = "response", newdata = pred_data2, se.fit = TRUE)
pred3 <- predict(m_DI0, type = "response", newdata = pred_data3, se.fit = TRUE)

# Plot des points observés
plot(BDD_esp_net3$LDMC, BDD_esp_net3$DI_test, 
     xlab = "LDMC", ylab = "Ignition Delay", 
     main = "Effet de LDMC sur DI selon SD")

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 
lines((foo*ecart+moy), pred2$fit, col = "red" ,lty=3) 
lines((foo*ecart+moy), pred3$fit, col = "red",lty=3) 

axis(2)

# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)




###### Score ########
#modèles
m_score0<-glm(score_normalise~SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr,data=BDD_esp,family="gaussian")
summary(m_score0)
m_score1<-glm(score_normalise~SD_cr+TD_cr+LA_cr+LMC_t24_cr+LT_cr,data=BDD_esp,family="gaussian")
summary(m_score1)


AIC(m_score0,m_score1)


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

###################### prédicion ############### mettre LMC_t24 ou LDMC
# Modèle GLM binomial (m_score0 ou m_score1)
moy <- mean(BDD_esp$LDMC, na.rm = TRUE)
ecart <- sd(BDD_esp$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_esp$LDMC_cr), max(BDD_esp$LDMC_cr),length.out = 100)

# on fait varier LMDC et on fixe les autres variables
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_esp$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp$LA_cr), length(foo)),
  LT_cr = rep(mean(BDD_esp$LT_cr), length(foo))
)

pred_data2 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp$SD_cr,na.rm = TRUE), length(foo)),
  LT_cr = rep(mean(BDD_esp$LT_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp$LA_cr), length(foo)),
  TD_cr = rep(quantile(BDD_esp$TD_cr,0.75), length(foo))
)

pred_data3 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_esp$SD_cr,na.rm = TRUE), length(foo)),
  LT_cr = rep(mean(BDD_esp$LT_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp$LA_cr), length(foo)),
  TD_cr = rep(quantile(BDD_esp$TD_cr,0.25), length(foo))
)

# Prédictions
pred1 <- predict(m_score0, type = "response", newdata = pred_data1, se.fit = TRUE)
pred2 <- predict(m_score0, type = "response", newdata = pred_data2, se.fit = TRUE)
pred3 <- predict(m_score0, type = "response", newdata = pred_data3, se.fit = TRUE)

# Plot des points observés
plot(BDD_esp$LDMC, BDD_esp$score, 
     xlab = "LDMC", ylab = "score", 
     main = "Effet de LDMC sur score",
     ylim=c(-3,3))

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 
lines((foo*ecart+moy), pred2$fit, col = "red" ,lty=3) 
lines((foo*ecart+moy), pred3$fit, col = "red",lty=3) 

axis(2)

# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)


###################### prédicion ############### mettre LMC_t24 ou LMC_t24
# Modèle GLM binomial (m_score0 ou m_score1)
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

pred_data2 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_esp$SD_cr,na.rm = TRUE), length(foo)),
  LT_cr = rep(mean(BDD_esp$LT_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp$LA_cr), length(foo)),
  TD_cr = rep(quantile(BDD_esp$TD_cr,0.75), length(foo))
)

pred_data3 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_esp$SD_cr,na.rm = TRUE), length(foo)),
  LT_cr = rep(mean(BDD_esp$LT_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp$LA_cr), length(foo)),
  TD_cr = rep(quantile(BDD_esp$TD_cr,0.25), length(foo))
)

# Prédictions
pred1 <- predict(m_score1, type = "response", newdata = pred_data1, se.fit = TRUE)
pred2 <- predict(m_score1, type = "response", newdata = pred_data2, se.fit = TRUE)
pred3 <- predict(m_score1, type = "response", newdata = pred_data3, se.fit = TRUE)

# Plot des points observés
plot(BDD_esp$LMC_t24, BDD_esp$score, 
     xlab = "LMC_t24", ylab = "Score", 
     main = "Effet de LMC_t24 sur score ",
     ylim=c(-3,3))

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 
lines((foo*ecart+moy), pred2$fit, col = "red" ,lty=3) 
lines((foo*ecart+moy), pred3$fit, col = "red",lty=3) 

axis(2)

# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)



