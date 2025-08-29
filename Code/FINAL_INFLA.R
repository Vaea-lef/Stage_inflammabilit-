#### ANALYSE DE DONNEE STAGE INFLAMMABILITE ####


#packages
library(FactoMineR)
library(factoextra)
library (corrplot)
library(lmerTest)
library(lme4)
library(MuMIn)


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






















################################################################################
############################ ECHELLE ESPECE ####################################
################################################################################



################## DISTRIBUTION DES DONNEES #####################
#Infla
hist(BDD_esp_net3$DI_test,xlab="DI",main="DI distribution",xlim=c(0,5),breaks=seq(0,5,0.5)) # asymétrique droite
hist(BDD_esp_net3$BT_test,xlab="BT",main="BT distribution",breaks=seq(0,120,10)) # asymétrique droite
hist(BDD_esp_net3$BB_test,,xlab="BB",main="BB distribution") # normale proportion
hist(BDD_echMT$MT,xlab="MT",main="MT distribution",breaks=seq(0,1000,100))      # normale 

#Infla ech
hist(BDD_ech3$DI_test,xlab="DI",main="DI distribution",xlim=c(0,6),breaks=seq(0,7,0.5)) # asymétrique droite
hist(BDD_ech3$BT_test,xlab="BT",main="BT distribution") # asymétrique droite
hist(BDD_ech3$BB_test,,xlab="BB",main="BB distribution",breaks=seq(0,100,10)) # normale proportion
hist(BDD_echMT$MT,xlab="MT",main="MT distribution") 


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





#################### ACP TRAITS (pour sélection)############################

# Sélection des colonnes des traits fonctionnels
colonnes_traits <- na.omit(BDD_esp[, setdiff(17:31, c(23,24, 27))])

#strandardiser les données
colonnes_traits_cr <- scale(log(colonnes_traits))

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















################################################################################
############################ ECHELLE ESPECE ###############################
################################################################################

#standardisation des données 
BDD_esp$SD_cr<-as.numeric(scale(BDD_esp$SD))
BDD_esp$TD_cr<-as.numeric(scale(BDD_esp$TD))
BDD_esp$LA_cr<-as.numeric(scale(BDD_esp$Surface_F))
BDD_esp$LDMC_cr<-as.numeric(scale(BDD_esp$LDMC))
BDD_esp$LMC_t24_cr<-as.numeric(scale(BDD_esp$LMC_t24))
BDD_esp$SLA_cr<-as.numeric(scale(BDD_esp$SLA))



############### MODELES #############################

###### FI ########
# Comptage du nombre d'essais par espèce dans BDD_ech
essais_par_espece <- table(BDD_ech$Nom_scientifique)
essais_par_espece

# Création d'une nouvelle colonne Nb_essais dans BDD_esp en s'appuyant sur le nom scientifique
BDD_esp$Nb_essais <- essais_par_espece[BDD_esp$Nom_scientifique]

# Vérification
head(BDD_esp)

#modèle
mFI3<-glm(cbind(Nb_FI,Nb_essais-Nb_FI)~SD_cr+TD_cr+SLA_cr+LA_cr+LDMC_cr,data=BDD_esp,family="binomial")
summary(mFI3)

mFI4<-glm(cbind(Nb_FI,Nb_essais-Nb_FI)~SD_cr+TD_cr+SLA_cr+LA_cr+LMC_t24_cr,data=BDD_esp,family="binomial")
summary(mFI4)

AIC(mFI3,mFI4)

###################### graph vert rouge ###############

# Extraire coefficients (sans l'intercept)
par(mar = c(5,10,5,5))
coefs1 <- summary(mFI3)$coefficients[-1, ]
coefs2 <- summary(mFI4)$coefficients[-1, ]

# Variables utiles
estimates1 <- coefs1[, "Estimate"]
stderr1 <- coefs1[, "Std. Error"]
pval1 <- coefs1[, "Pr(>|z|)"]
labels1 <- rownames(coefs1)

estimates2 <- coefs2[, "Estimate"]
stderr2 <- coefs2[, "Std. Error"]
pval2 <- coefs2[, "Pr(>|z|)"]
labels2 <- rownames(coefs2)

# Calcul des intervalles de confiance plus ou moins SE
ci1 <- cbind(estimates1 - stderr1, estimates1 + stderr1)
ci2 <- cbind(estimates2 - stderr2, estimates2 + stderr2)

# Étoiles de significativité
stars1 <- ifelse(pval1 < 0.001, "***",
                 ifelse(pval1 < 0.01, "**",
                        ifelse(pval1 < 0.05, "*",
                               ifelse(pval1 < 0.1, ".", ""))))

stars2 <- ifelse(pval2 < 0.001, "***",
                 ifelse(pval2 < 0.01, "**",
                        ifelse(pval2 < 0.05, "*",
                               ifelse(pval2 < 0.1, ".", ""))))


# Couleurs selon signe du coefficient
cols1 <- ifelse(estimates1 < 0, "red", "red")
cols2 <- ifelse(estimates2 < 0, "blue", "blue")

# Ordre des variables (du bas vers le haut)
y_pos1 <- length(estimates1):1
y_pos2 <- length(estimates2):1

# Plot de base
par(mar=c(4,6,2,2))
plot(estimates1, y_pos1,type = "n",
     xlim = c(-3,3),
     ylim = c(0.7,length(estimates1) +0.2 ),
     xlab = "Estimate",
     ylab = "",
     axes = FALSE,
     main = "Fréquence d'ignition",
     cex.main = 1.1)



# Ligne verticale à zéro
abline(v = 0, lty = 2)
abline(h = y_pos1 -0.2, lwd = 0.5, lty = 3, col="grey")

#points des coefs
points(estimates1, y_pos1 ,pch = 16, col = cols1,cex = 1.1)
points(estimates2, y_pos2 - 0.4,pch = 17, col = cols2,cex = 1.1)

# Barres d'erreur (IC)
segments(ci1[,1], y_pos1, ci1[,2], y_pos1, col = cols1)
segments(ci2[,1], y_pos2 - 0.4, ci2[,2], y_pos2 - 0.4, col = cols2)

# Axe Y avec noms des variables
axis(2, at = c(4.8,3.8,2.8,1.8,1,0.6), labels = c("SD","TD","SLA","LA","LDMC","LMC_t24"), las = 1,cex.axis = 0.9)


# Axe X
axis(1,cex.axis = 0.9)

# Valeurs des coefficients + étoiles
text(estimates1, y_pos1 + 0.16,
     labels = paste0(round(estimates1, 2), stars1),
     col = cols1, font = 2, cex =  0.9)
text(estimates2, y_pos2 -0.24,
     labels = paste0(round(estimates2, 2), stars2),
     col = cols2, font = 2, cex =  0.9)

###################### prédicion ############### mettre LMC_t24 ou LDMC
# Modèle GLM binomial (mFI3 ou mFI4)
moy <- mean(BDD_esp$LMC_t24, na.rm = TRUE)
moy
ecart <- sd(BDD_esp$LMC_t24, na.rm = TRUE)
ecart

# Valeurs de LMC_t24
foo <- seq(min(BDD_esp$LMC_t24_cr), max(BDD_esp$LMC_t24_cr),length.out = 100)

# on fait varier LMDC et on fixe les autres variables
pred_data1 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_esp$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_esp$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_esp$LA_cr), length(foo)),
 SLA_cr = rep(mean(BDD_esp$SLA_cr), length(foo))
)


# Prédictions et IC
pred1 <- predict(mFI4, type = "response", newdata = pred_data1, se.fit = TRUE)
pred_data1$LMC_t24_cr[pred1$fit<0.5]*ecart+moy
pred_data1$LMC_t24_cr[(pred1$fit+(1.96*pred1$se.fit))<0.5]*ecart+moy
pred_data1$LMC_t24_cr[(pred1$fit-(1.96*pred1$se.fit))<0.5]*ecart+moy

# Plot des points observés
plot(BDD_esp$LMC_t24, BDD_esp$Nb_FI / BDD_esp$Nb_essais, 
     xlab = "LMC_t24 (%)", ylab = "Fréquence d'ignition",ylim = c(0,1), 
     main = "Effet de LMC_t24 sur la fréquence d'ignition")

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 
abline(v = 329,  lty = 2, col = "blue")
abline(v = 382,  lty = 2,col="red")
abline(v = 462,  lty = 2, col = "blue")
segments(329, 0.5, 462,0.5, col = "blue",lwd = 2)
points(382,0.5,pch=19,col="red")

axis(2)

# Intervalle de confiance pour SD moyen
# Calcul borne supérieure (limitée à 1)
lines((foo * ecart + moy), pmin(pred1$fit + 1.96 * pred1$se.fit, 1), col = "blue", lty = 3)

# Calcul borne inférieure (limitée à 0)
lines((foo * ecart + moy), pmax(pred1$fit - 1.96 * pred1$se.fit, 0), col = "blue", lty = 3)



















################################################################################
############################ ECHELLE ECHANTILLON ###############################
################################################################################

############################### ACP INFLA (avec projection des traits) ######################################
head(BDD_ech)
BDD_ech$SD_mean <- BDD_ech$SD           # copier la colonne originale
BDD_ech$SD_mean[is.na(BDD_ech$SD_mean)] <- mean(BDD_ech$SD, na.rm = TRUE)
head(BDD_ech)
# Sélection des colonnes des composantes de l'inflammabilité
colonnes_infla <- BDD_ech[, c(10,12,13,15)]
head(colonnes_infla)

#sélection de variabes supplémentaires (traits)
colonnes_traits <- BDD_ech [, c(22,26,28,29,30,33)]
head(colonnes_traits)

# Centrage-réduction des données
colonnes_infla_cr <- scale(colonnes_infla)
colonnes_traits_cr <- scale(colonnes_traits)

# Nettoyage des données
colonnes_clean <- na.omit(colonnes_infla_cr)

# ACP centrée et réduite
res.pca <- prcomp(colonnes_clean, scale. = FALSE)
summary(res.pca)
res.pca

M_S1<-aggregate(res.pca$x[,1],by=list(BDD_ech$Nom_scientifique),mean)
SD_S1<-aggregate(res.pca$x[,1],by=list(BDD_ech$Nom_scientifique),sd)

M_S2<-aggregate(res.pca$x[,2],by=list(BDD_ech$Nom_scientifique),mean)
SD_S2<-aggregate(res.pca$x[,2],by=list(BDD_ech$Nom_scientifique),sd)

# Coordonnées
ind_coords <- res.pca$x                    # individus
var_coords <- res.pca$rotation             # variables
eig_vals <- res.pca$sdev^2                 # valeurs propres
explained_var <- round(100 * eig_vals / sum(eig_vals), 1)

colonnes_taits_coord <- cor(colonnes_traits_cr, ind_coords)

# Ajout d'un score d'inflammabilité basé sur la coordonnée de l'axe 1 
BDD_ech$score <- ind_coords[,1]


# Graphique de base
par(mar = c(5,5,5,5))  # marges

HC<-hclust(d=dist(cbind(M_S1$x,M_S2$x)),method="ward.D2")
plot(HC, hang = -1,labels=F, axes="n")
axis(2,cex.axis=0.6)
GR<-cutree(HC,k=5)
GR

COL<-character()
COL[GR==5]<-"#00610D"
COL[GR==4]<-"red"
COL[GR==3]<-"#63B802"
COL[GR==2]<-"#FFDD1F"
COL[GR==1]<-"#FF6D1F"

col=rgb()

# Tracer les individus
plot(ind_coords[,1], ind_coords[,2],
     xlim = range(ind_coords[,1]) * 1.2,
     ylim = c(-4,3.5),
     xlab = paste0("PC1 (", explained_var[1], "%)"),
     ylab = paste0("PC2 (", explained_var[2], "%)"),
     main = "ACP - Individus et variables d'inflammabilité",
     pch = 21, bg = rgb(0,0,0,0.5),col=NA,cex=0.8)

# Tracer les axes
abline(h = 0, v = 0, lty = 2)

# Ajouter les flèches des variables
arrows(0, 0,                        # départ
       var_coords[,1]*max(abs(ind_coords[,1])),
       var_coords[,2]*max(abs(ind_coords[,2])),
       length = 0.1, col = "red", lwd = 2)

# Variables supplémentaires en bleu
arrows(0, 0,
       colonnes_taits_coord[,1]*max(abs(ind_coords[,1])),
       colonnes_taits_coord[,2]*max(abs(ind_coords[,2])),
       length = 0.1, col = "#5490FF")

text(colonnes_taits_coord[,1]*max(abs(ind_coords[,1])),
     colonnes_taits_coord[,2]*max(abs(ind_coords[,2])),
     labels = colnames(colonnes_traits_cr),
     col = "blue", cex = 0.8)


# Tracer les individus
plot(ind_coords[,1], ind_coords[,2],
     xlim = range(ind_coords[,1]) * 1.2,
     ylim = c(-4,1.5),
     xlab = paste0("PC1 (", explained_var[1], "%)"),
     ylab = paste0("PC2 (", explained_var[2], "%)"),
     main = "ACP - Individus et variables d'inflammabilité",
     pch = 21, bg = rgb(190/255,190/255,190/255,0.5),col=NA,cex=0.8)

# Tracer les axes
abline(h = 0, v = 0, lty = 2)


segments(x0=M_S1$x,x1=M_S1$x,y0=M_S2$x-SD_S2$x,y1=M_S2$x+SD_S2$x)
segments(x0=M_S1$x-SD_S1$x,x1=M_S1$x+SD_S1$x,y0=M_S2$x,y1=M_S2$x)
points(M_S1$x,M_S2$x,pch=21,bg=COL,cex=1.5)




#ajout colonne groupe dans BDD
BDD_esp$groupe<-GR


# Reclassement de la variable groupe 
BDD_esp$groupe <- factor(BDD_esp$groupe, levels = c(5, 3, 2, 1, 4))
write.csv2(BDD_esp,"Data/BDD_esp_groupe.csv")

# boxplot MT
boxplot(MT ~ groupe, data = BDD_esp,
        main = "Distribution de MT par groupe",
        xlab = "Groupe",
        ylab = "MT (°C)",col=c("#00610D","#63B802","#FFDD1F","#FF6D1F", "red"))

# boxplot BT
boxplot(BT_test ~ groupe, data = BDD_esp,
        main = "Distribution de BT par groupe",
        xlab = "Groupe",
        ylab = "BT (s)",col=c("#00610D","#63B802","#FFDD1F","#FF6D1F", "red"))


# boxplot BB
boxplot(BB_test ~ groupe, data = BDD_esp,
        main = "Distribution de BB par groupe",
        xlab = "Groupe",
        ylab = "BB (%)",col=c("#00610D","#63B802","#FFDD1F","#FF6D1F", "red"))

# boxplot DI
boxplot(DI_test ~ groupe, data = BDD_esp,
        main = "Distribution de DI par groupe",
        xlab = "Groupe",
        ylab = "DI (s)",col=c("#00610D","#63B802","#FFDD1F","#FF6D1F", "red"))










#Calcul de la moyenne du score pour ajout dans BDD_esp
#création de table avec moyenne et sd pour chaque variable en fonction du nom de l'espèce
tem3<-BDD_ech[,5:34] ###sélection des colonnes comprenant les variables pour les intégrer dans la boucle
tem3
#création d'un bdd d'origine pour moyenne (sert pour merge)
BDD_moy_score <- aggregate(tem3[,1] ~ Nom_scientifique + ID_espece + Milieu_recolte, data = BDD_ech, FUN = mean, na.rm = TRUE)
BDD_moy_score[,4] <- round(BDD_moy_score[,4], 2)
colnames(BDD_moy_score)[4] <- colnames(tem3)[1]

#création d'un bdd d'origine pour sd (sert pour merge)
BDD_sd_score <- aggregate(tem3[,1] ~ Nom_scientifique + ID_espece + Milieu_recolte, data = BDD_ech, FUN = sd, na.rm = TRUE)
BDD_sd_score[,4] <- round(BDD_sd_score[,4], 2)
colnames(BDD_sd_score)[4] <- colnames(tem3)[1]

#Boucle pour les calcul des moyennes et écart-types
for (i in 2:ncol(tem3)) {
  
  # Moyenne
  tem3_moy_esp <- aggregate(tem3[, i] ~ Nom_scientifique + ID_espece+ Milieu_recolte, data = BDD_ech, FUN = mean, na.rm = TRUE)
  tem3_moy_esp[,4] <- round(tem3_moy_esp[,4], 2)
  colnames(tem3_moy_esp)[4] <- colnames(tem3)[i]
  BDD_moy_score <- merge(BDD_moy_score, tem3_moy_esp, by = c("Nom_scientifique", "ID_espece", "Milieu_recolte"), all = TRUE)
  
  # Ecart-type
  tem3_sd_esp <- aggregate(tem3[, i] ~ Nom_scientifique + ID_espece+ Milieu_recolte, data = BDD_ech, FUN = sd, na.rm = TRUE)
  tem3_sd_esp[,4] <- round(tem3_sd_esp[,4], 2)
  colnames(tem3_sd_esp)[4] <- colnames(tem3)[i]
  BDD_sd_score <- merge(BDD_sd_score, tem3_sd_esp, by = c("Nom_scientifique", "ID_espece" , "Milieu_recolte"), all = TRUE)
}

BDD_moy_score






####################### CLASSEMENT ESPECES ############################
par(mar = c(5,7,6,5))
# Palette de couleur
pal_vert_jaune <- colorRampPalette(c("darkgreen", "yellow"))
pal_jaune_rouge <- colorRampPalette(c("yellow", "red"))

# Définir les bornes du graph
mean_score <- mean(BDD_moy_score$score, na.rm = TRUE)
min_score <- min(BDD_moy_score$score, na.rm = TRUE)
max_score <- max(BDD_moy_score$score, na.rm = TRUE)

# Nombre de nuances de chaque côté
n_colors <- 100
n_left <- round((mean_score - min_score) / (max_score - min_score) * n_colors)
n_right <- n_colors - n_left

# Générer  couleurs
couleurs <- c(pal_vert_jaune(n_left), pal_jaune_rouge(n_right))

# Attribuer à chaque espèce une couleur selon sa position par rapport à la min–max
score_scaled <- round((BDD_moy_score$score - min_score) / (max_score - min_score) * (length(couleurs) - 1)) + 1

# Ordre des points
o <- order(BDD_moy_score$score)
couleur_points <- couleurs[score_scaled][o]

# Aligner les écarts-types sur l'ordre des espèces dans BDD_moy_score
# On suppose que les deux ont la même colonne Nom_scientifique
sd_aligned <- BDD_sd_score$score[match(BDD_moy_score$Nom_scientifique[o], BDD_sd_score$Nom_scientifique)]

# Y positions dans l'ordre trié
y_pos <- 1:length(BDD_moy_score$Nom_scientifique)

# Tracé
par(mar = c(4, 13, 0, 0))
plot(BDD_moy_score$score[o], 1:length(BDD_moy_score$Nom_scientifique), axes="n",xlim = c(-6,3))
# Axes
axis(2, at = 1:length(BDD_moy_score$Nom_scientifique), labels = BDD_moy_score$Nom_scientifique[o], las = 1, cex.axis = 0.85,font=3)
mtext("Espèces", side = 2, line = 11, cex = 1)
axis(1, at = seq(-6, 4, by = 1))
mtext("Score d'inflammabilité", side = 1, line = 3, cex = 1)

# Lignes de référence
abline(v = mean_score, lwd = 2)
abline(v = quantile(BDD_moy_score$score, na.rm = TRUE)[2], lwd = 2, lty = 3)
abline(v = quantile(BDD_moy_score$score, na.rm = TRUE)[4], lwd = 2, lty = 3)
abline(h = y_pos, lwd = 0.5, lty = 3, col="grey")


# Ajouter des segments horizontaux : de (score - sd) à (score + sd) pour chaque espèce
for (i in seq_along(y_pos)) {
  segments(
    x0 = BDD_moy_score$score[o][i] - sd_aligned[i],
    y0 = y_pos[i],
    x1 = BDD_moy_score$score[o][i] + sd_aligned[i],
    y1 = y_pos[i],
    col = "black", lwd = 2
  )
}


# Points colorés
points(BDD_moy_score$score[o], 1:length(BDD_moy_score$Nom_scientifique), 
       pch = 21, cex = 1.5, bg = couleur_points)





#################### ACP TRAITS ############################


# Sélection des colonnes des traits fonctionnels
colonnes_traits <- na.omit(BDD_ech[, setdiff(17:31, c(23,24, 27))])

# Vérifier les données
colonnes_traits


colonnes_traits_cr <- scale(log(colonnes_traits +0.01))
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
par(mar = c(10,10,10,10))
corrplot(mat_cor_trait, method = "color", tl.cex = 0.8, tl.col = "black", number.cex = 0.7, addCoef.col = "black")

library(caret)
traits_non_corrélés <- colonnes_traits[, -findCorrelation(mat_cor_trait, cutoff = 0.7)]
traits_non_corrélés




################## DISTRIBUTION DES DONNEES #####################

#Infla ech
hist(BDD_ech3$DI_test,xlab="DI",main="Distribution DI",xlim=c(0,7),breaks=seq(0,7,0.5)) # asymétrique droite
hist(BDD_ech3$BT_test,xlab="BT",main="Distribution BT") # asymétrique droite
hist(BDD_ech3$BB_test,,xlab="BB",main="Distribution BB",breaks=seq(0,100,10)) # normale proportion
hist(BDD_echMT$MT,xlab="MT",main="Distribution MT") 
hist(BDD_ech$score,xlab="Score d'inflammabilité",main="Distribution score")

#traits
hist(BDD_ech$Nb_rami)  
hist(BDD_ech$SD)       
hist(BDD_ech$TMC_t0)         
hist(BDD_ech$TMC_t24)  
hist(BDD_ech$TDMC)         
hist(BDD_ech$TD)            
hist(BDD_ech$TDIA)           
hist(BDD_ech$LMC_t0)         
hist(BDD_ech$LMC_t24)  
hist(BDD_ech$LDMC)          
hist(BDD_ech$Surface_F) 
hist(BDD_ech$SLA)            
hist(BDD_ech$LT)        










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

























############################ MODELES échelle ech ###############################

#standardisation des données 
BDD_ech$SD_cr<-as.numeric(scale(BDD_ech$SD))
BDD_ech$TD_cr<-as.numeric(scale(BDD_ech$TD))
BDD_ech$LA_cr<-as.numeric(scale(BDD_ech$Surface_F))
BDD_ech$LDMC_cr<-as.numeric(scale(BDD_ech$LDMC))
BDD_ech$LMC_t24_cr<-as.numeric(scale(BDD_ech$LMC_t24))
BDD_ech$VPD_cr<-as.numeric(scale(BDD_ech$VPD))
BDD_ech$SLA_cr<-as.numeric(scale(BDD_ech$SLA))

BDD_echMT$SD_cr<-as.numeric(scale(BDD_echMT$SD))
BDD_echMT$TD_cr<-as.numeric(scale(BDD_echMT$TD))
BDD_echMT$LA_cr<-as.numeric(scale(BDD_echMT$Surface_F))
BDD_echMT$LDMC_cr<-as.numeric(scale(BDD_echMT$LDMC))
BDD_echMT$LMC_t24_cr<-as.numeric(scale(BDD_echMT$LMC_t24))
BDD_echMT$VPD_cr<-as.numeric(scale(BDD_echMT$VPD))
BDD_echMT$SLA_cr<-as.numeric(scale(BDD_echMT$SLA))

BDD_ech3$SD_cr<-as.numeric(scale(BDD_ech3$SD))
BDD_ech3$TD_cr<-as.numeric(scale(BDD_ech3$TD))
BDD_ech3$LA_cr<-as.numeric(scale(BDD_ech3$Surface_F))
BDD_ech3$LDMC_cr<-as.numeric(scale(BDD_ech3$LDMC))
BDD_ech3$LMC_t24_cr<-as.numeric(scale(BDD_ech3$LMC_t24))
BDD_ech3$VPD_cr<-as.numeric(scale(BDD_ech3$VPD))
BDD_ech3$SLA_cr<-as.numeric(scale(BDD_ech3$SLA))


plot(BDD_ech$LDMC, BDD_ech$LMC_t24)
m<-glm(LMC_t24~LDMC,data=BDD_ech,family=Gamma(link="log"))
summary(m)


foo <- seq(min(BDD_ech$LDMC), max(BDD_ech$LDMC), length.out = 100)
pred_data1 <- data.frame(LDMC = foo)
pred1 <- predict(m, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
plot(BDD_ech$LDMC,BDD_ech$LMC_t24,main="Relation entre LDMC et LMC_t24",xlab="LDMC (mg/g)",ylab="LMC_t24 (%)")
lines(foo, pred1$fit, col = "black", lwd = 2) 
# Intervalle de confiance pour SD moyen
lines(foo, pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines(foo, pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)




#######################################################
###### MT #########
#######################################################

#modèles
m_MT1 <- lmer(MT ~ SD_cr+TD_cr+SLA_cr+LA_cr+LDMC_cr + (1 | Nom_scientifique), data = BDD_echMT)
summary(m_MT1)
r.squaredGLMM(m_MT1)
m_MT3 <- lmer(MT ~ SD_cr+TD_cr+SLA_cr+LA_cr+LMC_t24_cr + (1 | Nom_scientifique), data = BDD_echMT)
summary(m_MT3)
r.squaredGLMM(m_MT3)

AIC(m_MT1,m_MT3)

###################### graph coefs ###############

# Extraire coefficients (sans l'intercept)
par(mar = c(5,10,5,5))
coefs1 <- summary(m_MT1)$coefficients[-1, ]
coefs2 <- summary(m_MT3)$coefficients[-1, ]

# Variables utiles
estimates1 <- coefs1[, "Estimate"]
stderr1 <- coefs1[, "Std. Error"]
pval1 <- coefs1[, "Pr(>|t|)"]
labels1 <- rownames(coefs1)

estimates2 <- coefs2[, "Estimate"]
stderr2 <- coefs2[, "Std. Error"]
pval2 <- coefs2[, "Pr(>|t|)"]
labels2 <- rownames(coefs2)

# Calcul des intervalles de confiance plus ou moins SE
ci1 <- cbind(estimates1 - stderr1, estimates1 + stderr1)
ci2 <- cbind(estimates2 - stderr2, estimates2 + stderr2)

# Étoiles de significativité
stars1 <- ifelse(pval1 < 0.001, "***",
                ifelse(pval1 < 0.01, "**",
                       ifelse(pval1 < 0.05, "*",
                              ifelse(pval1 < 0.1, ".", ""))))

stars2 <- ifelse(pval2 < 0.001, "***",
                 ifelse(pval2 < 0.01, "**",
                        ifelse(pval2 < 0.05, "*",
                               ifelse(pval2 < 0.1, ".", ""))))


# Couleurs selon signe du coefficient
cols1 <- ifelse(estimates1 < 0, "red", "red")
cols2 <- ifelse(estimates2 < 0, "blue", "blue")

# Ordre des variables (du bas vers le haut)
y_pos1 <- length(estimates1):1
y_pos2 <- length(estimates2):1

# Plot de base
png("C:/IRD/Stage_inflammabilit-/Figures/FIN/Coef_MT.png",width = 500, height = 500,res=150)
par(mar=c(4,6,2,2))
plot(estimates1, y_pos1,type = "n",
     xlim = c(-100,100),
     ylim = c(0.7,length(estimates1) +0.2 ),
     xlab = "Estimate",
     ylab = "",
     axes = FALSE,
     main = "Température maximum",
     cex.main = 1.1)



# Ligne verticale à zéro
abline(v = 0, lty = 2)
abline(h = y_pos1 -0.2, lwd = 0.5, lty = 3, col="grey")


points(estimates1, y_pos1 ,pch = 16, col = cols1,cex = 1.1)
points(estimates2, y_pos2 - 0.4,pch = 17, col = cols2,cex = 1.1)

# Barres d'erreur (IC)
segments(ci1[,1], y_pos1, ci1[,2], y_pos1, col = cols1)
segments(ci2[,1], y_pos2 - 0.4, ci2[,2], y_pos2 - 0.4, col = cols2)

# Axe Y avec noms des variables
axis(2, at = c(4.8,3.8,2.8,1.8,1,0.6), labels = c("SD","TD","SLA","LA","LDMC","LMC_t24"), las = 1,cex.axis = 0.9)


# Axe X
axis(1,cex.axis = 0.9)

# Valeurs des coefficients + étoiles
text(estimates1, y_pos1 + 0.16,
     labels = paste0(round(estimates1, 2), stars1),
     col = cols1, font = 2, cex =  0.9)
text(estimates2, y_pos2 -0.24,
     labels = paste0(round(estimates2, 2), stars2),
     col = cols2, font = 2, cex =  0.9)


dev.off()

###################### prédicion ############### 

############  LDMC
moy <- mean(BDD_echMT$LDMC, na.rm = TRUE)
ecart <- sd(BDD_echMT$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_echMT$LDMC_cr), max(BDD_echMT$LDMC_cr),length.out = 100)

############## on fait varier LMDC
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_echMT$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_echMT$SLA_cr), length(foo))
)


# Prédictions
pred1 <- predict(m_MT1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_echMT$LDMC, BDD_echMT$MT,type="n", 
     xlab = "LDMC (mg/g)", ylab = "Température maximum (°C)", 
     main = "Effet de LDMC sur MT",ylim=c(500,850),xlim=c(175,550))

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 

# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)



############  LMC_t24
moy <- mean(BDD_echMT$LMC_t24, na.rm = TRUE)
ecart <- sd(BDD_echMT$LMC_t24, na.rm = TRUE)

# Valeurs de LMC_t24
foo <- seq(min(BDD_echMT$LMC_t24_cr), max(BDD_echMT$LMC_t24_cr),length.out = 100)

############## on fait varier LMDC
pred_data1 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_echMT$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_echMT$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_MT3, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_echMT$LMC_t24, BDD_echMT$MT,type="n", 
     xlab = "LMC_t24 (%)", ylab = "Température maximum (°C)", 
     main = "Effet de LMC_t24 sur MT",ylim=c(500,850),xlim=c(0,350))

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 

# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)


############## on fait varier LMDC et TD
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_echMT$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_echMT$SLA_cr), length(foo))
)

pred_data2 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(quantile(BDD_echMT$TD_cr,0.75), length(foo)),
  LA_cr = rep(mean(BDD_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_echMT$SLA_cr), length(foo))
)

pred_data3 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(quantile(BDD_echMT$TD_cr,0.25), length(foo)),
  LA_cr = rep(mean(BDD_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_echMT$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_MT1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred2 <- predict(m_MT1, type = "response", newdata = pred_data2, se.fit = TRUE, re.form = NA)
pred3 <- predict(m_MT1, type = "response", newdata = pred_data3, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_echMT$LDMC, BDD_echMT$MT,type="n", 
     xlab = "LDMC", ylab = "Maximum temperature", 
     main = "Effet de LDMC sur LDMC selon TD",ylim=c(500,850))

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 
lines((foo*ecart+moy), pred2$fit, col = "red" ,lty=3) 
lines((foo*ecart+moy), pred3$fit, col = "red",lty=3) 


# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)




##################################
###### BB #########
##################################

#modèles
BDD_ech3$BB_prop <- BDD_ech3$BB_test/100


m_BB0 <- lmer(BB_prop ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr+VPD_cr + (1 | Nom_scientifique), data = BDD_ech3)    #### meilleur modèle 
summary(m_BB0)
r.squaredGLMM(m_BB0)
m_BB1 <- lmer(BB_prop ~ SD_cr+TD_cr+SLA_cr+LA_cr+LDMC_cr + (1 | Nom_scientifique), data = BDD_ech3)    #### meilleur modèle 
summary(m_BB1)
r.squaredGLMM(m_BB1)
m_BB2 <- lmer(BB_prop ~ SD_cr+TD_cr+LA_cr+LMC_t24_cr+LT_cr+VPD_cr + (1 | Nom_scientifique), data = BDD_ech3)    #### meilleur modèle 
summary(m_BB2)
r.squaredGLMM(m_BB2)
m_BB3 <- lmer(BB_prop ~ SD_cr+TD_cr+SLA_cr+LA_cr+LMC_t24_cr + (1 | Nom_scientifique), data = BDD_ech3)    #### meilleur modèle 
summary(m_BB3)
r.squaredGLMM(m_BB3)

AIC(m_BB1,m_BB3)

###################### graph coefs ###############

# Extraire coefficients (sans l'intercept)
par(mar = c(5,10,5,5))
coefs1 <- summary(m_BB1)$coefficients[-1, ]
coefs2 <- summary(m_BB3)$coefficients[-1, ]

# Variables utiles
estimates1 <- coefs1[, "Estimate"]
stderr1 <- coefs1[, "Std. Error"]
pval1 <- coefs1[, "Pr(>|t|)"]
labels1 <- rownames(coefs1)

estimates2 <- coefs2[, "Estimate"]
stderr2 <- coefs2[, "Std. Error"]
pval2 <- coefs2[, "Pr(>|t|)"]
labels2 <- rownames(coefs2)

# Calcul des intervalles de confiance plus ou moins SE
ci1 <- cbind(estimates1 - stderr1, estimates1 + stderr1)
ci2 <- cbind(estimates2 - stderr2, estimates2 + stderr2)

# Étoiles de significativité
stars1 <- ifelse(pval1 < 0.001, "***",
                 ifelse(pval1 < 0.01, "**",
                        ifelse(pval1 < 0.05, "*",
                               ifelse(pval1 < 0.1, ".", ""))))

stars2 <- ifelse(pval2 < 0.001, "***",
                 ifelse(pval2 < 0.01, "**",
                        ifelse(pval2 < 0.05, "*",
                               ifelse(pval2 < 0.1, ".", ""))))


# Couleurs selon signe du coefficient
cols1 <- ifelse(estimates1 < 0, "red", "red")
cols2 <- ifelse(estimates2 < 0, "blue", "blue")

# Ordre des variables (du bas vers le haut)
y_pos1 <- length(estimates1):1
y_pos2 <- length(estimates2):1

# Plot de base
png("C:/IRD/Stage_inflammabilit-/Figures/FIN/Coef_BB.png",width = 500, height = 500,res=150)
par(mar=c(4,6,2,2))
plot(estimates1, y_pos1,type = "n",
     xlim = c(-0.3,0.3),
     ylim = c(0.7,length(estimates1) +0.2 ),
     xlab = "Estimate",
     ylab = "",
     axes = FALSE,
     main = "Biomasse brulée",
     cex.main = 1.1)



# Ligne verticale à zéro
abline(v = 0, lty = 2)
abline(h = y_pos1 -0.2, lwd = 0.5, lty = 3, col="grey")


points(estimates1, y_pos1 ,pch = 16, col = cols1,cex = 1.1)
points(estimates2, y_pos2 - 0.4,pch = 17, col = cols2,cex = 1.1)

# Barres d'erreur (IC)
segments(ci1[,1], y_pos1, ci1[,2], y_pos1, col = cols1)
segments(ci2[,1], y_pos2 - 0.4, ci2[,2], y_pos2 - 0.4, col = cols2)

# Axe Y avec noms des variables
axis(2, at = c(4.8,3.8,2.8,1.8,1,0.6), labels = c("SD","TD","SLA","LA","LDMC","LMC_t24"), las = 1,cex.axis = 0.9)

# Axe X
axis(1,cex.axis = 0.9)

# Valeurs des coefficients + étoiles
text(estimates1, y_pos1 + 0.16,
     labels = paste0(round(estimates1, 2), stars1),
     col = cols1, font = 2, cex =  0.9)
text(estimates2, y_pos2 -0.24,
     labels = paste0(round(estimates2, 2), stars2),
     col = cols2, font = 2, cex =  0.9)


dev.off()

################# Prédiction ########################
# LDMC
moy <- mean(BDD_ech3$LDMC, na.rm = TRUE)
ecart <- sd(BDD_ech3$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_ech3$LDMC_cr), max(BDD_ech3$LDMC_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_ech3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_ech3$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_ech3$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_ech3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_BB1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_ech3$LDMC, BDD_ech3$BB_prop,type="n", 
     xlab = "LDMC", ylab = "Biomasse brûlée (%)", 
     main = "Effet de LDMC sur BB ",ylim=c(0,1),xlim=c(175,550))

# Courbes de prédiction
lines((foo*ecart + moy), pmin(pmax(pred1$fit, 0), 1), col = "black", lwd = 2 )


# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)


########## LMC_t24
moy <- mean(BDD_ech3$LMC_t24, na.rm = TRUE)
ecart <- sd(BDD_ech3$LMC_t24, na.rm = TRUE)

# Valeurs de LMC_t24
foo <- seq(min(BDD_ech3$LMC_t24_cr), max(BDD_ech3$LMC_t24_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_ech3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_ech3$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_ech3$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_ech3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_BB3, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_ech3$LMC_t24, BDD_ech3$BB_prop,type="n", 
     xlab = "LMC_t24 (%)", ylab = "Biomasse brûlée (%)", 
     main = "Effet de LMC_t24 sur BB ",ylim=c(0,1),xlim=c(0,350))

# Courbes de prédiction
lines((foo*ecart + moy), pred1$fit, col = "black", lwd = 2 )


# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)


############## on fait varier LMDC et TD
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_echMT$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_echMT$SLA_cr), length(foo))
)

pred_data2 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(quantile(BDD_echMT$TD_cr,0.75), length(foo)),
  LA_cr = rep(mean(BDD_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_echMT$SLA_cr), length(foo))
)

pred_data3 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(quantile(BDD_echMT$TD_cr,0.25), length(foo)),
  LA_cr = rep(mean(BDD_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_echMT$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_BB1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred2 <- predict(m_BB1, type = "response", newdata = pred_data2, se.fit = TRUE, re.form = NA)
pred3 <- predict(m_BB1, type = "response", newdata = pred_data3, se.fit = TRUE, re.form = NA)


# Plot des points observés
plot(BDD_ech3$LDMC, BDD_ech3$BB_prop,type="n", 
     xlab = "LDMC", ylab = "Burnt Biomass", 
     main = "Effet de LDMC sur BB selon TD")

# Courbes de prédiction
lines((foo*ecart + moy), pmin(pmax(pred1$fit, 0), 1), col = "black", lwd = 2)
lines((foo*ecart+moy), pmin(pmax(pred1$fit, 0), 1), col = "red" ,lty=3) 
lines((foo*ecart+moy), pmin(pmax(pred1$fit, 0), 1), col = "red",lty=3) 


# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)







#######################################
###### BT #########
######################################

#modèles
m_BT0 <- glmer(BT_test ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr+VPD_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_ech3)  #meilleur modèle
summary(m_BT0)
r.squaredGLMM(m_BT0)
m_BT1 <- glmer(BT_test ~ SD_cr+TD_cr+SLA_cr+LA_cr+LDMC_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_ech3)  #meilleur modèle
summary(m_BT1)
r.squaredGLMM(m_BT1)
m_BT2 <- glmer(BT_test ~ SD_cr+TD_cr+LA_cr+LMC_t24_cr+LT_cr+VPD_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_ech3)  #meilleur modèle
summary(m_BT2)
r.squaredGLMM(m_BT2)
m_BT3 <- glmer(BT_test ~ SD_cr+TD_cr+SLA_cr+LA_cr+LMC_t24_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_ech3)  #meilleur modèle
summary(m_BT3)
r.squaredGLMM(m_BT3)

AIC(m_BT1,m_BT3)

###################### graph coefs ###############

# Extraire coefficients (sans l'intercept)
par(mar = c(5,10,5,5))
coefs1 <- summary(m_BT1)$coefficients[-1, ]
coefs2 <- summary(m_BT3)$coefficients[-1, ]

# Variables utiles
estimates1 <- coefs1[, "Estimate"]
stderr1 <- coefs1[, "Std. Error"]
pval1 <- coefs1[, "Pr(>|z|)"]
labels1 <- rownames(coefs1)

estimates2 <- coefs2[, "Estimate"]
stderr2 <- coefs2[, "Std. Error"]
pval2 <- coefs2[, "Pr(>|z|)"]
labels2 <- rownames(coefs2)

# Calcul des intervalles de confiance plus ou moins SE
ci1 <- cbind(estimates1 - stderr1, estimates1 + stderr1)
ci2 <- cbind(estimates2 - stderr2, estimates2 + stderr2)

# Étoiles de significativité
stars1 <- ifelse(pval1 < 0.001, "***",
                 ifelse(pval1 < 0.01, "**",
                        ifelse(pval1 < 0.05, "*",
                               ifelse(pval1 < 0.1, ".", ""))))

stars2 <- ifelse(pval2 < 0.001, "***",
                 ifelse(pval2 < 0.01, "**",
                        ifelse(pval2 < 0.05, "*",
                               ifelse(pval2 < 0.1, ".", ""))))


# Couleurs selon signe du coefficient
cols1 <- ifelse(estimates1 < 0, "red", "red")
cols2 <- ifelse(estimates2 < 0, "blue", "blue")

# Ordre des variables (du bas vers le haut)
y_pos1 <- length(estimates1):1
y_pos2 <- length(estimates2):1

# Plot de base
par(mar=c(4,6,2,2))
plot(estimates1, y_pos1,type = "n",
     xlim = c(-0.5,0.5),
     ylim = c(0.7,length(estimates1) +0.2 ),
     xlab = "Estimate",
     ylab = "",
     axes = FALSE,
     main = "Temps de combustion",
     cex.main = 1.1)



# Ligne verticale à zéro
abline(v = 0, lty = 2)
abline(h = y_pos1 -0.2, lwd = 0.5, lty = 3, col="grey")


points(estimates1, y_pos1 ,pch = 16, col = cols1,cex = 1.1)
points(estimates2, y_pos2 - 0.4,pch = 17, col = cols2,cex = 1.1)

# Barres d'erreur (IC)
segments(ci1[,1], y_pos1, ci1[,2], y_pos1, col = cols1)
segments(ci2[,1], y_pos2 - 0.4, ci2[,2], y_pos2 - 0.4, col = cols2)

# Axe Y avec noms des variables
axis(2, at = c(4.8,3.8,2.8,1.8,1,0.6), labels = c("SD","TD","SLA","LA","LDMC","LMC_t24"), las = 1,cex.axis = 0.9)

# Axe X
axis(1,cex.axis = 0.9)

# Valeurs des coefficients + étoiles
text(estimates1, y_pos1 + 0.16,
     labels = paste0(round(estimates1, 2), stars1),
     col = cols1, font = 2, cex =  0.9)
text(estimates2, y_pos2 -0.24,
     labels = paste0(round(estimates2, 2), stars2),
     col = cols2, font = 2, cex =  0.9)




################# Prédiction ########################
# LDMC
moy <- mean(BDD_ech3$LDMC, na.rm = TRUE)
ecart <- sd(BDD_ech3$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_ech3$LDMC_cr), max(BDD_ech3$LDMC_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_ech3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_ech3$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_ech3$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_ech3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_BT1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred1

# Plot des points observés
plot(BDD_ech3$LDMC, BDD_ech3$BT_test,type="n", 
     xlab = "LDMC", ylab = "Temps de combustion (s)", 
     main = "Effet de LDMC sur BT ",ylim=c(0,80),xlim=c(175,550))

# Courbes de prédiction
lines((foo*ecart + moy), pred1$fit, col = "black", lwd = 2 )

# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), (pred1$fit + 1.96 * pred1$se.fit), col = "blue", lty = 3)
lines((foo*ecart + moy), (pred1$fit - 1.96 * pred1$se.fit), col = "blue", lty = 3)

range(pred1$se.fit)
range(pred1$fit + 1.96 * pred1$se.fit - (pred1$fit - 1.96 * pred1$se.fit))



########## LMC_t24
moy <- mean(BDD_ech3$LMC_t24, na.rm = TRUE)
ecart <- sd(BDD_ech3$LMC_t24, na.rm = TRUE)

# Valeurs de LMC_t24
foo <- seq(min(BDD_ech3$LMC_t24_cr), max(BDD_ech3$LMC_t24_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_ech3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_ech3$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_ech3$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_ech3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_BT3, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_ech3$LMC_t24, BDD_ech3$BT,type="n", 
     xlab = "LMC_t24 (%)", ylab = "Temps de combustion (s)", 
     main = "Effet de LMC_t24 sur BT ",ylim=c(0,80),xlim=c(0,350))

# Courbes de prédiction
lines((foo*ecart + moy), pred1$fit, col = "black", lwd = 2 )


# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)





############## on fait varier LMDC et TD
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_echMT$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_echMT$SLA_cr), length(foo))
)

pred_data2 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(quantile(BDD_echMT$TD_cr,0.75), length(foo)),
  LA_cr = rep(mean(BDD_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_echMT$SLA_cr), length(foo))
)

pred_data3 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(quantile(BDD_echMT$TD_cr,0.25), length(foo)),
  LA_cr = rep(mean(BDD_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_echMT$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_BT1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred2 <- predict(m_BT1, type = "response", newdata = pred_data2, se.fit = TRUE, re.form = NA)
pred3 <- predict(m_BT1, type = "response", newdata = pred_data3, se.fit = TRUE, re.form = NA)


# Plot des points observés
plot(BDD_ech3$LDMC, BDD_ech3$BT,type="n",ylim=c(0,80), 
     xlab = "LDMC", ylab = "Temps de combustion", 
     main = "Effet de LDMC sur BT selon TD")

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 
lines((foo*ecart+moy), pred2$fit, col = "red" ,lty=3) 
lines((foo*ecart+moy), pred3$fit, col = "red",lty=3) 


# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), ((pred1$fit + 1.96 * pred1$se.fit)*ecart + moy), col = "blue", lty = 3)
lines((foo*ecart + moy), ((pred1$fit + 1.96 * pred1$se.fit)*ecart + moy), col = "blue", lty = 3)






###################
###### DI ######################################################################
###################

#modèles

m_DI0 <- glmer(DI_test ~ SD_cr+TD_cr+LA_cr+LDMC_cr+LT_cr+VPD_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_ech3)  #meilleur modèle
summary(m_DI0)
r.squaredGLMM(m_DI0)
m_DI1 <- glmer(DI_test ~ SD_cr+TD_cr+SLA_cr+LA_cr+LDMC_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_ech3)  #meilleur modèle
summary(m_DI1)
r.squaredGLMM(m_DI1)
m_DI2 <- glmer(DI_test ~ SD_cr+TD_cr+LA_cr+LMC_t24_cr+LT_cr+VPD_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_ech3)  #meilleur modèle
summary(m_DI2)
r.squaredGLMM(m_DI2)
m_DI3 <- glmer(DI_test ~ SD_cr+TD_cr+SLA_cr+LA_cr+LMC_t24_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_ech3)  #meilleur modèle
summary(m_DI3)
r.squaredGLMM(m_DI3)



AIC(m_DI1,m_DI3)

###################### graph coefs ###############

# Extraire coefficients (sans l'intercept)
par(mar = c(5,10,5,5))
coefs1 <- summary(m_DI1)$coefficients[-1, ]
coefs2 <- summary(m_DI3)$coefficients[-1, ]

# Variables utiles
estimates1 <- coefs1[, "Estimate"]
stderr1 <- coefs1[, "Std. Error"]
pval1 <- coefs1[, "Pr(>|z|)"]
labels1 <- rownames(coefs1)

estimates2 <- coefs2[, "Estimate"]
stderr2 <- coefs2[, "Std. Error"]
pval2 <- coefs2[, "Pr(>|z|)"]
labels2 <- rownames(coefs2)

# Calcul des intervalles de confiance plus ou moins SE
ci1 <- cbind(estimates1 - stderr1, estimates1 + stderr1)
ci2 <- cbind(estimates2 - stderr2, estimates2 + stderr2)

# Étoiles de significativité
stars1 <- ifelse(pval1 < 0.001, "***",
                 ifelse(pval1 < 0.01, "**",
                        ifelse(pval1 < 0.05, "*",
                               ifelse(pval1 < 0.1, ".", ""))))

stars2 <- ifelse(pval2 < 0.001, "***",
                 ifelse(pval2 < 0.01, "**",
                        ifelse(pval2 < 0.05, "*",
                               ifelse(pval2 < 0.1, ".", ""))))


# Couleurs selon signe du coefficient
cols1 <- ifelse(estimates1 < 0, "red", "red")
cols2 <- ifelse(estimates2 < 0, "blue", "blue")

# Ordre des variables (du bas vers le haut)
y_pos1 <- length(estimates1):1
y_pos2 <- length(estimates2):1

# Plot de base
png("C:/IRD/Stage_inflammabilit-/Figures/FIN/Coef_DI.png",width = 500, height = 500,res=150)
par(mar=c(4,6,2,2))
plot(estimates1, y_pos1,type = "n",
     xlim = c(-0.3,0.3),
     ylim = c(0.7,length(estimates1) +0.2 ),
     xlab = "Estimate",
     ylab = "",
     axes = FALSE,
     main = "Délai d'ignition",
     cex.main = 1.1)



# Ligne verticale à zéro
abline(v = 0, lty = 2)
abline(h = y_pos1 -0.2, lwd = 0.5, lty = 3, col="grey")


points(estimates1, y_pos1 ,pch = 16, col = cols1,cex = 1.1)
points(estimates2, y_pos2 - 0.4,pch = 17, col = cols2,cex = 1.1)

# Barres d'erreur (IC)
segments(ci1[,1], y_pos1, ci1[,2], y_pos1, col = cols1)
segments(ci2[,1], y_pos2 - 0.4, ci2[,2], y_pos2 - 0.4, col = cols2)

# Axe Y avec noms des variables
axis(2, at = c(4.8,3.8,2.8,1.8,1,0.6), labels = c("SD","TD","SLA","LA","LDMC","LMC_t24"), las = 1,cex.axis = 0.9)

# Axe X
axis(1,cex.axis = 0.9)

# Valeurs des coefficients + étoiles
text(estimates1, y_pos1 + 0.16,
     labels = paste0(round(estimates1, 2), stars1),
     col = cols1, font = 2, cex =  0.9)
text(estimates2, y_pos2 -0.24,
     labels = paste0(round(estimates2, 2), stars2),
     col = cols2, font = 2, cex =  0.9)


dev.off()


################# Prédiction ########################
# LDMC
moy <- mean(BDD_ech3$LDMC, na.rm = TRUE)
ecart <- sd(BDD_ech3$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_ech3$LDMC_cr), max(BDD_ech3$LDMC_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_ech3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_ech3$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_ech3$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_ech3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_DI1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_ech3$LDMC, BDD_ech3$DI_test,type="n", 
     xlab = "LDMC", ylab = "Délai d'ignition (s)", 
     main = "Effet de LDMC sur DI ",ylim=c(0.5,2),xlim=c(175,550))

# Courbes de prédiction
lines((foo*ecart + moy), pred1$fit, col = "black", lwd = 2 )


# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)


########## LMC_t24
moy <- mean(BDD_ech3$LMC_t24, na.rm = TRUE)
ecart <- sd(BDD_ech3$LMC_t24, na.rm = TRUE)

# Valeurs de LMC_t24
foo <- seq(min(BDD_ech3$LMC_t24_cr), max(BDD_ech3$LMC_t24_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_ech3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_ech3$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_ech3$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_ech3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_DI3, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_ech3$LMC_t24, BDD_ech3$DI_test,type="n", 
     xlab = "LMC_t24 (%)", ylab = "Délai d'ignition (s)", 
     main = "Effet de LMC_t24 sur DI ",ylim=c(0.5,2),xlim=c(0,350))

# Courbes de prédiction
lines((foo*ecart + moy), pred1$fit, col = "black", lwd = 2 )


# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)




# on fait varier LMDC et on fixe les autres variables
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_ech3$SD_cr,na.rm = TRUE), length(foo)),
  SLA_cr = rep(mean(BDD_ech3$SLA_cr), length(foo)),
  LA_cr = rep(mean(BDD_ech3$LA_cr), length(foo)),
  TD_cr = rep(mean(BDD_ech3$TD_cr), length(foo))
)


pred_data2 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_ech3$SD_cr,na.rm = TRUE), length(foo)),
  SLA_cr = rep(quantile(BDD_ech3$SLA_cr,0.75), length(foo)),
  LA_cr = rep(mean(BDD_ech3$LA_cr), length(foo)),
  TD_cr = rep(mean(BDD_ech3$TD_cr), length(foo))
)

pred_data3 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_ech3$SD_cr,na.rm = TRUE), length(foo)),
  SLA_cr = rep(quantile(BDD_ech3$SLA_cr,0.25), length(foo)),
  LA_cr = rep(mean(BDD_ech3$LA_cr), length(foo)),
  TD_cr = rep(mean(BDD_ech3$TD_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_DI1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred2 <- predict(m_DI1, type = "response", newdata = pred_data2, se.fit = TRUE, re.form = NA)
pred3 <- predict(m_DI1, type = "response", newdata = pred_data3, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_ech3$LDMC, BDD_ech3$DI_test,type="n", 
     xlab = "LDMC", ylab = "Ignition Delay", 
     main = "Effet de LDMC sur DI selon SLA",ylim=c(0.5,2))

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 
lines((foo*ecart+moy), pred2$fit, col = "red" ,lty=3) 
lines((foo*ecart+moy), pred3$fit, col = "red",lty=3)






#####################
###### Score ########
#####################

#modèles
m_score1 <- lmer(score ~ SD_cr+TD_cr+SLA_cr+LA_cr+LDMC_cr + (1 | Nom_scientifique), data = BDD_ech)
summary(m_score1)
r.squaredGLMM(m_score1)
m_score3 <- lmer(score ~ SD_cr+TD_cr+SLA_cr+LA_cr+LMC_t24_cr + (1 | Nom_scientifique), data = BDD_ech)
summary(m_score3)
r.squaredGLMM(m_score3)

AIC(m_score1,m_score3)



###################### graph coefs ###############

# Extraire coefficients (sans l'intercept)
par(mar = c(5,10,5,5))
coefs1 <- summary(m_score1)$coefficients[-1, ]
coefs2 <- summary(m_score3)$coefficients[-1, ]

# Variables utiles
estimates1 <- coefs1[, "Estimate"]
stderr1 <- coefs1[, "Std. Error"]
pval1 <- coefs1[, "Pr(>|t|)"]
labels1 <- rownames(coefs1)

estimates2 <- coefs2[, "Estimate"]
stderr2 <- coefs2[, "Std. Error"]
pval2 <- coefs2[, "Pr(>|t|)"]
labels2 <- rownames(coefs2)

# Calcul des intervalles de confiance plus ou moins SE
ci1 <- cbind(estimates1 - stderr1, estimates1 + stderr1)
ci2 <- cbind(estimates2 - stderr2, estimates2 + stderr2)

# Étoiles de significativité
stars1 <- ifelse(pval1 < 0.001, "***",
                 ifelse(pval1 < 0.01, "**",
                        ifelse(pval1 < 0.05, "*",
                               ifelse(pval1 < 0.1, ".", ""))))

stars2 <- ifelse(pval2 < 0.001, "***",
                 ifelse(pval2 < 0.01, "**",
                        ifelse(pval2 < 0.05, "*",
                               ifelse(pval2 < 0.1, ".", ""))))


# Couleurs selon signe du coefficient
cols1 <- ifelse(estimates1 < 0, "red", "red")
cols2 <- ifelse(estimates2 < 0, "blue", "blue")

# Ordre des variables (du bas vers le haut)
y_pos1 <- length(estimates1):1
y_pos2 <- length(estimates2):1

# Plot de base
par(mar=c(4,6,2,2))
plot(estimates1, y_pos1,type = "n",
     xlim = c(-2,2),
     ylim = c(0.7,length(estimates1) +0.2 ),
     xlab = "Estimate",
     ylab = "",
     axes = FALSE,
     main = "Score d'inflammabilité",
     cex.main = 1.1)



# Ligne verticale à zéro
abline(v = 0, lty = 2)
abline(h = y_pos1 -0.2, lwd = 0.5, lty = 3, col="grey")


points(estimates1, y_pos1 ,pch = 16, col = cols1,cex = 1.1)
points(estimates2, y_pos2 - 0.4,pch = 17, col = cols2,cex = 1.1)

# Barres d'erreur (IC)
segments(ci1[,1], y_pos1, ci1[,2], y_pos1, col = cols1)
segments(ci2[,1], y_pos2 - 0.4, ci2[,2], y_pos2 - 0.4, col = cols2)

# Axe Y avec noms des variables
axis(2, at = c(4.8,3.8,2.8,1.8,1,0.6), labels = c("SD","TD","SLA","LA","LDMC","LMC_t24"), las = 1,cex.axis = 0.9)

# Axe X
axis(1,cex.axis = 0.9)

# Valeurs des coefficients + étoiles
text(estimates1, y_pos1 + 0.16,
     labels = paste0(round(estimates1, 2), stars1),
     col = cols1, font = 2, cex =  0.9)
text(estimates2, y_pos2 -0.24,
     labels = paste0(round(estimates2, 2), stars2),
     col = cols2, font = 2, cex =  0.9)




################# Prédiction ########################
# LDMC
moy <- mean(BDD_ech$LDMC, na.rm = TRUE)
ecart <- sd(BDD_ech$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_ech$LDMC_cr), max(BDD_ech$LDMC_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_ech$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_ech$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_ech$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_ech$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_score1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_ech$LDMC, BDD_ech$score,type="n", 
     xlab = "LDMC", ylab = "Délai d'ignition (s)", 
     main = "Effet de LDMC sur score ",ylim=c(-5,3),xlim=c(100,550))

# Courbes de prédiction
lines((foo*ecart + moy), pred1$fit, col = "black", lwd = 2 )


# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)


########## LMC_t24
moy <- mean(BDD_ech$LMC_t24, na.rm = TRUE)
ecart <- sd(BDD_ech$LMC_t24, na.rm = TRUE)

# Valeurs de LMC_t24
foo <- seq(min(BDD_ech$LMC_t24_cr), 5,length.out = 100)
max(BDD_ech$LMC_t24_cr)
# on fait varier LMDC 
pred_data1 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_ech$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_ech$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_ech$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_ech$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_score3, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_ech$LMC_t24, BDD_ech$score,type="n", 
     xlab = "LMC_t24 (%)", ylab = "Délai d'ignition (s)", 
     main = "Effet de LMC_t24 sur score ",ylim=c(-5,3),xlim=c(0,700))

# Courbes de prédiction
lines((foo*ecart + moy), pred1$fit, col = "black", lwd = 2 )


# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)











###################################################################################
######################### ANALYSE GMIN #########################################
###################################################################################

hist(BDD_esp$LMC_t24)
hist(BDD_esp$LMC_t0)
hist(BDD_esp$Gmin)
BDD_esp$ratio <- BDD_esp$LMC_t24 / BDD_esp$LMC_t0
hist(BDD_esp$ratio)
BDD_esp$logLMC24<-log(BDD_esp$LMC_t24)
BDD_esp$logLMC0<-log(BDD_esp$LMC_t0)

BDD_esp$LMC_t0_cr<-as.numeric(scale(BDD_esp$LMC_t0))
BDD_esp$Gmin_cr<-as.numeric(scale(BDD_esp$Gmin))


m_gmin2<-glm(logLMC24~logLMC0+Gmin_cr,family=gaussian,data=BDD_esp)
summary(m_gmin2)
m_gmin2<-glm(LMC_t24~LMC_t0_cr+Gmin_cr,family=Gamma(link = "log"),data=BDD_esp)
summary(m_gmin2)
m_gmin3<-glm(ratio~LMC_t0_cr+Gmin_cr,family=gaussian,data=BDD_esp)
summary(m_gmin3)

AIC(m_gmin2,m_gmin3)

hist(resid(m_gmin2))

################# Prédiction ########################
# LMC_t0
moy <- mean(BDD_esp$LMC_t0, na.rm = TRUE)
ecart <- sd(BDD_esp$LMC_t0, na.rm = TRUE)

# Valeurs de LMC_t0
foo <- seq(min(BDD_esp$LMC_t0_cr), max(BDD_esp$LMC_t0_cr),length.out = 100)

# on fait varier LMDC et on fixe les autres variables
pred_data1 <- data.frame(
  LMC_t0_cr = foo,
  Gmin_cr = rep(mean(BDD_esp$Gmin_cr,na.rm = TRUE), length(foo))
)

pred_data2 <- data.frame(
  LMC_t0_cr = foo,
  Gmin_cr = rep(min(BDD_esp$Gmin_cr,na.rm = TRUE), length(foo))
)

pred_data3 <- data.frame(
  LMC_t0_cr = foo,
  Gmin_cr = rep(max(BDD_esp$Gmin_cr,na.rm = TRUE), length(foo))
)

# Prédictions
pred1 <- predict(m_gmin3, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred2 <- predict(m_gmin3, type = "response", newdata = pred_data2, se.fit = TRUE, re.form = NA)
pred3 <- predict(m_gmin3, type = "response", newdata = pred_data3, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_esp$LMC_t0, BDD_esp$LMC_t24, 
     xlab = "LMC_t0 (%)", ylab = "LMC_t24 (%)", 
     main = "Effet de LMC_t0 sur LMC_t24 selon gmin ",ylim=c(0,700))


# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit*(foo*ecart+moy), col = "#3B393B", lwd = 2) 
lines((foo*ecart+moy), pred2$fit*(foo*ecart+moy), col = "#005ED1" ,lwd = 2) 
lines((foo*ecart+moy), pred3$fit*(foo*ecart+moy), col = "#BD0000",lwd = 2)

# Prédictions
pred1 <- predict(m_gmin2, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred2 <- predict(m_gmin2, type = "response", newdata = pred_data2, se.fit = TRUE, re.form = NA)
pred3 <- predict(m_gmin2, type = "response", newdata = pred_data3, se.fit = TRUE, re.form = NA)


# Créer le graphique
plot(BDD_esp$LMC_t0, BDD_esp$LMC_t24, 
     xlab = "LMC_t0 (%)", ylab = "LMC_t24 (%)", 
     main = "Effet de LMC_t0 sur LMC_t24 selon gmin",
     xlim = c(100,800), ylim = c(0,800))

# Ajouter les courbes de prédiction
lines((foo * ecart + moy), pred1$fit, col = "#3B393B", lwd = 2) 
lines((foo * ecart + moy), pred2$fit, col = "#005ED1", lwd = 2) 
lines((foo * ecart + moy), pred3$fit, col = "#BD0000", lwd = 2)



polygon(x = c(-50, -50, 850), 
        y = c(850, -50, 850), 
        density = 20, angle = -45, col = "#615B5B", border = NA)

# Ajouter la ligne y = x
abline(a = 0, b = 1)

esp_a_annoter <- BDD_esp$Nom_scientifique == "Scaevola taccada" | 
  BDD_esp$Nom_scientifique == "Scaevola montana" |
  BDD_esp$Nom_scientifique== "Barringtonia asiatica"  

points (BDD_esp$LMC_t0[esp_a_annoter],
        BDD_esp$LMC_t24[esp_a_annoter],pch=21,bg = rgb(1,0,0),cex=1.2
)

text(BDD_esp$LMC_t0[esp_a_annoter] -10 ,
     BDD_esp$LMC_t24[esp_a_annoter] + 25,
     labels = c("B. asiatica","S. montana","S. taccada"),
     cex = 1)










# ZOOM
plot(BDD_esp$LMC_t0, BDD_esp$LMC_t24, 
     xlab = "LMC_t0 (%)", ylab = "LMC_t24 (%)", 
     main = "Effet de LMC_t0 sur LMC_t24 selon gmin",
     xlim = c(70,410), ylim = c(0,400))

# Ajouter les courbes de prédiction
lines((foo * ecart + moy), pred1$fit, col = "#3B393B", lwd = 2) 
lines((foo * ecart + moy), pred2$fit, col = "#005ED1", lwd = 2) 
lines((foo * ecart + moy), pred3$fit, col = "#BD0000", lwd = 2)

# Ajouter le polygone hachuré (zone où y > x)
polygon(x = c(-50, -50, 850), 
        y = c(850, -50, 850), 
        density = 15, angle = -45, col = "#615B5B", border = NA)

# Ajouter la ligne y = x
abline(a = 0, b = 1)

esp_a_annoter <- BDD_esp$Nom_scientifique == "Scaevola montana" | 
  BDD_esp$Nom_scientifique == "Hibiscus tiliaceus R" |
  BDD_esp$Nom_scientifique== "Acropogon bullatus" | 
  BDD_esp$Nom_scientifique== "Barringtonia asiatica" | 
  BDD_esp$Nom_scientifique== "Alphitonia neocaledonica" |
  BDD_esp$Nom_scientifique== "Leucaena leucocephala"

points (BDD_esp$LMC_t0[esp_a_annoter],
        BDD_esp$LMC_t24[esp_a_annoter],pch=21,bg = rgb(1,0,0),cex=1.2
)

text(BDD_esp$LMC_t0[esp_a_annoter] + 2,
     BDD_esp$LMC_t24[esp_a_annoter] + 15,
     labels = c("A. bullatus","A. neocaledonica","B. asiatica","H. tiliaceus R","L. leucocephala","S. montana"),
      cex = 1)






# Plot des points observés
plot(BDD_esp$LMC_t0, BDD_esp$LMC_t24, 
     xlab = "LMC_t0 (%)", ylab = "LMC_t24 (%)", 
     main = "Effet de LMC_t0 sur LMC_t24 selon gmin ",ylim=c(0,700))


# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "#3B393B", lwd = 2) 
lines((foo*ecart+moy), pred2$fit, col = "#005ED1" ,lwd = 2) 
lines((foo*ecart+moy), pred3$fit, col = "#BD0000",lwd = 2)

abline(a=0,b=1)
abline



########################## PLOT GMIN LMC_t0 ###########################

M_L2<-aggregate(BDD_ech$LMC_t0,by=list(BDD_ech$Nom_scientifique),mean)
SD_L2<-aggregate(BDD_ech$LMC_t0,by=list(BDD_ech$Nom_scientifique),sd)

M_L1<-aggregate(BDD_ech$Gmin,by=list(BDD_ech$Nom_scientifique),mean)
SD_L1<-aggregate(BDD_ech$Gmin,by=list(BDD_ech$Nom_scientifique),sd)



# Création du scatterplot
plot(BDD_esp$LMC_t0 ~ BDD_esp$Gmin, type="n",
     xlab = "Gmin", 
     ylab = "LMC_t0 (%)", 
     main = "Relation entre Gmin et LMC_t0",ylim=c(0,450),xlim=c(0,100))

esp_a_annoter <- BDD_esp$Nom_scientifique == "Scaevola montana" | 
  BDD_esp$Nom_scientifique == "Hibiscus tiliaceus R" |
  BDD_esp$Nom_scientifique== "Acropogon bullatus" |
  BDD_esp$Nom_scientifique== "Leucaena leucocephala" | 
  BDD_esp$Nom_scientifique== "Cerbera manghas" | 
  BDD_esp$Nom_scientifique== "Diospyros fasciculosa" |
  BDD_esp$Nom_scientifique== "Calophylum caledonicum" |
  BDD_esp$Nom_scientifique== "Alphitonia neocaledonica" |

segments(x0=BDD_esp$Gmin[esp_a_annoter],x1=BDD_esp$Gmin[esp_a_annoter],y0=BDD_esp$LMC_t0[esp_a_annoter]-SD_L2$x[esp_a_annoter],y1=BDD_esp$LMC_t0[esp_a_annoter]+SD_L2$x[esp_a_annoter])
segments(x0=BDD_esp$Gmin[esp_a_annoter]-SD_L1$x[esp_a_annoter],x1=BDD_esp$Gmin[esp_a_annoter]+SD_L1$x[esp_a_annoter],y0=BDD_esp$LMC_t0[esp_a_annoter],y1=BDD_esp$LMC_t0[esp_a_annoter])
points(M_L1$x,M_L2$x,pch=19,cex=0.8)
points (BDD_esp$Gmin[esp_a_annoter],
        BDD_esp$LMC_t0[esp_a_annoter],pch=21,bg = rgb(1,0,0),cex=1.2
)


text(BDD_esp$Gmin[esp_a_annoter],
     BDD_esp$LMC_t0[esp_a_annoter],
     labels = BDD_esp$Nom_scientifique[esp_a_annoter],
     pos = 4, offset = 0.5, cex = 0.8)

points (BDD_esp$Gmin[esp_a_annoter],
         BDD_esp$LMC_t0[esp_a_annoter],pch=21,bg = rgb(1,0,0)
         )

plot(BDD_esp$LMC_t0,BDD_esp$LMC_t24 ,xlim=c(0,300),ylim=c(0,300))
text(BDD_esp$LMC_t0,BDD_esp$LMC_t24,
     labels = BDD_esp$Nom_scientifique,cex=1,pos=2)
