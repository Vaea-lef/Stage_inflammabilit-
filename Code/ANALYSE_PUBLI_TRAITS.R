#### ANALYSE DE DONNEE PUBLICATION TRAITS ####





################# IMPORTATION ET CHARGEMENT PACKAGE ###########################


####PACKAGE

library(FactoMineR)
library(factoextra)
library (corrplot)
library(lmerTest)
library(lme4)
library(MuMIn)
library(ape)
library(nlme)
library(caper)
library(phytools)
library(geiger)
library(rotl)
library(dplyr)


#### IMPORT BDD

setwd("C:/IRD/Stage_inflammabilit-") #définition du répertoire de travail

##  importation BDD_ana_ech en format CSV
BDD_moy_ech<-read.csv2("Data/Publi/BDD_moy_ech_publi.csv", header = TRUE) #importation de l
names(BDD_moy_ech)[which(names(BDD_moy_ech) == "Nb_ramifications")] <- "Nb_rami"

BDD_moy_echMT<-read.csv2("Data/Publi/BDD_moy_echMT.csv", header = TRUE) #importation de l
names(BDD_moy_echMT)[which(names(BDD_moy_echMT) == "Nb_ramifications")] <- "Nb_rami"

BDD_moy_ech3<-read.csv2("Data/Publi/BDD_moy_ech3.csv", header = TRUE) #importation de l
names(BDD_moy_ech3)[which(names(BDD_moy_ech3) == "Nb_ramifications")] <- "Nb_rami"

##  importation BDD_sd_esp en format CSV
BDD_sd_esp<-read.csv2("Data/Publi/BDD_sd_esp.csv", header = TRUE)

##  importation BDD_moy_esp en format CSV
BDD_moy_esp<-read.csv2("Data/Publi/BDD_moy_esp.csv", header = TRUE) #importation de la base
names(BDD_moy_esp)[which(names(BDD_moy_esp) == "Nb_ramifications")] <- "Nb_rami"

##  importation BDD_moy_esp spéciale MT en format CSV
BDD_moy_espMT<-read.csv2("Data/Publi/BDD_moy_espMT.csv", header = TRUE) #importation de la base
names(BDD_moy_espMT)[which(names(BDD_moy_espMT) == "Nb_ramifications")] <- "Nb_rami"

##  importation BDD_moy_esp spéciale BT, BB, DI en format CSV
BDD_moy_esp3<-read.csv2("Data/Publi/BDD_moy_esp3.csv", header = TRUE) #importation de la base
names(BDD_moy_esp3)[which(names(BDD_moy_esp3) == "Nb_ramifications")] <- "Nb_rami"

##calcul FI

# Comptage du nombre d'essais par espèce dans BDD_moy_ech
essais_par_espece <- table(BDD_moy_ech$Nom_scientifique)
essais_par_espece

# Création d'une nouvelle colonne Nb_essais dans BDD_moy_esp en s'appuyant sur le nom scientifique
BDD_moy_esp$Nb_essais <- essais_par_espece[BDD_moy_esp$Nom_scientifique]

# Vérification
head(BDD_moy_esp)


head(BDD_moy_esp)
BDD_moy_esp$BD_mean <- BDD_moy_esp$BD           # copier la colonne originale
BDD_moy_esp$BD_mean[is.na(BDD_moy_esp$BD_mean)] <- mean(BDD_moy_esp$BD, na.rm = TRUE)
head(BDD_moy_esp)


head(BDD_moy_ech)
BDD_moy_ech$BD_mean <- BDD_moy_ech$BD           # copier la colonne originale
BDD_moy_ech$BD_mean[is.na(BDD_moy_ech$BD_mean)] <- mean(BDD_moy_ech$BD, na.rm = TRUE)
head(BDD_moy_ech)














############################################################################
#################### ACP TRAITS (pour sélection des traits modèle)############################
############################################################################


# Sélection des colonnes des traits fonctionnels
colonnes_traits <- BDD_moy_ech [, c("Nb_rami","Gmin","TD","TDMC","Surface_F","LDMC","SLA","LT","BD_mean")]

#strandardiser les données
colonnes_traits_cr <- scale(log(colonnes_traits)+1)

# Application de l'ACP
res.pca <- PCA(colonnes_traits, scale.unit = TRUE, graph = FALSE)

# Résumé des résultats
summary(res.pca)

# Graphique des contributions des variables aux composantes principales
fviz_pca_var(res.pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

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
############################ ECHELLE ECHANTILLON ###############################
################################################################################

############################### ACP INFLA (avec projection des traits) ######################################

# Sélection des colonnes des composantes de l'inflammabilité
colonnes_infla <- BDD_moy_ech[, c("score_DI", "BB_test", "BT", "MT")]
head(colonnes_infla)

#sélection de variabes supplémentaires (traits)
colonnes_traits <- BDD_moy_ech [, c("Nb_rami","Gmin","TD","TDMC","LDMC","SLA","BD_mean")]
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

M_S1<-aggregate(res.pca$x[,1],by=list(BDD_moy_ech$Nom_scientifique),mean)
SD_S1<-aggregate(res.pca$x[,1],by=list(BDD_moy_ech$Nom_scientifique),sd)

M_S2<-aggregate(res.pca$x[,2],by=list(BDD_moy_ech$Nom_scientifique),mean)
SD_S2<-aggregate(res.pca$x[,2],by=list(BDD_moy_ech$Nom_scientifique),sd)

# Coordonnées
ind_coords <- res.pca$x                    # individus
var_coords <- res.pca$rotation             # variables
eig_vals <- res.pca$sdev^2                 # valeurs propres
explained_var <- round(100 * eig_vals / sum(eig_vals), 1)

colonnes_taits_coord <- cor(colonnes_traits_cr, ind_coords)

# Ajout d'un score d'inflammabilité basé sur la coordonnée de l'axe 1 
BDD_moy_ech$score <- ind_coords[,1]


# Graphique de base
par(mar = c(5,5,5,5))  # marges

HC<-hclust(d=dist(cbind(M_S1$x,M_S2$x)),method="ward.D2")
plot(HC, hang = -1,labels=F, axes="n")
axis(2,cex.axis=0.6)
GR<-cutree(HC,k=5)
GR

COL<-character()
COL[GR==5]<-"#FF2A00"
COL[GR==4]<-"#CFF200"
COL[GR==3]<-"#FFE100"
COL[GR==2]<-"#6DC700"
COL[GR==1]<-"#FF8900"


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

text(var_coords[,1]*max(abs(ind_coords[,1])),
     var_coords[,2]*max(abs(ind_coords[,2])),labels = colnames(colonnes_infla_cr),
     col = "red", cex = 0.8)

# Variables supplémentaires en bleu
arrows(0, 0,
       colonnes_taits_coord[,1]*max(abs(ind_coords[,1])),
       colonnes_taits_coord[,2]*max(abs(ind_coords[,2])),
       length = 0.1, col = "#5490FF",lwd=0.5)

text(colonnes_taits_coord[,1]*max(abs(ind_coords[,1])),
     colonnes_taits_coord[,2]*max(abs(ind_coords[,2])),
     labels = colnames(colonnes_traits_cr),
     col = "blue", cex = 0.8)


# Tracer les individus
plot(ind_coords[,1], ind_coords[,2],
     xlim = range(ind_coords[,1]) * 1.2,
     ylim = c(-2,4),
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
BDD_moy_esp$groupe<-GR


# Reclassement de la variable groupe 
BDD_moy_esp$groupe <- factor(BDD_moy_esp$groupe, levels = c(2,4, 3, 1,5))
write.csv2(BDD_moy_esp,"Data/BDD_moy_esp_groupe.csv")

# boxplot MT
boxplot(MT ~ groupe, data = BDD_moy_esp,
        main = "Distribution de MT par groupe",
        xlab = "Groupe",
        ylab = "MT (°C)",col=c("#025E00","#6DC700","#CFF200","#FFE100","#FF8900","#FF2A00"))

# boxplot BT
boxplot(BT_test ~ groupe, data = BDD_moy_esp,
        main = "Distribution de BT par groupe",
        xlab = "Groupe",
        ylab = "BT (s)",col=c("#025E00","#6DC700","#CFF200","#FFE100","#FF8900","#FF2A00"))


# boxplot BB
boxplot(BB_test ~ groupe, data = BDD_moy_esp,
        main = "Distribution de BB par groupe",
        xlab = "Groupe",
        ylab = "BB (%)",col=c("#025E00","#6DC700","#CFF200","#FFE100","#FF8900","#FF2A00"))

# boxplot DI
boxplot(DI_test ~ groupe, data = BDD_moy_esp,
        main = "Distribution de DI par groupe",
        xlab = "Groupe",
        ylab = "DI (s)",col=c("#025E00","#6DC700","#CFF200","#FFE100","#FF8900","#FF2A00"))










#Calcul de la moyenne du score pour ajout dans BDD_moy_esp
#création de table avec moyenne et sd pour chaque variable en fonction du nom de l'espèce
tem3<-BDD_moy_ech[,7:36] ###sélection des colonnes comprenant les variables pour les intégrer dans la boucle
tem3
#création d'un bdd d'origine pour moyenne (sert pour merge)
BDD_moy_score <- aggregate(tem3[,1] ~ Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_moy_ech, FUN = mean, na.rm = TRUE)
BDD_moy_score[,6] <- round(BDD_moy_score[,6], 2)
colnames(BDD_moy_score)[6] <- colnames(tem3)[1]

#création d'un bdd d'origine pour sd (sert pour merge)
BDD_sd_score <- aggregate(tem3[,1] ~ Nom_scientifique + Genre + Espece + ID_espece + Milieu_recolte, data = BDD_moy_ech, FUN = sd, na.rm = TRUE)
BDD_sd_score[,6] <- round(BDD_sd_score[,6], 2)
colnames(BDD_sd_score)[6] <- colnames(tem3)[1]

#Boucle pour les calcul des moyennes et écart-types
for (i in 2:ncol(tem3)) {
  
  # Moyenne
  tem3_moy_esp <- aggregate(tem3[, i] ~ Nom_scientifique + Genre + Espece + ID_espece+ Milieu_recolte, data = BDD_moy_ech, FUN = mean, na.rm = TRUE)
  tem3_moy_esp[,6] <- round(tem3_moy_esp[,6], 2)
  colnames(tem3_moy_esp)[6] <- colnames(tem3)[i]
  BDD_moy_score <- merge(BDD_moy_score, tem3_moy_esp, by = c("Nom_scientifique", "Genre" ,"Espece","ID_espece", "Milieu_recolte"), all = TRUE)
  
  # Ecart-type
  tem3_sd_esp <- aggregate(tem3[, i] ~ Nom_scientifique + Genre + Espece + ID_espece+ Milieu_recolte, data = BDD_moy_ech, FUN = sd, na.rm = TRUE)
  tem3_sd_esp[,6] <- round(tem3_sd_esp[,6], 2)
  colnames(tem3_sd_esp)[6] <- colnames(tem3)[i]
  BDD_sd_score <- merge(BDD_sd_score, tem3_sd_esp, by = c("Nom_scientifique", "Genre" ,"Espece","ID_espece", "Milieu_recolte"), all = TRUE)
}

BDD_moy_score
write.csv2(BDD_moy_score,"Data/Publi/BDD_moy_score.csv")


BDD_moy_score




####################### CLASSEMENT ESPECES ############################
par(mar = c(4,13,0,3))
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
plot(BDD_moy_score$score[o], 1:length(BDD_moy_score$Nom_scientifique), axes="n",xlim = c(-6,4))
# Axes
axis(2, at = 1:length(BDD_moy_score$Nom_scientifique), labels = BDD_moy_score$Nom_scientifique[o], las = 1, cex.axis = 0.65,font=3)
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














###########################################
############# ARBRE PHYLO #################
###########################################


# Chargement du fichier de l'arbre paftol
tree <- read.tree("arbre.tree")

# 1. création de nouveaux noms SANS encore les assigner à l'arbre
nouveaux_noms <- sapply(strsplit(tree$tip.label, "_"), function(x) paste(x[3], x[4], sep = " "))

# 2. identification des tips en doublons à supprimer (garder le premier, supprimer les suivants)
tips_a_supprimer <- tree$tip.label[duplicated(nouveaux_noms)]

# 3. supression des doublons de l'arbre
tree <- drop.tip(tree, tips_a_supprimer)

# 4. renome les noms
tree$tip.label <- sapply(strsplit(tree$tip.label, "_"), function(x) paste(x[3], x[4], sep = " "))

# 5. Vérification
sum(duplicated(tree$tip.label))  # doit être 0
length(tree$tip.label)           # nombre de tips restants
head(tree$tip.label)

# Remettre les rownames
BDD_moy_esp_phylo <- BDD_moy_score
rownames(BDD_moy_esp_phylo) <- BDD_moy_esp_phylo$Nom_scientifique

# Vérifier la concordance arbre <-> BDD
check <- name.check(tree, BDD_moy_esp_phylo)

length(check$data_not_tree)  # espèces dans BDD mais pas dans l'arbre

# Trouver lignes dont  Nom_scientifique est dans check$data_not_tree
lignes <- c()
for (i in 1:nrow(BDD_moy_esp_phylo)) {
  if (length(which(check$data_not_tree == BDD_moy_esp_phylo$Nom_scientifique[i])) > 0) {
    lignes <- c(lignes, i)
  }
}

# Extraire les genres de ces lignes
genres_manquants <- unique(BDD_moy_esp_phylo$Genre[lignes])
genres_manquants

# Créer genres_arbre
genres_arbre <- unique(sapply(strsplit(tree$tip.label, " "), function(x) x[1]))

# Vérifier
head(genres_arbre)

# Genres qui ont au moins un représentant dans l'arbre
presents <- c()
for (g in genres_manquants) {
  if (length(which(genres_arbre == g)) > 0) {
    presents <- c(presents, g)
  }
}
presents

# Genres totalement absents de l'arbre
absents <- c()
for (g in genres_manquants) {
  if (length(which(genres_arbre == g)) == 0) {
    absents <- c(absents, g)
  }
}
absents


BDD_moy_esp_phylo <- BDD_moy_score[BDD_moy_score$Genre != "Agathis" & BDD_moy_score$Genre != "Archidendropsis", ]
nrow(BDD_moy_score)           # doit toujours être 67
nrow(BDD_moy_esp_phylo)  # doit être 65

rownames(BDD_moy_esp_phylo) <- BDD_moy_esp_phylo$Nom_scientifique
check <- name.check(tree, BDD_moy_esp_phylo)
length(check$data_not_tree)  # doit être 49



missing_species <- check$data_not_tree

for(sp in missing_species){
  
  genus <- BDD_moy_esp_phylo$Genre[BDD_moy_esp_phylo$Nom_scientifique == sp]
  
  if(length(genus) == 1 && !is.na(genus)){
    
    genus_tips <- tree$tip.label[grepl(paste0("^", genus, " "), tree$tip.label)]
    
    # CAS 1 : plusieurs espèces du genre → on greffe au MRCA
    if(length(genus_tips) >= 2){
      node <- getMRCA(tree, genus_tips)
      tree <- bind.tip(tree, sp, where = node, position = 0)
    }
    
    # CAS 2 : une seule espèce du genre → on greffe dessus
    else if(length(genus_tips) == 1){
      ref <- genus_tips[1]
      tree <- bind.tip(tree, tip.label = sp, where = which(tree$tip.label == ref), position = 0)
      message(paste("Attaché à espèce unique du genre :", genus))
    }
    
    # CAS 3 : genre absent (ne devrait plus arriver)
    else {
      message(paste("Genre absent de l'arbre :", genus))
    }
  }
}

# Vérifier
check2 <- name.check(tree, BDD_moy_esp_phylo)
length(check2$data_not_tree)  # doit être 0

# Résoudre les polytomies créées par le greffage
tree <- multi2di(tree)
is.binary(tree)  # doit être TRUE

# Élaguer l'arbre pour ne garder que tes 65 espèces
tree_final <- keep.tip(tree, BDD_moy_esp_phylo$Nom_scientifique)
length(tree_final$tip.label)  # doit être 65


#Suppression des 102 doublons dans l'arbre source
#Renommage des tip labels en Genre espece
#Exclusion des 2 espèces sans genre dans l'arbre (Agathis, Archidendropsis)
#Greffage des 49 espèces manquantes
#Résolution des polytomies
#Élagage final à 65 espèces



#################SIGNAL PHYLO###################################
library(phytools)
tree_final <- compute.brlen(tree_final, method = "Grafen") #### pour redonner longueur au branches de l'arbre
is.ultrametric(tree_final)  # doit être TRUE

variables_toutes <- c("DI_test", "BB_test", "BT", "MT", "FI","score", 
                      "BD_mean", "Nb_rami", "TD", "Gmin", "LDMC", "SLA", "TDMC")

resultats_tous <- data.frame(
  variable = variables_toutes,
  lambda = NA,
  lambda_p = NA,
  K = NA,
  K_p = NA
)

for(i in 1:length(variables_toutes)){
  
  trait <- setNames(BDD_moy_esp_phylo[[variables_toutes[i]]], BDD_moy_esp_phylo$Nom_scientifique)
  
  res_lambda <- phylosig(tree_final, trait, method = "lambda", test = TRUE)
  resultats_tous$lambda[i] <- res_lambda$lambda
  resultats_tous$lambda_p[i] <- res_lambda$P
  
  res_K <- phylosig(tree_final, trait, method = "K", test = TRUE)
  resultats_tous$K[i] <- res_K$K
  resultats_tous$K_p[i] <- res_K$P
}

resultats_tous





################################################################################
############################ ECHELLE ESPECE ####################################
################################################################################



################## DISTRIBUTION DES DONNEES #####################
#Infla
hist(BDD_moy_esp3$DI_test,xlab="DI",main="DI distribution",xlim=c(0,10),breaks=seq(0,10,1)) # asymétrique droite
hist(BDD_moy_esp3$BT_test,xlab="BT",main="BT distribution",breaks=seq(0,120,10)) # asymétrique droite
hist(BDD_moy_esp3$BB_test,,xlab="BB",main="BB distribution") # normale proportion
hist(BDD_moy_espMT$MT,xlab="MT",main="MT distribution",breaks=seq(0,1000,100))      # normale 

#Infla ech
hist(BDD_moy_ech3$DI_test,xlab="DI",main="DI distribution",xlim=c(0,10),breaks=seq(0,10,0.5)) # asymétrique droite
hist(BDD_moy_ech3$BT_test,xlab="BT",main="BT distribution") # asymétrique droite
hist(BDD_moy_ech3$BB_test,,xlab="BB",main="BB distribution",breaks=seq(0,100,10)) # normale proportion
hist(BDD_moy_echMT$MT,xlab="MT",main="MT distribution") 


#traits
hist(BDD_moy_esp$Nb_rami)  
hist(BDD_moy_esp$BD_mean)       
hist(BDD_moy_esp$TMC_t0)         
hist(BDD_moy_esp$TMC_t24)  
hist(BDD_moy_esp$TDMC)         
hist(BDD_moy_esp$TD)            
hist(BDD_moy_esp$TDIA)           
hist(BDD_moy_esp$LMC_t0)         
hist(BDD_moy_esp$LMC_t24)  
hist(BDD_moy_esp$LDMC)          
hist(BDD_moy_esp$Surface_F) 
hist(BDD_moy_esp$SLA)            
hist(BDD_moy_esp$LT)        


#standardisation des données 
BDD_moy_esp$BD_mean_cr<-as.numeric(scale(BDD_moy_esp$BD_mean))
BDD_moy_esp$TD_cr<-as.numeric(scale(BDD_moy_esp$TD))
BDD_moy_esp$LA_cr<-as.numeric(scale(BDD_moy_esp$Surface_F))
BDD_moy_esp$LDMC_cr<-as.numeric(scale(BDD_moy_esp$LDMC))
BDD_moy_esp$LMC_t24_cr<-as.numeric(scale(BDD_moy_esp$LMC_t24))
BDD_moy_esp$SLA_cr<-as.numeric(scale(BDD_moy_esp$SLA))
BDD_moy_esp$Nb_rami_cr<-as.numeric(scale(BDD_moy_esp$Nb_rami))
BDD_moy_esp$Gmin_cr<-as.numeric(scale(BDD_moy_esp$Gmin))






############### MODELES classiques sans pgls #############################
##########################################################################



###### FI ########

#modèle
mFI3<-glm(cbind(Nb_FI,Nb_essais-Nb_FI)~Gmin_cr+BD_mean_cr+TD_cr+SLA_cr+Nb_rami_cr+LDMC_cr,data=BDD_moy_esp,family="binomial")
summary(mFI3)

mFI4<-glm(cbind(Nb_FI,Nb_essais-Nb_FI)~Gmin_cr+BD_mean_cr+TD_cr+SLA_cr+Nb_rami_cr+LMC_t24_cr,data=BDD_moy_esp,family="binomial")
summary(mFI4)

AIC(mFI3,mFI4)

###################### graph bleu rouge avec LDMC+LMC_t24 ###############

# Extraire coefficients (sans l'intercept)
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


#####  Plot comparaison modèle 1 et 2 ####
par(mar=c(4,6,2,2))
plot(estimates1, y_pos1,type = "n",
     xlim = c(-4,4),
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
axis(2, at = c(5.8,4.8,3.8,2.8,1.8,1,0.6), labels = c("Gmin","BD","TD","SLA","Nb_Rami","LDMC","LMC_t24"), las = 1,cex.axis = 0.9)

# Axe X
axis(1,cex.axis = 0.9)

# Valeurs des coefficients + étoiles
text(estimates1, y_pos1 + 0.16,
     labels = paste0(round(estimates1, 2), stars1),
     col = cols1, font = 2, cex =  0.9)
text(estimates2, y_pos2 -0.24,
     labels = paste0(round(estimates2, 2), stars2),
     col = cols2, font = 2, cex =  0.9)



#####  Plot seulement modèle 1 ####
plot(estimates1, y_pos1,type = "n",
     xlim = c(-4,4),
     ylim = c(0.7,length(estimates1) +0.2 ),
     xlab = "Estimate",
     ylab = "",
     axes = FALSE,
     main = "Fréquence d'ignition",
     cex.main = 1.1)

# Ligne verticale à zéro
abline(v = 0, lty = 2)
abline(h = y_pos1, lwd = 0.5, lty = 3, col="grey")

#points des coefs
points(estimates1, y_pos1 ,pch = 16, col = cols1,cex = 1.1)

# Barres d'erreur (IC)
segments(ci1[,1], y_pos1, ci1[,2], y_pos1, col = cols1)

# Axe Y avec noms des variables
axis(2, at = c(6,5,4,3,2,1), labels = c("Gmin","BD","TD","SLA","Nb_Rami","LDMC"), las = 1,cex.axis = 0.9)

# Axe X
axis(1,cex.axis = 0.9)

# Valeurs des coefficients + étoiles
text(estimates1, y_pos1 + 0.16,
     labels = paste0(round(estimates1, 2), stars1),
     col = cols1, font = 2, cex =  0.9)








###################### prédicion ############### mettre LMC_t24 ou LDMC
# Modèle GLM binomial (mFI3 ou mFI4)
moy <- mean(BDD_moy_esp$LMC_t24, na.rm = TRUE)
moy
ecart <- sd(BDD_moy_esp$LMC_t24, na.rm = TRUE)
ecart

# Valeurs de LMC_t24
foo <- seq(min(BDD_moy_esp$LMC_t24_cr), max(BDD_moy_esp$LMC_t24_cr),length.out = 100)

# on fait varier LMDC et on fixe les autres variables
pred_data1 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_moy_esp$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_esp$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_esp$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_esp$SLA_cr), length(foo))
)


# Prédictions et IC
pred1 <- predict(mFI4, type = "response", newdata = pred_data1, se.fit = TRUE)
pred_data1$LMC_t24_cr[pred1$fit<0.5]*ecart+moy
pred_data1$LMC_t24_cr[(pred1$fit+(1.96*pred1$se.fit))<0.5]*ecart+moy
pred_data1$LMC_t24_cr[(pred1$fit-(1.96*pred1$se.fit))<0.5]*ecart+moy

# Plot des points observés
plot(BDD_moy_esp$LMC_t24, BDD_moy_esp$Nb_FI / BDD_moy_esp$Nb_essais, 
     xlab = "LMC_t24 (%)", ylab = "Fréquence d'ignition",ylim = c(0,1), 
     main = "Effet de LMC_t24 sur la fréquence d'ignition")

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 
abline(v = 335,  lty = 2, col = "blue")
abline(v = 397,  lty = 2,col="red")
abline(v = 479,  lty = 2, col = "blue")
segments(335, 0.5, 479,0.5, col = "blue",lwd = 2)
points(397,0.5,pch=19,col="red")

axis(2)

# Intervalle de confiance pour SD moyen
# Calcul borne supérieure (limitée à 1)
lines((foo * ecart + moy), pmin(pred1$fit + 1.96 * pred1$se.fit, 1), col = "blue", lty = 3)

# Calcul borne inférieure (limitée à 0)
lines((foo * ecart + moy), pmax(pred1$fit - 1.96 * pred1$se.fit, 0), col = "blue", lty = 3)











#################################################################
################### ECHELLE ECHANTILLON #########################
#################################################################


################## DISTRIBUTION DES DONNEES #####################

#Infla ech
hist(BDD_moy_ech3$DI_test,xlab="DI",main="Distribution DI",xlim=c(0,7),breaks=seq(0,9,0.5)) # asymétrique droite
hist(BDD_moy_ech3$BT_test,xlab="BT",main="Distribution BT") # asymétrique droite
hist(BDD_moy_ech3$BB_test,,xlab="BB",main="Distribution BB",breaks=seq(0,100,10)) # normale proportion
hist(BDD_moy_echMT$MT,xlab="MT",main="Distribution MT") 
hist(BDD_moy_ech$score,xlab="Score d'inflammabilité",main="Distribution score")

#traits
hist(BDD_moy_ech$Nb_rami)  
hist(BDD_moy_ech$SD)       
hist(BDD_moy_ech$TMC_t0)         
hist(BDD_moy_ech$TMC_t24)  
hist(BDD_moy_ech$TDMC)         
hist(BDD_moy_ech$TD)            
hist(BDD_moy_ech$TDIA)           
hist(BDD_moy_ech$LMC_t0)         
hist(BDD_moy_ech$LMC_t24)  
hist(BDD_moy_ech$LDMC)          
hist(BDD_moy_ech$Surface_F) 
hist(BDD_moy_ech$SLA)            
hist(BDD_moy_ech$LT)        



#standardisation des données 
BDD_moy_ech$BD_mean_cr<-as.numeric(scale(BDD_moy_ech$BD_mean))
BDD_moy_ech$TD_cr<-as.numeric(scale(BDD_moy_ech$TD))
BDD_moy_ech$LA_cr<-as.numeric(scale(BDD_moy_ech$Surface_F))
BDD_moy_ech$LDMC_cr<-as.numeric(scale(BDD_moy_ech$LDMC))
BDD_moy_ech$LMC_t24_cr<-as.numeric(scale(BDD_moy_ech$LMC_t24))
BDD_moy_ech$Gmin_cr<-as.numeric(scale(BDD_moy_ech$Gmin))
BDD_moy_ech$SLA_cr<-as.numeric(scale(BDD_moy_ech$SLA))
BDD_moy_ech$Nb_rami_cr<-as.numeric(scale(BDD_moy_ech$Nb_rami))


BDD_moy_echMT$BD_mean_cr<-as.numeric(scale(BDD_moy_echMT$BD_mean))
BDD_moy_echMT$TD_cr<-as.numeric(scale(BDD_moy_echMT$TD))
BDD_moy_echMT$LA_cr<-as.numeric(scale(BDD_moy_echMT$Surface_F))
BDD_moy_echMT$LDMC_cr<-as.numeric(scale(BDD_moy_echMT$LDMC))
BDD_moy_echMT$LMC_t24_cr<-as.numeric(scale(BDD_moy_echMT$LMC_t24))
BDD_moy_echMT$Gmin_cr<-as.numeric(scale(BDD_moy_echMT$Gmin))
BDD_moy_echMT$SLA_cr<-as.numeric(scale(BDD_moy_echMT$SLA))
BDD_moy_echMT$Nb_rami_cr<-as.numeric(scale(BDD_moy_echMT$Nb_rami))

BDD_moy_ech3$BD_mean_cr<-as.numeric(scale(BDD_moy_ech3$BD_mean))
BDD_moy_ech3$TD_cr<-as.numeric(scale(BDD_moy_ech3$TD))
BDD_moy_ech3$LA_cr<-as.numeric(scale(BDD_moy_ech3$Surface_F))
BDD_moy_ech3$LDMC_cr<-as.numeric(scale(BDD_moy_ech3$LDMC))
BDD_moy_ech3$LMC_t24_cr<-as.numeric(scale(BDD_moy_ech3$LMC_t24))
BDD_moy_ech3$Gmin_cr<-as.numeric(scale(BDD_moy_ech3$Gmin))
BDD_moy_ech3$SLA_cr<-as.numeric(scale(BDD_moy_ech3$SLA))
BDD_moy_ech3$Nb_rami_cr<-as.numeric(scale(BDD_moy_ech3$Nb_rami))

plot(BDD_moy_ech$LDMC, BDD_moy_ech$LMC_t24)
m<-glm(LMC_t24~LDMC,data=BDD_moy_ech,family=Gamma(link="log"))
summary(m)


foo <- seq(min(BDD_moy_ech$LDMC), max(BDD_moy_ech$LDMC), length.out = 100)
pred_data1 <- data.frame(LDMC = foo)
pred1 <- predict(m, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
plot(BDD_moy_ech$LDMC,BDD_moy_ech$LMC_t24,main="Relation entre LDMC et LMC_t24",xlab="LDMC (mg/g)",ylab="LMC_t24 (%)")
lines(foo, pred1$fit, col = "black", lwd = 2) 
# Intervalle de confiance pour SD moyen
lines(foo, pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines(foo, pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)




#######################################################
###### MT #########
#######################################################

#modèles
m_MT1 <- lmer(MT ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LDMC_cr + (1 | Nom_scientifique), data = BDD_moy_echMT)
summary(m_MT1)
r.squaredGLMM(m_MT1)
m_MT3 <- lmer(MT ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LMC_t24_cr + (1 | Nom_scientifique), data = BDD_moy_echMT)
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
axis(2, at = c(4.8,3.8,2.8,1.8,1,0.6), labels = c("SD","TD","SLA","Nb_rami","LDMC","LMC_t24"), las = 1,cex.axis = 0.9)


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
moy <- mean(BDD_moy_echMT$LDMC, na.rm = TRUE)
ecart <- sd(BDD_moy_echMT$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_moy_echMT$LDMC_cr), max(BDD_moy_echMT$LDMC_cr),length.out = 100)

############## on fait varier LMDC
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_echMT$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_echMT$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_echMT$SLA_cr), length(foo))
)


# Prédictions
pred1 <- predict(m_MT1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
par(mar=c(4,5,4,4))
plot(BDD_moy_echMT$LDMC, BDD_moy_echMT$MT,type="n", 
     xlab = "LDMC (mg/g)", ylab = "Température maximum (°C)", 
     main = "Effet de LDMC sur MT",ylim=c(500,850),xlim=c(175,550))

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 

# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)



############  LMC_t24
moy <- mean(BDD_moy_echMT$LMC_t24, na.rm = TRUE)
ecart <- sd(BDD_moy_echMT$LMC_t24, na.rm = TRUE)

# Valeurs de LMC_t24
foo <- seq(min(BDD_moy_echMT$LMC_t24_cr), max(BDD_moy_echMT$LMC_t24_cr),length.out = 100)

############## on fait varier LMDC
pred_data1 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_moy_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_echMT$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_echMT$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_echMT$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_MT3, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_moy_echMT$LMC_t24, BDD_moy_echMT$MT,type="n", 
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
  SD_cr = rep(mean(BDD_moy_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_echMT$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_echMT$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_echMT$SLA_cr), length(foo))
)

pred_data2 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(quantile(BDD_moy_echMT$TD_cr,0.75), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_echMT$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_echMT$SLA_cr), length(foo))
)

pred_data3 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(quantile(BDD_moy_echMT$TD_cr,0.25), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_echMT$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_echMT$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_MT1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred2 <- predict(m_MT1, type = "response", newdata = pred_data2, se.fit = TRUE, re.form = NA)
pred3 <- predict(m_MT1, type = "response", newdata = pred_data3, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_moy_echMT$LDMC, BDD_moy_echMT$MT,type="n", 
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
BDD_moy_ech3$BB_prop <- BDD_moy_ech3$BB_test/100


m_BB0 <- lmer(BB_prop ~ SD_cr+TD_cr+Nb_rami_cr+LDMC_cr+LT_cr+VPD_cr + (1 | Nom_scientifique), data = BDD_moy_ech3)    #### meilleur modèle 
summary(m_BB0)
r.squaredGLMM(m_BB0)
m_BB1 <- lmer(BB_prop ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LDMC_cr + (1 | Nom_scientifique), data = BDD_moy_ech3)    #### meilleur modèle 
summary(m_BB1)
r.squaredGLMM(m_BB1)
m_BB2 <- lmer(BB_prop ~ SD_cr+TD_cr+Nb_rami_cr+LMC_t24_cr+LT_cr+VPD_cr + (1 | Nom_scientifique), data = BDD_moy_ech3)    #### meilleur modèle 
summary(m_BB2)
r.squaredGLMM(m_BB2)
m_BB3 <- lmer(BB_prop ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LMC_t24_cr + (1 | Nom_scientifique), data = BDD_moy_ech3)    #### meilleur modèle 
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
axis(2, at = c(4.8,3.8,2.8,1.8,1,0.6), labels = c("SD","TD","SLA","Nb_rami","LDMC","LMC_t24"), las = 1,cex.axis = 0.9)

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
moy <- mean(BDD_moy_ech3$LDMC, na.rm = TRUE)
ecart <- sd(BDD_moy_ech3$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_moy_ech3$LDMC_cr), max(BDD_moy_ech3$LDMC_cr),length.out = 100)

# on fait varier LMDC 
par(mar=c(4,5,4,4))

pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_ech3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_ech3$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_ech3$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_ech3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_BB1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_moy_ech3$LDMC, BDD_moy_ech3$BB_prop,type="n", 
     xlab = "LDMC", ylab = "Biomasse brûlée (%)", 
     main = "Effet de LDMC sur BB ",ylim=c(0,1),xlim=c(175,550))

# Courbes de prédiction
lines((foo*ecart + moy), pmin(pmax(pred1$fit, 0), 1), col = "black", lwd = 2 )


# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)


########## LMC_t24
moy <- mean(BDD_moy_ech3$LMC_t24, na.rm = TRUE)
ecart <- sd(BDD_moy_ech3$LMC_t24, na.rm = TRUE)

# Valeurs de LMC_t24
foo <- seq(min(BDD_moy_ech3$LMC_t24_cr), max(BDD_moy_ech3$LMC_t24_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_moy_ech3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_ech3$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_ech3$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_ech3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_BB3, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_moy_ech3$LMC_t24, BDD_moy_ech3$BB_prop,type="n", 
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
  SD_cr = rep(mean(BDD_moy_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_echMT$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_moy_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_echMT$SLA_cr), length(foo))
)

pred_data2 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(quantile(BDD_moy_echMT$TD_cr,0.75), length(foo)),
  LA_cr = rep(mean(BDD_moy_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_echMT$SLA_cr), length(foo))
)

pred_data3 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(quantile(BDD_moy_echMT$TD_cr,0.25), length(foo)),
  LA_cr = rep(mean(BDD_moy_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_echMT$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_BB1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred2 <- predict(m_BB1, type = "response", newdata = pred_data2, se.fit = TRUE, re.form = NA)
pred3 <- predict(m_BB1, type = "response", newdata = pred_data3, se.fit = TRUE, re.form = NA)


# Plot des points observés
plot(BDD_moy_ech3$LDMC, BDD_moy_ech3$BB_prop,type="n", 
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
m_BT0 <- glmer(BT_test ~ SD_cr+TD_cr+Nb_rami_cr+LDMC_cr+LT_cr+VPD_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_moy_ech3)  #meilleur modèle
summary(m_BT0)
r.squaredGLMM(m_BT0)
m_BT1 <- glmer(BT_test ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LDMC_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_moy_ech3)  #meilleur modèle
summary(m_BT1)
r.squaredGLMM(m_BT1)
m_BT2 <- glmer(BT_test ~ SD_cr+TD_cr+Nb_rami_cr+LMC_t24_cr+LT_cr+VPD_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_moy_ech3)  #meilleur modèle
summary(m_BT2)
r.squaredGLMM(m_BT2)
m_BT3 <- glmer(BT_test ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LMC_t24_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_moy_ech3)  #meilleur modèle
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
axis(2, at = c(4.8,3.8,2.8,1.8,1,0.6), labels = c("SD","TD","SLA","Nb_rami","LDMC","LMC_t24"), las = 1,cex.axis = 0.9)

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
moy <- mean(BDD_moy_ech3$LDMC, na.rm = TRUE)
ecart <- sd(BDD_moy_ech3$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_moy_ech3$LDMC_cr), max(BDD_moy_ech3$LDMC_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_ech3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_ech3$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_ech3$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_ech3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_BT1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred1

# Plot des points observés
plot(BDD_moy_ech3$LDMC, BDD_moy_ech3$BT_test,type="n", 
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
moy <- mean(BDD_moy_ech3$LMC_t24, na.rm = TRUE)
ecart <- sd(BDD_moy_ech3$LMC_t24, na.rm = TRUE)

# Valeurs de LMC_t24
foo <- seq(min(BDD_moy_ech3$LMC_t24_cr), max(BDD_moy_ech3$LMC_t24_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_moy_ech3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_ech3$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_ech3$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_ech3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_BT3, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_moy_ech3$LMC_t24, BDD_moy_ech3$BT,type="n", 
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
  SD_cr = rep(mean(BDD_moy_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_echMT$TD_cr), length(foo)),
  LA_cr = rep(mean(BDD_moy_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_echMT$SLA_cr), length(foo))
)

pred_data2 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(quantile(BDD_moy_echMT$TD_cr,0.75), length(foo)),
  LA_cr = rep(mean(BDD_moy_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_echMT$SLA_cr), length(foo))
)

pred_data3 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_echMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(quantile(BDD_moy_echMT$TD_cr,0.25), length(foo)),
  LA_cr = rep(mean(BDD_moy_echMT$LA_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_echMT$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_BT1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred2 <- predict(m_BT1, type = "response", newdata = pred_data2, se.fit = TRUE, re.form = NA)
pred3 <- predict(m_BT1, type = "response", newdata = pred_data3, se.fit = TRUE, re.form = NA)


# Plot des points observés
plot(BDD_moy_ech3$LDMC, BDD_moy_ech3$BT,type="n",ylim=c(0,80), 
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

m_DI0 <- glmer(DI_test ~ SD_cr+TD_cr+Nb_rami_cr+LDMC_cr+SLA_cr+VPD_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_moy_ech3)  #meilleur modèle
summary(m_DI0)
r.squaredGLMM(m_DI0)
m_DI1 <- glmer(DI_test ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LDMC_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_moy_ech3)  #meilleur modèle
summary(m_DI1)
r.squaredGLMM(m_DI1)
m_DI2 <- glmer(DI_test ~ SD_cr+TD_cr+Nb_rami_cr+LMC_t24_cr+LT_cr+VPD_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_moy_ech3)  #meilleur modèle
summary(m_DI2)
r.squaredGLMM(m_DI2)
m_DI3 <- glmer(DI_test ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LMC_t24_cr + (1 | Nom_scientifique), family=Gamma(link="log"), data = BDD_moy_ech3)  #meilleur modèle
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
axis(2, at = c(4.8,3.8,2.8,1.8,1,0.6), labels = c("SD","TD","SLA","Nb_rami","LDMC","LMC_t24"), las = 1,cex.axis = 0.9)

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
moy <- mean(BDD_moy_ech3$LDMC, na.rm = TRUE)
ecart <- sd(BDD_moy_ech3$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_moy_ech3$LDMC_cr), max(BDD_moy_ech3$LDMC_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_ech3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_ech3$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_ech3$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_ech3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_DI1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_moy_ech3$LDMC, BDD_moy_ech3$DI_test,type="n", 
     xlab = "LDMC", ylab = "Délai d'ignition (s)", 
     main = "Effet de LDMC sur DI ",ylim=c(0.5,2),xlim=c(175,550))

# Courbes de prédiction
lines((foo*ecart + moy), pred1$fit, col = "black", lwd = 2 )


# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)


########## LMC_t24
moy <- mean(BDD_moy_ech3$LMC_t24, na.rm = TRUE)
ecart <- sd(BDD_moy_ech3$LMC_t24, na.rm = TRUE)

# Valeurs de LMC_t24
foo <- seq(min(BDD_moy_ech3$LMC_t24_cr), max(BDD_moy_ech3$LMC_t24_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_moy_ech3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_ech3$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_ech3$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_ech3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_DI3, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_moy_ech3$LMC_t24, BDD_moy_ech3$DI_test,type="n", 
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
  SD_cr = rep(mean(BDD_moy_ech3$SD_cr,na.rm = TRUE), length(foo)),
  SLA_cr = rep(mean(BDD_moy_ech3$SLA_cr), length(foo)),
  LA_cr = rep(mean(BDD_moy_ech3$LA_cr), length(foo)),
  TD_cr = rep(mean(BDD_moy_ech3$TD_cr), length(foo))
)


pred_data2 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_ech3$SD_cr,na.rm = TRUE), length(foo)),
  SLA_cr = rep(quantile(BDD_moy_ech3$SLA_cr,0.75), length(foo)),
  LA_cr = rep(mean(BDD_moy_ech3$LA_cr), length(foo)),
  TD_cr = rep(mean(BDD_moy_ech3$TD_cr), length(foo))
)

pred_data3 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_ech3$SD_cr,na.rm = TRUE), length(foo)),
  SLA_cr = rep(quantile(BDD_moy_ech3$SLA_cr,0.25), length(foo)),
  LA_cr = rep(mean(BDD_moy_ech3$LA_cr), length(foo)),
  TD_cr = rep(mean(BDD_moy_ech3$TD_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_DI1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred2 <- predict(m_DI1, type = "response", newdata = pred_data2, se.fit = TRUE, re.form = NA)
pred3 <- predict(m_DI1, type = "response", newdata = pred_data3, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_moy_ech3$LDMC, BDD_moy_ech3$DI_test,type="n", 
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
m_score1 <- lmer(score ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LDMC_cr + (1 | Nom_scientifique), data = BDD_moy_ech)
summary(m_score1)
r.squaredGLMM(m_score1)
m_score3 <- lmer(score ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LMC_t24_cr + (1 | Nom_scientifique), data = BDD_moy_ech)
summary(m_score3)
m_score2 <- lmer(score ~ LDMC_cr + (1 | Nom_scientifique), data = BDD_moy_ech)
summary(m_score2)
r.squaredGLMM(m_score2)

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
axis(2, at = c(4.8,3.8,2.8,1.8,1,0.6), labels = c("SD","TD","SLA","Nb_rami","LDMC","LMC_t24"), las = 1,cex.axis = 0.9)

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
moy <- mean(BDD_moy_ech$LDMC, na.rm = TRUE)
ecart <- sd(BDD_moy_ech$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_moy_ech$LDMC_cr), max(BDD_moy_ech$LDMC_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_ech$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_ech$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_ech$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_ech$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_score1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

#background transparent 
par(bg = NA)

# Plot des points observés
plot(BDD_moy_ech$LDMC, BDD_moy_ech$score,type="n",cex.axis=1.2, 
     xlab = "LDMC", ylab = "Flammability score", 
     main = "Effet de LDMC sur score ",ylim=c(-4,4),xlim=c(150,600))
box(lwd = 2)
points(BDD_moy_ech$LDMC, BDD_moy_ech$score, pch=16, col="#634124")

col_ic <- rgb(109/255, 168/255, 111/255)

# Polygone de l'intervalle de confiance
x_vals <- (foo * ecart + moy)

polygon(
  x = c(x_vals, rev(x_vals)),
  y = c(pred1$fit + 1.96 * pred1$se.fit,
        rev(pred1$fit - 1.96 * pred1$se.fit)),
  col = col_ic,
  border = NA
)




# Courbes de prédiction
lines((foo*ecart + moy), pred1$fit, col = "#C94802", lwd = 3 )
lines((foo*ecart + moy), pred1$fit,  lwd = 3 )

# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)





########## LMC_t24
moy <- mean(BDD_moy_ech$LMC_t24, na.rm = TRUE)
ecart <- sd(BDD_moy_ech$LMC_t24, na.rm = TRUE)

# Valeurs de LMC_t24
foo <- seq(min(BDD_moy_ech$LMC_t24_cr), 5,length.out = 100)
max(BDD_moy_ech$LMC_t24_cr)
# on fait varier LMDC 
pred_data1 <- data.frame(
  LMC_t24_cr = foo,
  SD_cr = rep(mean(BDD_moy_ech$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_ech$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_ech$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_ech$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_score3, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_moy_ech$LMC_t24, BDD_moy_ech$score,type="n", 
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

hist(BDD_moy_esp$LMC_t24)
hist(BDD_moy_esp$LMC_t0)
hist(BDD_moy_esp$Gmin)
BDD_moy_esp$ratio <- BDD_moy_esp$LMC_t24 / BDD_moy_esp$LMC_t0
hist(BDD_moy_esp$ratio)
BDD_moy_esp$logLMC24<-log(BDD_moy_esp$LMC_t24)
BDD_moy_esp$logLMC0<-log(BDD_moy_esp$LMC_t0)

BDD_moy_esp$LMC_t0_cr<-as.numeric(scale(BDD_moy_esp$LMC_t0))
BDD_moy_esp$Gmin_cr<-as.numeric(scale(BDD_moy_esp$Gmin))


m_gmin2<-glm(logLMC24~logLMC0+Gmin_cr,family=gaussian,data=BDD_moy_esp)
summary(m_gmin2)
m_gmin2<-glm(LMC_t24~LMC_t0_cr+Gmin_cr,family=Gamma(link = "log"),data=BDD_moy_esp)
summary(m_gmin2)
m_gmin3<-glm(ratio~LMC_t0_cr+Gmin_cr,family=gaussian,data=BDD_moy_esp)
summary(m_gmin3)

AIC(m_gmin2,m_gmin3)

hist(resid(m_gmin3)) ### gmin3 retenu!!!
plot(m_gmin3)
hist

################# Prédiction ########################
# LMC_t0
moy <- mean(BDD_moy_esp$LMC_t0, na.rm = TRUE)
ecart <- sd(BDD_moy_esp$LMC_t0, na.rm = TRUE)

# Valeurs de LMC_t0
foo <- seq(min(BDD_moy_esp$LMC_t0_cr), max(BDD_moy_esp$LMC_t0_cr),length.out = 100)

# on fait varier LMDC et on fixe les autres variables
pred_data1 <- data.frame(
  LMC_t0_cr = foo,
  Gmin_cr = rep(mean(BDD_moy_esp$Gmin_cr,na.rm = TRUE), length(foo))
)

pred_data2 <- data.frame(
  LMC_t0_cr = foo,
  Gmin_cr = rep(min(BDD_moy_esp$Gmin_cr,na.rm = TRUE), length(foo))
)

pred_data3 <- data.frame(
  LMC_t0_cr = foo,
  Gmin_cr = rep(max(BDD_moy_esp$Gmin_cr,na.rm = TRUE), length(foo))
)




# Prédictions avec gmin3
pred1 <- predict(m_gmin3, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred2 <- predict(m_gmin3, type = "response", newdata = pred_data2, se.fit = TRUE, re.form = NA)
pred3 <- predict(m_gmin3, type = "response", newdata = pred_data3, se.fit = TRUE, re.form = NA)


# Plot des points observés
plot(BDD_moy_esp$LMC_t0, BDD_moy_esp$LMC_t24, 
     xlab = "LMC_t0 (%)", ylab = "LMC_t24 (%)", 
     main = "Effet de LMC_t0 sur LMC_t24 selon gmin ",ylim=c(0,700))

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit*(foo*ecart+moy), col = "#3B393B", lwd = 2) 
lines((foo*ecart+moy), pred2$fit*(foo*ecart+moy), col = "#005ED1" ,lwd = 2) 
lines((foo*ecart+moy), pred3$fit*(foo*ecart+moy), col = "#BD0000",lwd = 2)


polygon(x = c(-50, -50, 850), 
        y = c(850, -50, 850), 
        density = 20, angle = -45, col = "#615B5B", border = NA)

# Ajouter la ligne y = x
abline(a = 0, b = 1)

esp_a_annoter <- BDD_moy_esp$Nom_scientifique == "Scaevola taccada" | 
  BDD_moy_esp$Nom_scientifique == "Scaevola montana" |
  BDD_moy_esp$Nom_scientifique== "Barringtonia asiatica" |
  BDD_moy_esp$Nom_scientifique== "Hibiscus tiliaceus R" 

points (BDD_moy_esp$LMC_t0[esp_a_annoter],
        BDD_moy_esp$LMC_t24[esp_a_annoter],pch=21,bg = rgb(1,0,0),cex=1.2
)

text(BDD_moy_esp$LMC_t0[esp_a_annoter] -10 ,
     BDD_moy_esp$LMC_t24[esp_a_annoter] + 25,
     labels = c("B. asiatica","H. tiliaceus R","S. montana","S. taccada"),
     cex = 1)







# Prédictions avec gmin2
pred1 <- predict(m_gmin2, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred2 <- predict(m_gmin2, type = "response", newdata = pred_data2, se.fit = TRUE, re.form = NA)
pred3 <- predict(m_gmin2, type = "response", newdata = pred_data3, se.fit = TRUE, re.form = NA)


# Créer le graphique
par(mar = c(5,5,5,5))
plot(BDD_moy_esp$LMC_t0, BDD_moy_esp$LMC_t24, 
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

esp_a_annoter <- BDD_moy_esp$Nom_scientifique == "Scaevola taccada" | 
  BDD_moy_esp$Nom_scientifique == "Scaevola montana" |
  BDD_moy_esp$Nom_scientifique== "Barringtonia asiatica"  |
  BDD_moy_esp$Nom_scientifique== "Morinda citrifolia" |
  BDD_moy_esp$Nom_scientifique== "Acalypha sp" 

points (BDD_moy_esp$LMC_t0[esp_a_annoter],
        BDD_moy_esp$LMC_t24[esp_a_annoter],pch=21,bg = rgb(1,0,0),cex=1.2
)

text(BDD_moy_esp$LMC_t0[esp_a_annoter] -10 ,
     BDD_moy_esp$LMC_t24[esp_a_annoter] + 25,
     labels = c("B. asiatica","S. montana","S. taccada","M. citrifolia", "Acalypha"),
     cex = 1)









# ZOOM
plot(BDD_moy_esp$LMC_t0, BDD_moy_esp$LMC_t24, 
     xlab = "LMC_t0 (%)", ylab = "LMC_t24 (%)", 
     main = "Effet de LMC_t0 sur LMC_t24 selon gmin",
     xlim = c(70,410), ylim = c(0,400))


# Ajouter les courbes de prédiction
lines((foo*ecart+moy), pred1$fit*(foo*ecart+moy), col = "#3B393B", lwd = 2) 
lines((foo*ecart+moy), pred2$fit*(foo*ecart+moy), col = "#005ED1" ,lwd = 2) 
lines((foo*ecart+moy), pred3$fit*(foo*ecart+moy), col = "#BD0000",lwd = 2)


# Ajouter le polygone hachuré (zone où y > x)
polygon(x = c(-50, -50, 850), 
        y = c(850, -50, 850), 
        density = 15, angle = -45, col = "#615B5B", border = NA)

# Ajouter la ligne y = x
abline(a = 0, b = 1)

esp_a_annoter <- BDD_moy_esp$Nom_scientifique == "Scaevola montana" | 
  BDD_moy_esp$Nom_scientifique == "Hibiscus tiliaceus R" |
  BDD_moy_esp$Nom_scientifique== "Acropogon bullatus" | 
  BDD_moy_esp$Nom_scientifique== "Barringtonia asiatica" | 
  BDD_moy_esp$Nom_scientifique== "Alphitonia neocaledonica" |
  BDD_moy_esp$Nom_scientifique== "Leucaena leucocephala" |
  BDD_moy_esp$Nom_scientifique== "Morinda citrifolia" |
  BDD_moy_esp$Nom_scientifique== "Cerbera manghas" |
  BDD_moy_esp$Nom_scientifique== "Acalypha sp"

points (BDD_moy_esp$LMC_t0[esp_a_annoter],
        BDD_moy_esp$LMC_t24[esp_a_annoter],pch=21,bg = rgb(1,0,0),cex=1.2
)

text(BDD_moy_esp$LMC_t0[esp_a_annoter] + 2,
     BDD_moy_esp$LMC_t24[esp_a_annoter] + 15,
     labels = c("A. sp","A. bullatus","A. neocaledonica","B. asiatica","C. manghas","H. tiliaceus R","L. leucocephala","M. citrifolia","S. montana"),
     cex = 1)
























########################## PLOT GMIN LMC_t0 ###########################

M_L2<-aggregate(BDD_moy_ech$LMC_t0,by=list(BDD_moy_ech$Nom_scientifique),mean)
SD_L2<-aggregate(BDD_moy_ech$LMC_t0,by=list(BDD_moy_ech$Nom_scientifique),sd)

M_L1<-aggregate(BDD_moy_ech$Gmin,by=list(BDD_moy_ech$Nom_scientifique),mean)
SD_L1<-aggregate(BDD_moy_ech$Gmin,by=list(BDD_moy_ech$Nom_scientifique),sd)



# Création du scatterplot
plot(BDD_moy_esp$LMC_t0 ~ BDD_moy_esp$Gmin, type="n",
     xlab = "Gmin", 
     ylab = "LMC_t0 (%)", 
     main = "Relation entre Gmin et LMC_t0",ylim=c(0,450),xlim=c(0,100))

esp_a_annoter <- BDD_moy_esp$Nom_scientifique == "Scaevola montana" | 
  BDD_moy_esp$Nom_scientifique == "Hibiscus tiliaceus R" |
  BDD_moy_esp$Nom_scientifique== "Acropogon bullatus" |
  BDD_moy_esp$Nom_scientifique== "Leucaena leucocephala" | 
  BDD_moy_esp$Nom_scientifique== "Cerbera manghas" | 
  BDD_moy_esp$Nom_scientifique== "Diospyros fasciculosa" |
  BDD_moy_esp$Nom_scientifique== "Calophylum caledonicum" |
  BDD_moy_esp$Nom_scientifique== "Alphitonia neocaledonica" |
  
  segments(x0=BDD_moy_esp$Gmin[esp_a_annoter],x1=BDD_moy_esp$Gmin[esp_a_annoter],y0=BDD_moy_esp$LMC_t0[esp_a_annoter]-SD_L2$x[esp_a_annoter],y1=BDD_moy_esp$LMC_t0[esp_a_annoter]+SD_L2$x[esp_a_annoter])
segments(x0=BDD_moy_esp$Gmin[esp_a_annoter]-SD_L1$x[esp_a_annoter],x1=BDD_moy_esp$Gmin[esp_a_annoter]+SD_L1$x[esp_a_annoter],y0=BDD_moy_esp$LMC_t0[esp_a_annoter],y1=BDD_moy_esp$LMC_t0[esp_a_annoter])
points(M_L1$x,M_L2$x,pch=19,cex=0.8)
points (BDD_moy_esp$Gmin[esp_a_annoter],
        BDD_moy_esp$LMC_t0[esp_a_annoter],pch=21,bg = rgb(1,0,0),cex=1.2
)


text(BDD_moy_esp$Gmin[esp_a_annoter],
     BDD_moy_esp$LMC_t0[esp_a_annoter],
     labels = BDD_moy_esp$Nom_scientifique[esp_a_annoter],
     pos = 4, offset = 0.5, cex = 0.8)

points (BDD_moy_esp$Gmin[esp_a_annoter],
        BDD_moy_esp$LMC_t0[esp_a_annoter],pch=21,bg = rgb(1,0,0)
)

plot(BDD_moy_esp$LMC_t0,BDD_moy_esp$LMC_t24 ,xlim=c(0,300),ylim=c(0,300))
text(BDD_moy_esp$LMC_t0,BDD_moy_esp$LMC_t24,
     labels = BDD_moy_esp$Nom_scientifique,cex=1,pos=2)
















































##### TEST LMC_t24 ###########
hist(BDD_moy_esp$LMC_t24)

BDD_moy_esp$LDMC_cr<-as.numeric(scale(BDD_moy_esp$LDMC))
BDD_moy_esp$Gmin_cr<-as.numeric(scale(BDD_moy_esp$Gmin))
BDD_moy_esp_netMT$LDMC_cr<-as.numeric(scale(BDD_moy_esp_netMT$LDMC))
BDD_moy_esp_netMT$Gmin_cr<-as.numeric(scale(BDD_moy_esp_netMT$Gmin))
BDD_moy_esp_net3$LDMC_cr<-as.numeric(scale(BDD_moy_esp_net3$LDMC))
BDD_moy_esp_net3$Gmin_cr<-as.numeric(scale(BDD_moy_esp_net3$Gmin))

#modèles
m<-glm(BDD_moy_esp$LMC_t24~LDMC_cr+Gmin_cr,data=BDD_moy_esp,family=Gamma(link="log"))
summary(m)


################# Prédiction ########################
# LDMC
moy <- mean(BDD_moy_esp$LDMC, na.rm = TRUE)
ecart <- sd(BDD_moy_esp$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_moy_esp$LDMC_cr), max(BDD_moy_esp$LDMC_cr),length.out = 100)

# on fait varier Gmin et on fixe LDMC
pred_data1 <- data.frame(
  LDMC_cr = foo,
  Gmin_cr = rep(mean(BDD_moy_esp$Gmin_cr,na.rm = TRUE), length(foo))
)

pred_data2 <- data.frame(
  LDMC_cr = foo,
  Gmin_cr = rep(min(BDD_moy_esp$Gmin_cr,na.rm = TRUE), length(foo))
)

pred_data3 <- data.frame(
  LDMC_cr = foo,
  Gmin_cr = rep(max(BDD_moy_esp$Gmin_cr,na.rm = TRUE), length(foo))
)




# Prédictions avec m
pred1 <- predict(m, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred2 <- predict(m, type = "response", newdata = pred_data2, se.fit = TRUE, re.form = NA)
pred3 <- predict(m, type = "response", newdata = pred_data3, se.fit = TRUE, re.form = NA)


# Plot des points observés
plot(BDD_moy_esp$LDMC, BDD_moy_esp$LMC_t24, 
     xlab = "LDMC (mg/g)", ylab = "LMC_t24 (%)", 
     main = "Effet de LDMC sur LMC_t24 selon gmin ",ylim=c(0,700))

# foo (centré-réduit) en valeurs originales de LDMC pour l'axe X
foo_original <- foo * ecart + moy

# Courbes de prédiction 
lines(foo_original, pred1$fit, col = "#3B393B", lwd = 2)  # Gmin moyen
lines(foo_original, pred2$fit, col = "#005ED1", lwd = 2)  # Gmin min
lines(foo_original, pred3$fit, col = "#BD0000", lwd = 2)  # Gmin max













############################ MODELES échelle esp ###############################

#standardisation des données 
BDD_moy_esp_net$SD_cr<-as.numeric(scale(BDD_moy_esp_net$SD))
BDD_moy_esp_net$TD_cr<-as.numeric(scale(BDD_moy_esp_net$TD))
BDD_moy_esp_net$LA_cr<-as.numeric(scale(BDD_moy_esp_net$Surface_F))
BDD_moy_esp_net$LDMC_cr<-as.numeric(scale(BDD_moy_esp_net$LDMC))
BDD_moy_esp_net$LMC_t24_cr<-as.numeric(scale(BDD_moy_esp_net$LMC_t24))
BDD_moy_esp_net$Gmin_cr<-as.numeric(scale(BDD_moy_esp_net$Gmin))
BDD_moy_esp_net$SLA_cr<-as.numeric(scale(BDD_moy_esp_net$SLA))
BDD_moy_esp_net$Nb_rami_cr<-as.numeric(scale(BDD_moy_esp_net$Nb_rami))


BDD_moy_esp_netMT$SD_cr<-as.numeric(scale(BDD_moy_esp_netMT$SD))
BDD_moy_esp_netMT$TD_cr<-as.numeric(scale(BDD_moy_esp_netMT$TD))
BDD_moy_esp_netMT$LA_cr<-as.numeric(scale(BDD_moy_esp_netMT$Surface_F))
BDD_moy_esp_netMT$LDMC_cr<-as.numeric(scale(BDD_moy_esp_netMT$LDMC))
BDD_moy_esp_netMT$LMC_t24_cr<-as.numeric(scale(BDD_moy_esp_netMT$LMC_t24))
BDD_moy_esp_netMT$Gmin_cr<-as.numeric(scale(BDD_moy_esp_netMT$Gmin))
BDD_moy_esp_netMT$SLA_cr<-as.numeric(scale(BDD_moy_esp_netMT$SLA))
BDD_moy_esp_netMT$Nb_rami_cr<-as.numeric(scale(BDD_moy_esp_netMT$Nb_rami))

BDD_moy_esp_net3$SD_cr<-as.numeric(scale(BDD_moy_esp_net3$SD))
BDD_moy_esp_net3$TD_cr<-as.numeric(scale(BDD_moy_esp_net3$TD))
BDD_moy_esp_net3$LA_cr<-as.numeric(scale(BDD_moy_esp_net3$Surface_F))
BDD_moy_esp_net3$LDMC_cr<-as.numeric(scale(BDD_moy_esp_net3$LDMC))
BDD_moy_esp_net3$LMC_t24_cr<-as.numeric(scale(BDD_moy_esp_net3$LMC_t24))
BDD_moy_esp_net3$Gmin_cr<-as.numeric(scale(BDD_moy_esp_net3$Gmin))
BDD_moy_esp_net3$SLA_cr<-as.numeric(scale(BDD_moy_esp_net3$SLA))
BDD_moy_esp_net3$Nb_rami_cr<-as.numeric(scale(BDD_moy_esp_net3$Nb_rami))

BDD_moy_score$SD_cr<-as.numeric(scale(BDD_moy_score$SD))
BDD_moy_score$TD_cr<-as.numeric(scale(BDD_moy_score$TD))
BDD_moy_score$LA_cr<-as.numeric(scale(BDD_moy_score$Surface_F))
BDD_moy_score$LDMC_cr<-as.numeric(scale(BDD_moy_score$LDMC))
BDD_moy_score$LMC_t24_cr<-as.numeric(scale(BDD_moy_score$LMC_t24))
BDD_moy_score$Gmin_cr<-as.numeric(scale(BDD_moy_score$Gmin))
BDD_moy_score$SLA_cr<-as.numeric(scale(BDD_moy_score$SLA))
BDD_moy_score$Nb_rami_cr<-as.numeric(scale(BDD_moy_score$Nb_rami))



#######################################################
###### MT #########
#######################################################

#modèles
MT1 <- glm (BDD_moy_esp_netMT$MT ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LDMC_cr*Gmin_cr, data = BDD_moy_esp_netMT)
summary(MT1)


###################### prédicion ############### 

############  LDMC
moy <- mean(BDD_moy_esp_netMT$LDMC, na.rm = TRUE)
ecart <- sd(BDD_moy_esp_netMT$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_moy_esp_netMT$LDMC_cr), max(BDD_moy_esp_netMT$LDMC_cr),length.out = 100)

############## on fait varier LMDC
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_esp_netMT$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_esp_netMT$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_esp_netMT$Nb_rami_cr), length(foo)),
  Gmin_cr = rep(mean(BDD_moy_esp_netMT$Gmin_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_esp_netMT$SLA_cr), length(foo))
)


# Prédictions
pred1 <- predict(MT1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
par(mar=c(4,5,4,4))
plot(BDD_moy_esp_netMT$LDMC, BDD_moy_esp_netMT$MT,type="n", 
     xlab = "LDMC (mg/g)", ylab = "Température maximum (°C)", 
     main = "Effet de LDMC sur MT",ylim=c(500,850),xlim=c(175,550))

# Courbes de prédiction
lines((foo*ecart+moy), pred1$fit, col = "black", lwd = 2) 

# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)








##################################
###### BB #########
##################################

#modèles
BDD_moy_esp_net3$BB_prop <- BDD_moy_esp_net3$BB_test/100


BB0 <- glm(BB_prop ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LDMC_cr*Gmin_cr, data = BDD_moy_esp_net3)    #### meilleur modèle 
summary(BB0)
r.squaredGLMM(BB0)



################# Prédiction ########################
# LDMC
moy <- mean(BDD_moy_esp_net3$LDMC, na.rm = TRUE)
ecart <- sd(BDD_moy_esp_net3$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_moy_esp_net3$LDMC_cr), max(BDD_moy_esp_net3$LDMC_cr),length.out = 100)

# on fait varier LMDC 
par(mar=c(4,5,4,4))

pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_esp_net3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_esp_net3$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_esp_net3$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_esp_net3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(BB1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_moy_esp_net3$LDMC, BDD_moy_esp_net3$BB_prop,type="n", 
     xlab = "LDMC", ylab = "Biomasse brûlée (%)", 
     main = "Effet de LDMC sur BB ",ylim=c(0,1),xlim=c(175,550))

# Courbes de prédiction
lines((foo*ecart + moy), pmin(pmax(pred1$fit, 0), 1), col = "black", lwd = 2 )


# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)







#######################################
###### BT #########
######################################

#modèles
BT0 <- glm(BT_test ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LDMC_cr*Gmin_cr , family=Gamma(link="log"), data = BDD_moy_esp_net3)  #meilleur modèle
summary(BT0)



################# Prédiction ########################
# LDMC
moy <- mean(BDD_moy_esp_net3$LDMC, na.rm = TRUE)
ecart <- sd(BDD_moy_esp_net3$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_moy_esp_net3$LDMC_cr), max(BDD_moy_esp_net3$LDMC_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_esp_net3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_esp_net3$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_esp_net3$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_esp_net3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_BT1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)
pred1

# Plot des points observés
plot(BDD_moy_esp_net3$LDMC, BDD_moy_esp_net3$BT_test,type="n", 
     xlab = "LDMC", ylab = "Temps de combustion (s)", 
     main = "Effet de LDMC sur BT ",ylim=c(0,80),xlim=c(175,550))

# Courbes de prédiction
lines((foo*ecart + moy), pred1$fit, col = "black", lwd = 2 )

# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), (pred1$fit + 1.96 * pred1$se.fit), col = "blue", lty = 3)
lines((foo*ecart + moy), (pred1$fit - 1.96 * pred1$se.fit), col = "blue", lty = 3)

range(pred1$se.fit)
range(pred1$fit + 1.96 * pred1$se.fit - (pred1$fit - 1.96 * pred1$se.fit))



###################
###### DI ######################################################################
###################

#modèles

DI0 <- glm(DI_test ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LDMC_cr+Gmin_cr , family=Gamma(link="log"), data = BDD_moy_esp_net3)  #meilleur modèle
summary(DI0)



################# Prédiction ########################
# LDMC
moy <- mean(BDD_moy_esp_net3$LDMC, na.rm = TRUE)
ecart <- sd(BDD_moy_esp_net3$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_moy_esp_net3$LDMC_cr), max(BDD_moy_esp_net3$LDMC_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_esp_net3$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_esp_net3$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_esp_net3$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_esp_net3$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_DI1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

# Plot des points observés
plot(BDD_moy_esp_net3$LDMC, BDD_moy_esp_net3$DI_test,type="n", 
     xlab = "LDMC", ylab = "Délai d'ignition (s)", 
     main = "Effet de LDMC sur DI ",ylim=c(0.5,2),xlim=c(175,550))

# Courbes de prédiction
lines((foo*ecart + moy), pred1$fit, col = "black", lwd = 2 )


# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)




#####################
###### Score ########
#####################

#modèles
score1 <- glm(score ~ SD_cr+TD_cr+SLA_cr+Nb_rami_cr+LDMC_cr+Gmin_cr, data = BDD_moy_score)
summary(score1)



################# Prédiction ########################
# LDMC
moy <- mean(BDD_moy_ech$LDMC, na.rm = TRUE)
ecart <- sd(BDD_moy_ech$LDMC, na.rm = TRUE)

# Valeurs de LDMC
foo <- seq(min(BDD_moy_ech$LDMC_cr), max(BDD_moy_ech$LDMC_cr),length.out = 100)

# on fait varier LMDC 
pred_data1 <- data.frame(
  LDMC_cr = foo,
  SD_cr = rep(mean(BDD_moy_ech$SD_cr,na.rm = TRUE), length(foo)),
  TD_cr = rep(mean(BDD_moy_ech$TD_cr), length(foo)),
  Nb_rami_cr = rep(mean(BDD_moy_ech$Nb_rami_cr), length(foo)),
  SLA_cr = rep(mean(BDD_moy_ech$SLA_cr), length(foo))
)

# Prédictions
pred1 <- predict(m_score1, type = "response", newdata = pred_data1, se.fit = TRUE, re.form = NA)

#background transparent 
par(bg = NA)

# Plot des points observés
plot(BDD_moy_ech$LDMC, BDD_moy_ech$score,type="n",cex.axis=1.2, 
     xlab = "LDMC", ylab = "Flammability score", 
     main = "Effet de LDMC sur score ",ylim=c(-4,4),xlim=c(150,600))
box(lwd = 2)
points(BDD_moy_ech$LDMC, BDD_moy_ech$score, pch=16, col="#634124")

col_ic <- rgb(109/255, 168/255, 111/255)

# Polygone de l'intervalle de confiance
x_vals <- (foo * ecart + moy)

polygon(
  x = c(x_vals, rev(x_vals)),
  y = c(pred1$fit + 1.96 * pred1$se.fit,
        rev(pred1$fit - 1.96 * pred1$se.fit)),
  col = col_ic,
  border = NA
)




# Courbes de prédiction
lines((foo*ecart + moy), pred1$fit, col = "#C94802", lwd = 3 )
lines((foo*ecart + moy), pred1$fit,  lwd = 3 )

# Intervalle de confiance pour SD moyen
lines((foo*ecart + moy), pred1$fit + 1.96 * pred1$se.fit, col = "blue", lty = 3)
lines((foo*ecart + moy), pred1$fit - 1.96 * pred1$se.fit, col = "blue", lty = 3)








