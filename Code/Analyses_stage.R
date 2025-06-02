
#### ANALYSE DE DONNEE STAGE INFLAMMABILITE ####

##  importation BDD_ana_ech en format CSV
setwd("C:/IRD/Stage_inflammabilit-") #définition du répertoire de travail
BDD_ech<-read.csv2("Data/BDD_ana_ech.csv", header = TRUE) #importation de la base
BDD_ech
dim(BDD_ech)

BDD_sd_esp<-read.csv2("Data/BDD_sd_esp.csv", header = TRUE)
BDD_sd_esp

##  importation BDD_moy_esp en format CSV
BDD_esp<-read.csv2("Data/BDD_moy_esp.csv", header = TRUE) #importation de la base
BDD_esp
names(BDD_esp)[which(names(BDD_esp) == "Nb_ramifications")] <- "Nb_rami"
dim(BDD_esp)



#ACP choix assollement 

# Charger les bibliothèques nécessaires
library(FactoMineR)
library(factoextra)


################  ACP INFLAMMABILITE ####################

# Sélectionner uniquement les colonnes des pourcentages
colonnes_infla <- BDD_esp[, c(3:6)]

# Vérifier les données
head(colonnes_infla)

# Appliquer l'ACP
res.pca <- PCA(colonnes_infla, scale.unit = TRUE, graph = FALSE)

# Résumé des résultats
summary(res.pca)

# Graphique des contributions des variables aux composantes principales
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#6fec00", "#ff9e00", "#Ff0000"),
             repel = TRUE)


# Graphique des individus 
fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)

# Graphique combiné des variables et des individus
fviz_pca_biplot(res.pca,col.var = "contrib",
                gradient.cols =c("#00AFBB", "#E7B800", "#FC4E07") , repel = TRUE)

# Afficher l'ébouli
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50), main="Graphique de l'ébouli")

# Ajout du score d'inflammabilité basé sur la coordonnée de l'axe 1 ou axe 1 et 2
coord <- res.pca$ind$coord  # coordonnées des individus
BDD_esp$score <- coord[, 1]  # Dim 1 = axe 1
BDD_esp$score2 <- (0.8*coord[,1] + 0.2*coord[,2]) / 2

min_score <- min(BDD_esp$score2, na.rm = TRUE)
max_score <- max(BDD_esp$score2, na.rm = TRUE)

BDD_esp$score_normalise <- -1 + (BDD_esp$score2 - min_score) * 2 / (max_score - min_score)


View(BDD_esp)






#################### ACP TRAITS ################
# Sélectionner uniquement les colonnes des pourcentages
colonnes_traits <- na.omit(BDD_esp[,c(7:22)])

# Vérifier les données
colonnes_traits

# Appliquer l'ACP
res.pca <- PCA(colonnes_traits, scale.unit = TRUE, graph = FALSE)

# Résumé des résultats
summary(res.pca)

# Graphique des contributions des variables aux composantes principales
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)


# Graphique des individus 
fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)

# Graphique combiné des variables et des individus
fviz_pca_biplot(res.pca, repel = TRUE)

# Afficher l'ébouli
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50), main="Graphique de l'ébouli")




############ ACP de tout ###########

# Sélectionner uniquement les colonnes des pourcentages
colonnes_all <-na.omit(BDD_esp[,c(3:21)])

# Vérifier les données
colonnes_all

# Appliquer l'ACP
res.pca <- PCA(colonnes_all, scale.unit = TRUE, graph = FALSE)

# Résumé des résultats
summary(res.pca)

# Graphique des contributions des variables aux composantes principales
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)


# Graphique des individus (facultatif, si vous analysez chaque observation)
fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)

# Graphique combiné des variables et des individus
fviz_pca_biplot(res.pca, repel = TRUE)

# Afficher l'ébouli
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50), main="Graphique de l'ébouli")





#################### PLOT ##################################
hist(BDD_esp$score2)
text(BDD_esp$Nom_scientifique)

### LMC-t24 avec les 4 composantes infla ####
plot(BDD_esp$LMC_t24,BDD_esp$BT)
text(BDD_esp$LMC_t24,BDD_esp$BT,BDD_esp$Nom_scientifique)

plot(BDD_esp$LMC_t24,BDD_esp$BB)
text(BDD_esp$LMC_t24,BDD_esp$BB,BDD_esp$Nom_scientifique)

plot(BDD_esp$LMC_t24,BDD_esp$MT)
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






####################### matrice de corrélation ######################

library (corrplot)
mat_cor<-cor(colonnes_all)
mat_sub <- mat_cor[1:4, ]
corrplot(mat_sub)
## afficher les valeurs
round(mat_sub, 2)






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

# Inverser DI 
df_prep$DI <- max(df_prep$DI, na.rm = TRUE) - df_prep$DI

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

