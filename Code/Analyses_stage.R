
#### ANALYSE DE DONNEE STAGE INFLAMMABILITE ####

##  importation BDD_ana_ech en format CSV
setwd("C:/IRD/Stage_inflammabilit-") #définition du répertoire de travail
BDD_ech<-read.csv2("Data/BDD_ana_ech.csv", header = TRUE) #importation de la base
BDD_ech
dim(BDD_ech)

##  importation BDD_moy_esp en format CSV
BDD_esp<-read.csv2("Data/BDD_moy_esp.csv", header = TRUE) #importation de la base
BDD_esp
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


# Graphique des individus (facultatif, si vous analysez chaque observation)
fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)

# Graphique combiné des variables et des individus
fviz_pca_biplot(res.pca, repel = TRUE)

# Afficher l'ébouli
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50), main="Graphique de l'ébouli")

# Ajouter un score d'inflammabilité basé sur la coordonnée de l'axe 1 ou axe 1 et 2
coord <- res.pca$ind$coord  # coordonnées des individus
BDD_esp$score <- coord[, 1]  # Dim 1 = axe 1
BDD_esp$score2 <- (0.63*coord[,1] + 0.25*coord[,2]) / 2
View(BDD_esp)

#################### ACP TRAITS ################
# Sélectionner uniquement les colonnes des pourcentages
colonnes_traits <- na.omit(BDD_esp[,c(7:21)])

# Vérifier les données
colonnes_traits

# Appliquer l'ACP
res.pca <- PCA(colonnes_traits, scale.unit = TRUE, graph = FALSE)

# Résumé des résultats
summary(res.pca)

# Graphique des contributions des variables aux composantes principales
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#0536d0", "#05d096", "#08d005"),
             repel = TRUE)


# Graphique des individus 
fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)

# Graphique combiné des variables et des individus
fviz_pca_biplot(res.pca, repel = TRUE)

# Afficher l'ébouli
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50), main="Graphique de l'ébouli")




############ ACP de tout ###########

# Sélectionner uniquement les colonnes des pourcentages
colonnes_all <- BDD_esp[,c(3:21)]

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





#### PLOT ######
plot(BDD_esp$score,BDD_esp$SLA)
m<-lm(BT~LT,data=BDD_ech)
summary(m)
plot(m)
anova(m)

hist(BDD_ech$MT)
hist(BDD_ech$BT)
hist(BDD_ech$DI)
hist(BDD_ech$BB)
