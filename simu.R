library(Rhabit)
library(viridis)
library(ggplot2)
library(gridExtra)
library(abind)
library(dplyr)
library(RandomFields)
library(RandomFieldsUtils)
library(markovchain)
library(tidyverse)
library(RColorBrewer)
library(MASS)
library(mvtnorm)
library(psych)
library(FactoMineR)
library(stringr)


################################################################################
###                Générateur de la chaîne des états cachés                  ###                               
################################################################################
# Fonction permettant de générer une chaîne de Markov. Elle est généralisée de 
# façon à ce qu'on n'ait pas à donner la suite de "noms" des états.
# Paramètres : 
#           - A, la matrice de transition de la chaîne au format matrix. 
#           - T, le nombre souhaité d'itérations de la chaîne.
#
# Retourne une suite numérique des états successifs.

CM_generateur = function(A, nbr_obs){
  K = dim(A)[1]
  CM<-new("markovchain", states = as.character(lapply(1:K, as.character)),
          transitionMatrix = A, 
          name = "Q")
  return(as.numeric(markovchainSequence( nbr_obs, CM)))
}


################################################################################
###                        Création des paramètres                           ###                               
################################################################################

BETA = function(K, J){
  # K : nb états dans la chaine.
  # J : nombre de covariables.
  n = K*J
  liste = runif(n,-5,5)
  return(matrix(liste, ncol = K, nrow = J))
}
Nu = function(BETA, vit) { 
  # BETA : une matrice contenant les paramètres beta du modèle. De taille 
  #        J * K.
  # vit : liste contenant les vitesses du processus dans chacun des états. Donc 
  #       de taille K.
  l = sweep(BETA, 2, STATS = vit**2, FUN =  "*")
  return(l)
}
Nu(matrix(1,2,2),c(0.4,0.2))
BetaToNu = function(Beta, Vit){
  # Fonction qui prend en paramètre les betas et les gamma selon les etats et qui 
  # doit renvoyer la matrice theta correspondante.
  Nu =  sweep(list_to_matrix(Beta), 2, STATS = Vit**2, FUN =  "*")
  return(Nu)
}
b = matrix_to_list(Nu(BETA(2,2),0.4))
#BetaToNu(matrix_to_list(Nu(BETA(2,2),0.4)),c(0.4,0.3))
#theta = Nu(BETA(K,J), vit = c(.4,.38))

################################################################################
###                     Probabilités des apparitions                         ###                               
################################################################################
# Fonction calculant la probabilité des émissions, c'est à dire la probabilité 
# de chaque "partie" du déplacement de l'animal. Il prend en compte la 
# dimension 1 et 2.

# Paramètres :
#         - obs, le data.frame contenant les informations sur le déplacement.
#         - C, la matrice des covariables environnementales.
## C est utilisé juste au desuus pouyr Beta * Vit
#         - theta, le paramètre de lien.
#         - Delta, la suite des pas de temps.
#         - vit, la "vitesse" du processus aléatoire.
#         - dimension, la dimension considérée, elle est de base égale à 2. 
#
# Retourne un dataframe avec deux composantes :
#         - B, la matrice des probabilités des déplacements, T lignes, K 
#           colonnes.
#
# Remarques :
#         - Seules les dimensions 1 et 2 sont prises en compte ici.


proba_emission = function(obs, C, theta, Delta, Vits, dimension = 2){

  nbr_obs = length(obs$Z1) - 1  # On enleve 1 pour ne pas prendre le NA.
  K = dim(theta)[2]
  
  # Création de la matrice.
  B = matrix(0, ncol = K, nrow = nbr_obs)
  if (dimension == 2){
    # On implémente Z sous un format pratique pour la suite.
    Z = cbind(obs$Z1[1:nbr_obs], obs$Z2[1:nbr_obs])
    
    # Remplissage de la matrice.
    for (t in 1:nbr_obs){
      C_tilde = matrix( c(C[t,],C[nbr_obs + t,]), nrow=2, byrow=TRUE)
      mu =  (C_tilde %*% theta) / sqrt(Delta[t]) 
      for (k in 1:K){B[t,k] = dmvnorm( Z[t,], mu[,k], Vits[k]**2 * diag(2))}
    }
    # lapply( 1:K, function(d){dmvnorm(Z[t,], mu[,d], Vits[d]**2 * diag(2))}) ?
    # Peut être un lapply de lapply ?
    # lapply(1:nbr_obs, function(t){ C_tilde; Mu; lapply(1:K, function(d){dmvnorm(Z[t,], mu[,d], Vits[d]**2 * diag(2))}) 
    #})
  }
  else {
    # Remplissage de la matrice.
    for (t in 1:nbr_obs){
      for (k in 1:K){B[t,k] = dnorm(obs$Z[t],Delta[t]*(C[t,] %*% theta[,k]),vit)}
    }
  }
  return(B)
}
#proba_emission(Obs, C, theta, incr, )
################################################################################
###                      Génération des observations                         ###                               
################################################################################

Obs_generateur = function(Q, theta, C, vit, delta, dimension = 2){
  nbr_obs = length(Q) - 2 
  
  if (dimension == 2){
    # On considère que l'animal part de l'origine.
    X1 = c(0)
    X2 = c(0)
    Z1 = c(0)
    Z2 = c(0)
    
    for (i in 2:nbr_obs) {
      # MU.
      C_tilde = matrix(c(C[i,],C[nbr_obs+i,]),nrow=2,byrow=TRUE) # Matrice C modifiée.
      mu = C_tilde%*%theta[,Q[i]]
      
      # Sigma.
      #Sigma = (vit**2) * delta[i-1]* diag(2) 
      Sigma = (vit**2) * diag(2) 
      
      # Calcul du déplacement.
      e = mvrnorm(1,  mu/sqrt(delta[i-1]), Sigma)
      #e = mvrnorm(1,  mu * delta[i-1], Sigma)
      if(abs(e[1]) > 4 ){ browser()}
      # Stockage des données.
      X1 = c(X1, X1[i-1] + e[1])
      X2 = c(X2, X2[i-1] + e[2])
      Z1 = c(Z1, e[1])
      Z2 = c(Z2, e[2])
    }
    return(data.frame('X1' = X1, 'X2' = X2, 'Z1' = Z1, 'Z2' = Z2))
  }
  
  else {
    # On considère que l'animal part de l'origine.
    X = c(0) # Départ déterministe en l'origine.
    Z = c(0)
    
    # Sigma. 
    Sigma = delta * vit**2
    
    for (i in 2:nbr_obs) {
      # Mu.
      mu = C[i,]%*%theta[,Q[i]]
      
      # Calcul du déplacement.
      e = rnorm(1, delta[i-1]*mu, Sigma[i-1])
      
      # Stockage des données.
      X = c(X, X[i-1] + e)
      Z = c(Z, e)
    }
    return(data.frame('X' = X, 'Z' = Z))
  }
}

Generation_observation = function(theta, Q, liste_cov, vit, tps, loc0 = c(0,0)){
  K = dim(theta)[2]
  nbr_obs = length(Q)
  Obs <- matrix(NA, nbr_obs, 2)
  Obs[1,] = loc0
  for (t in (1:(nbr_obs-1))) {
    print(t)
    # On regarde l'état dans lequel on se trouve.
    etat = Q[t]
    print(etat)
    # On simule le déplacement.
    xy = simLangevinMM(matrix_to_list(theta)[[etat]], vit, c(tps[t],tps[t+1]),
                       loc0 = c(0,0), cov_list = liste_cov, keep_grad = FALSE)
    Obs[t,] = as.vector(as.numeric(xy[2,][,1:2])) # On ne prend que les coordonnées du déplacement.
    loc0 = as.vector(as.numeric(xy[2,][,1:2]))
  }
  return(Obs)}

Generation_observation2.0 = function(beta, Q, C, vit, time, loc0 = c(0,0), affichage = TRUE){
  K = dim(theta)[2]
  t = length(Q)
  print(t)
  Obs <- matrix(NA, t-1, 2)
  Obs[1,] = loc0
  for (t in 1:(t-1)) {
    xy = simLangevinMM(beta[[Q[t]]], vit, c(tps[t],tps[t+1]), c(0,0), liste_cov, keep_grad = FALSE)
    Obs[t,] = as.vector(as.numeric(xy[2,][,1:2])) # On ne prend que les coordonnées du déplacement.
    loc0 = as.vector(as.numeric(xy[2,][,1:2]))
  }
  
  # On s'occupe du format de retour. 
  
  # Les Z. 
  
  
  # Les X. 
  Observations = data.frame('X1' = Obs[,1], 'X2' = Obs[,2],
                            'Z1' = c(diff(Obs[,1]), NA), 
                            'Z2' = c(diff(Obs[,2]), NA))
  # if (affichage){
  #   Obs_affichage = Observations
  #   Obs_affichage$'etats' = as.factor(Q[1:t-1])
  #   
  #   ggplot(Obs_affichage, aes(X1, X2)) +
  #     geom_path()+
  #     geom_point(aes(colour = etats)) 
  # }
  
  return(Observations)}

Generation_observation3.0 = function(liste_theta, etats_caches, liste_cov, Vits, tps, loc0 = c(0,0), affichage = TRUE, epsilon = 0.5){
  
  minx = min(liste_cov[[1]]$x) + 4 *epsilon
  miny = min(liste_cov[[1]]$y) + 4 *epsilon
  maxx = max(liste_cov[[1]]$x) - 4 *epsilon
  maxy = min(liste_cov[[1]]$y) - 4 *epsilon
  K = dim(theta)[2]
  J = length(liste_cov)
  nbr_obs = length(etats_caches) 
  simu <-  simLangevinMM(liste_theta[[ etats_caches[1] ]], Vits[etats_caches[1]], c(tps[1],tps[2]), loc0 = loc0, liste_cov, keep_grad = TRUE)
  for (t in 2:nbr_obs) {
    prov <- simLangevinMM(liste_theta[[ etats_caches[t] ]], Vits[etats_caches[t]], c(tps[t-1],tps[t]), loc0 = as.numeric(simu[t-1,1:2]), liste_cov, keep_grad = TRUE)
    if(prov[2,1]<minx) {
      prov[2,1]= minx + epsilon
    }
    if(prov[2,1]>maxx)  {
      prov[2,1] = maxx - epsilon
    }
    if(prov[2,2]<miny)  {
      prov[2,2] = miny + epsilon
    }
    if(prov[2,2]<minx)  {
      prov[2,2] = maxy - epsilon
    }
    
    simu[t-1, 4:(2*J+3)] <-prov[1,4:(2*J+3)]
    simu[t, ] <- prov[2,]
  }
  
  # Potentiellement mettre un 0 au début comme avant plutôt qu'un NA à la fin, ce serait plus propre. 
  simu <- simu %>% mutate(Z1 = x-lag(x), Z2 = y -lag(y)) %>% mutate(etats_caches = etats_caches)
  return(simu)}




