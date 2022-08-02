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



################################################################################
#                                                                              #
#                         Les différentes fonctions                            #
#                                                                              #
################################################################################       

################################################################################
###                       Quelques fonctions pratiques                       ###                               
################################################################################
# Fonction permettant de rendre stochastique une matrice devant l'être.
# Utilisée dans l'EM.
format_mat = function(A){
  l = ncol(A)
  for (i in 1:l){
    A[i,] = A[i,]/(sum(A[i,]))
  }
  return(A)
}

increments = function(liste){
  l = c()
  for (i in 1:(length(liste)-1)){
    l = c(l,liste[i+1] - liste[i])
  }
  return(l)
}

matrix_to_list = function(mat){
  l = list()
  N = ncol(mat)
  for (i in 1:N){
    l[[i]] = c(mat[,i])
  }
  return(l)
}


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

CM_generateur = function(A, T){
  
  # On commence par générer la suite "states".
  K = dim(A)[1]
  states = c()
  for (i in 1:K){states = c(states, as.character(i))}
  
  # On génère maintenant la chaîne de Markov. 
  CM<-new("markovchain", states = states,
          transitionMatrix = A, 
          name = "Q")
  return(as.numeric(markovchainSequence(T,CM)))
}





################################################################################
###                        Création des paramètres                           ###                               
################################################################################

BETA = function(K, J){
  n = K*J
  liste = runif(n,-5,5)
  return(matrix(liste, ncol = K, nrow = J))
}
Nu = function(BETA, delta, vit) { 
  
  # Les constantes utiles.
  dim = dim(BETA)
  J = dim[1]
  K = dim[2]
  
  # Création de la matrice.
  L = matrix(1, ncol = K, nrow = J)
  
  # Remplissage.
  for (k in 1:K){
    for (j in 1:J){
      L[j,k] = (BETA[j,k]*sqrt(delta)*vit**2)
    }
  }
  return(L)
}



################################################################################
###                     Probabilités des apparitions                         ###                               
################################################################################
# Fonction calculant la probabilité des émissions, c'est à dire la probabilité 
# de chaque "partie" du déplacement de l'animal. Il prend en compte la 
# dimension 1 et 2.

# Paramètres :
#         - obs, le data.frame contenant les informations sur le déplacement.
#         - C, la matrice des covariables environnementales.
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


proba_emission = function(obs, C, theta, Delta, vit, dimension = 2){
  T = dim(obs)[1]
  K = dim(theta)[2]
  
  # Création de la matrice.
  B = matrix(0, ncol = K, nrow = T)
  
  if (dimension == 2){
    # On implémente Z sous un format pratique pour la suite.
    Z = cbind(obs$Z1, obs$Z2)
    
    # Remplissage de la matrice.
    for (t in 1:T){
      C_tilde = matrix(c(C[t,],C[T+t,]),nrow=2,byrow=TRUE)
      #mu = Delta[t] * C_tilde %*% theta
      mu =  C_tilde %*% theta / sqrt(Delta[t]) 
      #Sigma = vit**2 * Delta[t] * diag(2)
      Sigma = vit**2 * diag(2)
      for (k in 1:K){B[t,k] = dmvnorm(Z[t,], mu[,k], Sigma)}
      # if (B[t,] = matrix(0,1,K)) {B[t,] = matrix(3.53e-300,1,K) }
    }
  }
  else {
    # Remplissage de la matrice.
    for (t in 1:T){
      for (k in 1:K){B[t,k] = dnorm(obs$Z[t],Delta[t]*(C[t,] %*% theta[,k]),vit)}
    }
  }
  
  
  
  return(B)
}



################################################################################
###                      Génération des observations                         ###                               
################################################################################



Generation_observation = function(beta, Q, C, vit, Delta, loc0 = c(0,0)){
  K = dim(theta)[2]
  nb_obs = length(Q)
  Obs <- matrix(NA, nb_obs, 2)
  Obs[1,] = loc0
  for (t in 1:nb_obs) {
    
    # On regarde l'état dans lequel on se trouve.
    q = Q[t]
    
    # On simule le déplacement.
    xy = simLangevinMM(beta[q], vit, c(Delta[t-1],Delta[t]), loc0, cov_list = C, keep_grad = FALSE)
    print(xy)
    Obs[t,] = as.vector(as.numeric(xy[2,][,1:2])) # On ne prend que les coordonnées du déplacement.
    loc0 = as.vector(as.numeric(xy[2,][,1:2]))
  }
  return(Obs)}







################################################################################
###                     Procédures forward et backward                       ###                               
################################################################################

forward_2.0 = function( A, B, PI){
  T = dim(B)[1]
  K = dim(B)[2]
  
  alp = matrix(1, ncol = K, nrow = T)
  somme = c()
  
  # On initialise.
  alp[1,] = PI * B[1,]
  
  for (t in 1:(T-1)){
    a = (alp[t,] %*% A) %*% diag(B[t+1,])
    somme = c(somme,sum(a))
    alp[t+1,] = a / sum(a)
  }
  return(list(alp,somme,sum(alp[T,])))
}
backward_2.0 = function( A, B){
  T = dim(B)[1]
  K = dim(B)[2]
  
  # Création de la matrice initialisée (dernière ligne égale à 1).
  bet = matrix(1, ncol = K, nrow = T)
  somme = c()
  
  for (t in (T-1):1){
    b = A %*% (B[t+1,] * bet[t+1,])
    somme = c(somme,sum(b))
    bet[t,] = b / sum(b)
  }
  return(list(bet,somme))
}







################################################################################
###                   Initialisation et Résultat optimal                     ###                               
################################################################################

Result_opt = function(obs, N_etats, C, Q, J){
  T = dim(Obs)[[1]]
  Coef = c()
  
  # On gère les différentes dimensions.
  if (dimension == 2){Z = data.frame(obs$Z1,obs$Z2)}
  else {Z = obs$Z}
  
  # On découpe Z selon les différents états. Je vais simplement stocker les 
  # indices et non pas directement les valeurs des Z_i. 
  
  # Répartition.
  Z_rep = list()
  model = list()
  for (k in 1:N_etats){
    Z_rep[k] = c(0)
  }
  
  # On ne parcourt qu'une seule fois les données et on les répartit dans les 
  # différentes listes selon l'état prédit pour le déplacement. 
  for (i in 1:T){
    Z_rep[[Q[i]]] = c(Z_rep[[Q[i]]],i)   # On ajoute la coordonnée dans la bonne liste.
  }
  
  # On construit les différentes sous-parties de Z (autant que d'états).  
  for (i in 1:N_etats){
    l1 = c()
    l2 = c()
    C_nv = c()
    lgr = length(Z_rep[[i]])
    for (t in Z_rep[[i]][2:lgr]){
      l1 = c(l1,Z[t,1])
      l2 = c(l2,Z[t,2])
      C_nv = c(C_nv,C[t,])
    }
    C_nv = matrix(c(C_nv,C_nv), nrow = 2*length(l1),byrow = TRUE)
    model = lm(c(l1,l2) ~ C_nv)
    coef = coef(model)[1:J+1]
    Coef = c(Coef,coef)
  }
  
  # On s'occupe maintenant d'estimer la matrice de transition. 
  mcFitMLE <- markovchainFit(data = Q)
  A = matrix(mcFitMLE$estimate[1:N_etats],N_etats,N_etats)
  
  return(list(A,matrix(Coef, ncol = N_etats)))
  
}

initialisation = function(obs, N_etats, C, J, methode = 'kmeans'){
  T = dim(Obs)[[1]]
  Coef = c()
  
  if (methode == 'kmeans'){
    # On gère les différentes dimensions.
    if (dimension == 2){Z = data.frame(obs$Z1,obs$Z2)}
    else {Z = obs$Z}
    print(Z)
    # On effectue le kmeans sur les observations.
    km = kmeans(Z, N_etats)
    
    # On extrait la suite prédite des états. 
    Q_km = km$cluster
    
    # On découpe Z selon les différents états. Je vais simplement stocker les 
    # indices et non pas directement les valeurs des Z_i. 
    
    # Répartition.
    Z_rep = list()
    model = list()
    for (k in 1:N_etats){
      Z_rep[k] = c(0)
    }
    
    # On ne parcourt qu'une seule fois les données et on les répartit dans les 
    # différentes listes selon l'état prédit pour le déplacement. 
    for (i in 1:T){
      Z_rep[[Q_km[i]]] = c(Z_rep[[Q_km[i]]],i)   # On ajoute la coordonnée dans la bonne liste.
    }
    
    # On construit les différentes sous-parties de Z (autant que d'états).  
    for (i in 1:N_etats){
      l1 = c()
      l2 = c()
      C_nv = c()
      lgr = length(Z_rep[[i]])
      for (t in Z_rep[[i]][2:lgr]){
        l1 = c(l1,Z[t,1])
        l2 = c(l2,Z[t,2])
        C_nv = c(C_nv,C[t,])
      }
      C_nv = matrix(c(C_nv,C_nv), nrow = 2*length(l1),byrow = TRUE)
      model = lm(c(l1,l2) ~ C_nv)
      coef = coef(model)[1:J+1]
      Coef = c(Coef,coef)
    }
  }
  
  # On s'occupe maintenant d'estimer la matrice de transition. 
  #print(matrix(Q_km,1,T))
  mcFitMLE <- markovchainFit(data = Q_km)
  A = matrix(mcFitMLE$estimate[1:N_etats],N_etats,N_etats)
  
  return(list(A,matrix(Coef, ncol = N_etats)))
}

initialisation_modif = function(obs, N_etats, C, J, methode = 'kmeans'){
  nbr_obs = dim(Obs)[[1]]
  Coef = c()
  
  if (methode == 'kmeans'){
    # On gère les différentes dimensions.
    if (dimension == 2){Z = data.frame(obs$Z1[1:(nbr_obs-2)],obs$Z2[1:(nbr_obs-2)])}
    else {Z = obs$Z}
    print(Z)
    # On effectue le kmeans sur les observations.
    km = kmeans(Z, N_etats)
    
    # On extrait la suite prédite des états. 
    Q_km = km$cluster
    
    # On découpe Z selon les différents états. Je vais simplement stocker les 
    # indices et non pas directement les valeurs des Z_i. 
    
    # Répartition.
    Z_rep = list()
    model = list()
    for (k in 1:N_etats){
      Z_rep[k] = c(0)
    }
    
    # On ne parcourt qu'une seule fois les données et on les répartit dans les 
    # différentes listes selon l'état prédit pour le déplacement. 
    for (i in 1:nbr_obs){
      Z_rep[[Q_km[i]]] = c(Z_rep[[Q_km[i]]],i)   # On ajoute la coordonnée dans la bonne liste.
    }
    
    # On construit les différentes sous-parties de Z (autant que d'états).  
    for (i in 1:N_etats){
      l1 = c()
      l2 = c()
      C_nv = c()
      lgr = length(Z_rep[[i]])
      for (t in Z_rep[[i]][2:lgr]){
        l1 = c(l1,Z[t,1])
        l2 = c(l2,Z[t,2])
        C_nv = c(C_nv,C[t,])
      }
      C_nv = matrix(c(C_nv,C_nv), nrow = 2*length(l1),byrow = TRUE)
      model = lm(c(l1,l2) ~ C_nv)
      coef = coef(model)[1:J+1]
      Coef = c(Coef,coef)
    }
  }
  
  # On s'occupe maintenant d'estimer la matrice de transition. 
  #print(matrix(Q_km,1,T))
  mcFitMLE <- markovchainFit(data = Q_km)
  A = matrix(mcFitMLE$estimate[1:N_etats],N_etats,N_etats)
  
  return(list(A,matrix(Coef, ncol = N_etats)))
}


#Init = initialisation(Obs, K, C, J)





################################################################################
###                               Viterbi                                    ###                               
################################################################################

retourner = function(liste){
  l = length(liste)
  liste_retournee = c()
  for (i in l:1){
    liste_retournee = c(liste_retournee, liste[i])
  }
  return(liste_retournee)
}

Viterbi = function(A,B,PI){
  T = dim(B)[1]
  K = dim(A)[1]
  
  # Création des matrices.
  delta = matrix(1,T,K)
  phi = matrix(0,T,K)
  
  # Initialisation de la matrice delta.
  delta[1,] = PI * B[1,]
  # Récurrence.
  for (t in 2:T){
    for (j in 1:K){
      
      # On calcule toutes les transitions possibles arrivant dans l'état j.
      dA = delta[t-1,] * A[,j] * 4
      
      # On trouve la transition la plus probable.
      Ind_max = which.max(dA)
      C_max = dA[Ind_max]
      
      # On modifie delta.
      delta[t,j] = C_max * B[t,j]
      
      # On modifie phi.
      phi[t,j] = Ind_max 
    }
  }
  
  # On construit maintenant la suite d'état la plus probable.
  
  # On prend l'état final le plus probable.
  fin_max = which.max(delta[T,])
  P_et = delta[T,fin_max]
  Q_et = c(fin_max)  # Initialisation de la suite d'état.
  for (t in (T-1):1){
    new_etat = phi[t+1,Q_et[T - t]]
    #print(new_etat)
    Q_et = c(Q_et, new_etat)
  }
  return(retourner(Q_et))
}



################################################################################
###                                 EM                                       ###                               
################################################################################

EM_Langevin_modif_A = function(obs, Lambda, delta, vit, C, G = 10, moyenne = FALSE){
  
  compteur = 0
  
  # On gère la dimension du modèle. 
  Dim = dim(obs)[2]
  if (Dim == 4){
    Z = cbind( obs$Z1, obs$Z2)
    dimension = 2
  }
  else {Z = obs$Z}
  
  # Extraction des paramètres du modèle.
  A = Lambda$A
  B = Lambda$B
  PI = Lambda$PI
  
  # On gère l'option moyenne si besoin.
  somme_theta = matrix(0,J,K)
  somme_A = matrix(0,K,K)
  
  
  # Construction de la matrice C avec la division temporelle.
  C_temp = C
  for (i in 1:T){
    C_temp[i] = C_temp[i]/sqrt(delta[i])
    C_temp[T + i] = C_temp[T+i]/sqrt(delta[i])
  }
  
  
  while (compteur < G){
    print(compteur)
    ### EXPECTATION.
    
    # GAMMA.
    F = forward_2.0( A, B, PI)
    alp = F[[1]]
    S_alp = F[[2]]
    proba_obs = F[[3]]
    
    #print(F)
    
    Back = backward_2.0( A, B)
    bet = Back[[1]]
    S_bet = Back[[2]]
    gam = alp * bet
    
    for (t in 1:T){gam[t,] = gam[t,]/(sum(gam[t,]))}
    
    #print(gam)
    
    # On gère la potentiel présence de NA dans la matrice gam.
    if (any(is.na(gam))){
      warning("Il y a présence d'au moins un NA dans la matrice gam, voici le dernier résultat")
      if (moyenne){return(list(somme_A/G,somme_theta/G,sqrt(vit)))}
      else {return(list(A,theta_nv,sqrt(vit)))}
    }
    
    
    ## CALCUL DE A.
    
    Xi = array(1,dim = c(K,K,T),)
    for (t in 1:(T-1)){
      #print(diag(gam[t,] * (1/bet[t,])) %*% A %*% diag(B[t+1,] * bet[t+1,]))
      Xi[,,t] = diag(gam[t,] * (1/bet[t,])) %*% A %*% diag(B[t+1,] * bet[t+1,])
    }
    
    # Je fais la somme des Xi pour t allant de 1 à T-1.
    somme_Xi = matrix(0,K,K)
    for (t in 1:(T-1)){somme_Xi = somme_Xi + Xi[,,t]}
    
    # Je fais la somme des gamma pour t allant de 1 à T-1.
    somme_gam = matrix(0, nrow = K, ncol = K, byrow = TRUE)
    for (k in 1:K){
      sg = sum(gam[1:T-1,k])
      #print(matrix(sg, nrow= 1, ncol = K, byrow = TRUE))
      somme_gam[k,] = matrix(sg, nrow= 1, ncol = K, byrow = TRUE)
    }
    
    
    # On obtient l'estimateur de la matrice A.
    A = format_mat(somme_Xi * somme_gam)
    somme_A = somme_A + A
    
    
    # THETA.
    theta_nv = matrix(1,J,K)
    for (k in 1:K){
      
      # On gère les deux cas différents selon la dimension.
      if (dimension == 2){
        model = lm(c(Z[,1],Z[,2]) ~ C_temp, weights= c(gam[,k],gam[,k]))
      }
      else {
        model = lm(Z ~ C, weights=gam[,k])
      }
      
      # On récupère les coefficients.
      vit = summary(model)$sigma
      theta_nv[,k] = coef(model)[2:(J+1)]
    }
    #print(theta_nv)
    # On gère la potentielle moyenne à calculer.
    somme_theta = somme_theta + theta_nv
    
    # On met à jour la matrice des probabilités des émissions.
    B = proba_emission(obs, C, theta_nv, delta, vit, dimension)
    # On met à jour le compteur.
    compteur = compteur + 1
  }
  
  # On gère la moyenne si nécessaire.
  if (moyenne){return(list(somme_A/G,somme_theta/G,sqrt(vit)))}
  else {return(list(A,theta_nv,sqrt(vit)))}
}














################################################################################
#                                                                              #
#                         Simulation des observations                          #
#                                                                              #
################################################################################       

nbr_obs = 1000      
K = 3       
J = 2        
dimension = 2  
vit = 0.4            
PI = c(.5,.3,.2)               


# Creation de la suite des instants.
tps_final = 1000
ano = 50 # nombre d'anomalies.

instants = seq(1, tps_final, length.out = (nbr_obs + ano))
anomalies = sample(1:(nbr_obs + ano),ano)
tps = instants[-anomalies]
incr = increments(tps)

# Creation de la liste des covariables via Rhabit.

liste_cov = list()
for (i in 1:J){
  liste_cov[[i]] = simSpatialCov(lim, nu, rho, sigma2, resol = resol,
                                 mean_function = mean_function,
                                 raster_like = t)
}

# Creation de la suite des etats caches.

A = matrix(c(.85,.05,.1,.03,.91,.06,.03,.07,.9),
           ncol = K,
           nrow = K,
           byrow = TRUE)

Q = CM_generateur( A, nbr_obs)

# Le parametre de lien. 

theta = Nu(BETA(K,J), 1, vit)
theta

# Simulation des observations en utilisant Rhabit. 

Obs = Generation_observation(beta = matrix_to_list(theta), 
                             Q, C = liste_cov, vit, tps)

# On calcule les valeurs du gradient des covariables en les observations et 
# les met sous le bon format. 

MatObs = matrix(c(Obs$X1,Obs$X2),nbr_obs -1 ,2,byrow = FALSE)
CovGradLocs = covGradAtLocs(MatObs, liste_cov)
C = matrix(NA, 2*(nbr_obs-1), J)
for (t in 1:(t-1)){
  for (j in 1:J){
    print(t)
    C[t,j] = CovGradLocs[t, j, 1]
    C[t - 1 + nbr_obs,j] = CovGradLocs[t, j, 2]
  }
}



# On initialise les parametres. 

Init = initialisation_modif(Obs, K, C, J)
A_init = Init[[1]]; theta_init = Init[[2]]


vit_initiale = 0.2  
theta_initial = Nu(BETA(K,J), 1, vit_initiale)   

Lambda = list('A' = A,
              'B' = proba_emission(Obs, C, theta, Delta_erreurs,  vit_initiale, 
                                   dimension),
              'PI' = PI)


E = EM_Langevin_modif_A( Obs, Lambda, incr, vit, C, G = 20, moyenne = FALSE)
print(list(E[[1]],E[[2]],E[[3]]))
theta



