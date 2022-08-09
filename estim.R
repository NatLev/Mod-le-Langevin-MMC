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


proba_emission = function(increments, param){
  K = length(param)
  cov_index <- str_detect(colnames(increments), "cov")
  # Création de la matrice.
  B = matrix(0, ncol = K, nrow = nbr_obs)
  for (k in 1:K){
    prov <- matrix(dnorm( increments$deplacement, mean=  as.matrix(increments[, cov_index]) %*% matrix(param[[k]]$nu, nrow=2), sd = param[[k]]$vitesse), ncol= 2) 
    B[,k] <- prov[,1]* prov[,2]    
  }
  return(B)
}



################################################################################
###                       Procedures forward et backward                                
### A : matrice de transition.
### B : matrice des probabilités d emission sous le modele.
### PI : vecteur des probabilites initiales des etats.
###
### La fonction forward2.0 applique la procedure forward au modele. Elle renvoie
### une matrice contenant les probabilites. La matrice a nbr_obs lignes et 
### autant de colonnes qu il y a d etats.
################################################################################

forward_2.0 = function( A, B, PI){
  
  nbr_obs = dim(B)[1]
  K = dim(B)[2]
  alp = matrix(1, ncol = K, nrow = nbr_obs)
  
  # On initialise.
  alp[1,] = PI * B[1,]
  
  for (t in 1:(nbr_obs-1)){
    a = (alp[t,] %*% A) %*% diag(B[t+1,])
    alp[t+1,] = a / sum(a)
  }
  return(alp)
}
# alp = forward_2.0( Lambda$A, Lambda$B, Lambda$PI)


backward_2.0 = function( A, B){
  nbr_obs = dim(B)[1]
  K = dim(B)[2]
  
  # Création de la matrice initialisée (dernière ligne égale à 1).
  bet = matrix(1, ncol = K, nrow = nbr_obs)
  
  for (t in (nbr_obs-1):1){
    b = A %*% (B[t+1,] * bet[t+1,])
    bet[t,] = b / sum(b)
  }
  return(bet)
}


################################################################################
###                   Initialisation et Résultat optimal                     
###
### deplacements : vecteur 2*(n-1) lignes, les n-1 premieres deplacements en x,
###     et les n -1 suivants deplacements en y (divises par sqrt(accroisements))
### accroissements : vecteur (n-1) lignes, le temps entre deux acquisitions de 
###     positions.
### etats : vecteur (n-1) lignes, les états fixés.
### C : matrice (2(n-1), J) les n-1 premiere colonnes dérivées des covariables 
###     en x, les suivantes en y.
### 
### La fonction renvoie une liste d'autant de liste qu'il y a d'etats differents
### dans etats. Chaque sous liste est composee de $theta qui est le vecteur des 
### coefficients associes a l etat et $vitesse qui est la vitesse du processus 
### dans cet etat. 
## En tuilisant les notations du papier 
## Y les déplacements normalisés par sqrt(delta_t) 2n lignes
## Z lla matrice des gradients * 0.5 * multipliée par Tdelta,
## Etats la liste des etats
################################################################################

estim_etatsconnus = function(increments, etats){
  Y=increments_dta$deplacement 
  cov_index <- str_detect(colnames(increments, 'cov'))
  Z = increments_dta[, cov_index]
  K <-  n_distinct(etats)
  etats_2n <- c(etats, etats)
  
  estim_list <- lapply(1:K, function(k){
    index <- which(etats_2n == k)
    Y_k <- Y[index]
    Z_k <-  as.matrix(Z[index, ]  )
    res <- lm(Y_k ~ -1 + Z_k)
    return(list(nu = coef(res), vitesse= summary(res)$sigma, beta= coef(res) / summary(res)$sigma^2))

  })
  return(estim_list)
}


initialisation2.0 = function(increments, K,  methode = 'kmeans'){
  nbr_obs = nrow(increments)/2
  if (methode == 'kmeans'){
    Z <- matrix(increments$deplacement, nrow = nbr_obs)
    km = kmeans(Z, K) 
  }
  return(estim_etatsconnus(increments, km$cluster))}

#initialisation2.0(Obs, K, C, incr)

################################################################################
###                       Algorithme de Viterbi                                
### A : matrice de transition.
### B : matrice des probabilités d'emission sous le modele.
### PI : vecteur des probabilites initiales des etats.
###
### La fonction applique l'algorithme de Viterbi et renvoie la suite d etats 
### caches la plus probable selon le modele. La suite est renvoyee sous la forme
### d une liste. 
################################################################################





Viterbi = function(A,B,PI){
  nbr_obs = dim(B)[1]
  K = dim(A)[1]
  
  # Création des matrices.
  delta = matrix(1,nbr_obs,K)
  phi = matrix(0,nbr_obs,K)
  
  # Initialisation de la matrice delta.
  delta[1,] = PI * B[1,]
  # Récurrence.
  for (t in 2:nbr_obs){
    for (j in 1:K){
      
      # On calcule toutes les transitions possibles arrivant dans l'état j.
      dA = delta[t-1,] * A[,j] * 4  # Pourquoi 4 ???
      
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
  fin_max = which.max(delta[nbr_obs,])
  P_et = delta[nbr_obs,fin_max]
  Q_et = c(fin_max)  # Initialisation de la suite d'état.
  for (t in (nbr_obs-1):1){
    new_etat = phi[t+1,Q_et[nbr_obs - t]]
    #print(new_etat)
    Q_et = c(Q_et, new_etat)
  }
  return(retourner(Q_et))
}

##############################################################################
################################################################################
###                            Algorithme EM                     
###
### obs : dataframe de taille nbr_obs-1 avec les déplacements selon les deux 
###       variables rangés dans chaque colonne. Il n'y a pas de division par 
###       sqrt(incr) ici. 
### Lambda : liste à trois elements contenant les paramètres du modèle initial. 
###     $A : la matrice de transition initiale.
###     $B : la matrice de probabilité d'émission initiale. 
###     $PI : Vecteur contenant les probabilités initiales de commencer dans 
###           chaque état.
### delta : liste des pas de temps de taille nbr_obs-1.  
### vit : liste des vitesses initiales pour chaque etat. 
### G : Nombre d iteration de l'EM (par défaut 20).
### moyenne : booléen (par défaut FALSE) régissant l'utilisation d une moyenne
###           ou non.
### dimension : nombre de dimension (par défaut 2).
### 
### La fonction renvoie une liste : 
###       $A : matrice de transition des etats caches. 
###       $Nu : paramètre Nu des etats caches.
###       $Vitesses : liste des vitesses de chaque etat cache.
################################################################################

# Y = c(Obs$Z1[1:nbr_obs-1],Obs$Z2[1:nbr_obs-1])/sqrt(incr)
# EM_Langevin_modif_A = function(obs, Lambda, delta, vit, C, G = 10, 
#                                moyenne = FALSE, dimension = 2){
#   
#   compteur = 0
#   
#   # On gère la dimension du modèle. 
# 
#   if (dimension == 2){
#     Z = cbind( obs$Z1, obs$Z2)
#     dimension = 2
#   }
#   
#   else {Z = obs$Z}
#   nbr_obs = dim(Z)[1] - 1
#   
#   # Extraction des paramètres du modèle.
#   A = Lambda$A
#   B = Lambda$B[1:nbr_obs,]
#   PI = Lambda$PI
#   
#   # On gère l'option moyenne si besoin.
#   somme_theta = matrix(0,J,K)
#   somme_A = matrix(0,K,K)
#   
#   
#   # # Construction de la matrice C avec la division temporelle.
#   # C_temp = C
#   # for (i in 1:nbr_obs){
#   #   C_temp[i] = C_temp[i]/sqrt(delta[i])
#   #   C_temp[nbr_obs + i] = C_temp[nbr_obs+i]/sqrt(delta[i])
#   # }
#   
#   
#   while (compteur < G){
#     print(paste('Tour',compteur))
#     ### EXPECTATION.
#     
#     # GAMMA.
#     alp = forward_2.0( A, B, PI)
#     bet = backward_2.0( A, B)
#     
#     gam = alp * bet
#     for (t in 1:dim(gam)[1]){gam[t,] = gam[t,]/(sum(gam[t,]))}
#     print(gam)
#     
#     
#     # On gère la potentiel présence de NA dans la matrice gam.
#     if (any(is.na(gam))){
#       warning("Il y a présence d'au moins un NA dans la matrice gam, voici le dernier résultat")
#       if (moyenne){return(list(somme_A/G,somme_theta/G,sqrt(vit)))}
#       else {return(list(A,theta_nv,sqrt(vit)))}
#     }
#     
#     
#     ## CALCUL DE A.
#     
#     Xi = array(1,dim = c( K, K, nbr_obs-1),)
#     for (t in 1:(nbr_obs-2)){
#       Xi[,,t] = diag(gam[t,] * (1/bet[t,])) %*% A %*% diag(B[t+1,] * bet[t+1,])
#     }
#     
#     # Je fais la somme des Xi pour t allant de 1 à T-1.
#     somme_Xi = matrix(0,K,K)
#     for (t in 1:(nbr_obs-1)){somme_Xi = somme_Xi + Xi[,,t]}
#     
#     # Je fais la somme des gamma pour t allant de 1 à T-1.
#     somme_gam = matrix(0, nrow = K, ncol = K, byrow = TRUE)
#     for (k in 1:K){
#       sg = sum(gam[1:nbr_obs-1,k])
#       #print(matrix(sg, nrow= 1, ncol = K, byrow = TRUE))
#       somme_gam[k,] = matrix(1/sg, nrow= 1, ncol = K, byrow = TRUE)
#     }
#     
#     
#     # On obtient l'estimateur de la matrice A.
#     A = format_mat(somme_Xi * somme_gam)
#     somme_A = somme_A + A
#     
#     
#     # THETA.
#     theta_nv = matrix(1,J,K)
#     Vits = c()
#     print(length(c(gam[,k],gam[,k])))
#     for (k in 1:K){
#       # On gère les deux cas différents selon la dimension.
#       if (dimension == 2){
#         model = lm(c(Z[1:nbr_obs,1],Z[1:nbr_obs,2]) ~ C, weights= c(gam[,k],gam[,k]))
#         
#       }
#       else {
#         model = lm(Z ~ C, weights=gam[,k])
#       }
#       
#       # On récupère les coefficients.
#       Vits = c(Vits, summary(model)$sigma)
#       theta_nv[,k] = coef(model)[2:(J+1)]
#     }
#     # On gère la potentielle moyenne à calculer.
#     somme_theta = somme_theta + theta_nv
#     print(theta_nv)
#     # On met à jour la matrice des probabilités des émissions.
#     B = proba_emission(obs, C, theta_nv, delta, Vits)
#     #browser()
#     # On met à jour le compteur.
#     compteur = compteur + 1
#   }
#   
#   # On gère la moyenne si nécessaire.
#   if (moyenne){return(list(somme_A/G,somme_theta/G,sqrt(vit)))}
#   else {return(list(A,theta_nv,sqrt(Vits)))}
# }
# 
# 
# E = EM_Langevin_modif_A( Obs, Lambda, incr, Vits_init, C, G = 10, moyenne = FALSE)
# 
# 

Z = cbind(increments_dta$deplacement[1:l],increments_dta$deplacement[(l+1):(2*l)])


EM_Langevin_modif_A = function(obs, Lambda, delta, vit, C, G = 10, 
                               moyenne = FALSE, dimension = 2){
  
  compteur = 0
  
  # On gère la dimension du modèle. 

  if (dimension == 2){
    Z = cbind( obs$Z1, obs$Z2)
    dimension = 2
  }
  
  else {Z = obs$Z}
  nbr_obs = dim(Z)[1] - 1
  
  # Extraction des paramètres du modèle.
  A = Lambda$A
  B = Lambda$B[1:nbr_obs,]
  PI = Lambda$PI
  
  # On gère l'option moyenne si besoin.
  somme_theta = matrix(0,J,K)
  somme_A = matrix(0,K,K)
  
  
  # # Construction de la matrice C avec la division temporelle.
  # C_temp = C
  # for (i in 1:nbr_obs){
  #   C_temp[i] = C_temp[i]/sqrt(delta[i])
  #   C_temp[nbr_obs + i] = C_temp[nbr_obs+i]/sqrt(delta[i])
  # }
  
  
  while (compteur < G){
    print(paste('Tour',compteur))
    ### EXPECTATION.
    
    # GAMMA.
    alp = forward_2.0( A, B, PI)
    bet = backward_2.0( A, B)
    
    gam = alp * bet
    for (t in 1:dim(gam)[1]){gam[t,] = gam[t,]/(sum(gam[t,]))}
    print(gam)
    
    
    # On gère la potentiel présence de NA dans la matrice gam.
    if (any(is.na(gam))){
      warning("Il y a présence d'au moins un NA dans la matrice gam, voici le dernier résultat")
      if (moyenne){return(list(somme_A/G,somme_theta/G,sqrt(vit)))}
      else {return(list(A,theta_nv,sqrt(vit)))}
    }
    
    
    ## CALCUL DE A.
    
    Xi = array(1,dim = c( K, K, nbr_obs-1),)
    for (t in 1:(nbr_obs-2)){
      Xi[,,t] = diag(gam[t,] * (1/bet[t,])) %*% A %*% diag(B[t+1,] * bet[t+1,])
    }
    
    # Je fais la somme des Xi pour t allant de 1 à T-1.
    somme_Xi = matrix(0,K,K)
    for (t in 1:(nbr_obs-1)){somme_Xi = somme_Xi + Xi[,,t]}
    
    # Je fais la somme des gamma pour t allant de 1 à T-1.
    somme_gam = matrix(0, nrow = K, ncol = K, byrow = TRUE)
    for (k in 1:K){
      sg = sum(gam[1:nbr_obs-1,k])
      #print(matrix(sg, nrow= 1, ncol = K, byrow = TRUE))
      somme_gam[k,] = matrix(1/sg, nrow= 1, ncol = K, byrow = TRUE)
    }
    
    
    # On obtient l'estimateur de la matrice A.
    A = format_mat(somme_Xi * somme_gam)
    somme_A = somme_A + A
    
    
    # THETA.
    theta_nv = matrix(1,J,K)
    Vits = c()
    print(length(c(gam[,k],gam[,k])))
    for (k in 1:K){
      # On gère les deux cas différents selon la dimension.
      if (dimension == 2){
        model = lm(c(Z[1:nbr_obs,1],Z[1:nbr_obs,2]) ~ C, weights= c(gam[,k],gam[,k]))
        
      }
      else {
        model = lm(Z ~ C, weights=gam[,k])
      }
      
      # On récupère les coefficients.
      Vits = c(Vits, summary(model)$sigma)
      theta_nv[,k] = coef(model)[2:(J+1)]
    }
    # On gère la potentielle moyenne à calculer.
    somme_theta = somme_theta + theta_nv
    print(theta_nv)
    # On met à jour la matrice des probabilités des émissions.
    B = proba_emission(obs, C, theta_nv, delta, Vits)
    #browser()
    # On met à jour le compteur.
    compteur = compteur + 1
  }
  
  # On gère la moyenne si nécessaire.
  if (moyenne){return(list(A = somme_A/G,Nu = somme_theta/G,
                           Vitesses = sqrt(Vits)))}
  else {return(list(A = A, Nu = theta_nv, Vitesses = sqrt(Vits)))}
}


EM_Langevin = function(increments_dta, Lambda, vit, C, G = 10, moyenne = FALSE, 
                       dimension = 2){
  
  compteur = 0
  
  # On gère la dimension du modèle. 
  
  if (dimension == 2){
    Z = cbind( obs$Z1, obs$Z2)
    Y = c(obs$Z1[1:nbr_obs-1],obs$Z2[1:nbr_obs-1])/sqrt(incr)
  }
  
  else {Z = obs$Z}
  nbr_obs = dim(Z)[1] - 1
  
  # Extraction des paramètres du modèle.
  A = Lambda$A
  B = Lambda$B[1:nbr_obs,]
  PI = Lambda$PI
  
  # On gère l'option moyenne si besoin.
  somme_theta = matrix(0,J,K)
  somme_A = matrix(0,K,K)
  
  # Construction de la matrice C avec la division temporelle.
  C_temp = C
  for (i in 1:nbr_obs){
    C_temp[i] = C_temp[i]/sqrt(delta[i])
    C_temp[nbr_obs + i] = C_temp[nbr_obs+i]/sqrt(delta[i])
  }

  while (compteur < G){
    print(paste('Tour',compteur))
    ### EXPECTATION.
    
    # GAMMA.
    alp = forward_2.0( A, B, PI)
    bet = backward_2.0( A, B)
    
    gam = alp * bet
    for (t in 1:dim(gam)[1]){gam[t,] = gam[t,]/(sum(gam[t,]))}
    print(gam)
    
    
    # On gère la potentiel présence de NA dans la matrice gam.
    if (any(is.na(gam))){  # anyNA marche aussi je crois bien.
      warning("Il y a présence d'au moins un NA dans la matrice gam, voici le dernier résultat")
      if (moyenne){return(list(somme_A/G,somme_theta/G,sqrt(vit)))}
      else {return(list(A,theta_nv,sqrt(vit)))}
    }
    
    ## CALCUL DE A.
    
    # # Formule Bilmes.
    # Xi = array(1,dim = c( K, K, nbr_obs-1),)
    # for (t in 1:(nbr_obs-2)){
    #   Xi[,,t] = diag(gam[t,] * (1/bet[t,])) %*% A %*% diag(B[t+1,] * bet[t+1,])
    # }
    
    # Formule Rabiner. 
    Xi = array(1,dim = c( K, K, nbr_obs-1),)
    for (t in 1:(nbr_obs-2)){
      Xi[,,t] = diag(gam[t,]) %*% A %*% diag(B[t+1,] * bet[t+1,])
      Xi[,,t] = Xi[,,t]/sum(Xi[,,t])
    }
    
    # Je fais la somme des Xi pour t allant de 1 à T-1.
    somme_Xi = matrix(0,K,K)
    for (t in 1:(nbr_obs-1)){somme_Xi = somme_Xi + Xi[,,t]}
    
    # Je fais la somme des gamma pour t allant de 1 à T-1.
    somme_gam = matrix(0, nrow = K, ncol = K, byrow = TRUE)
    for (k in 1:K){
      sg = sum(gam[1:nbr_obs-1,k])
      somme_gam[k,] = matrix(1/sg, nrow= 1, ncol = K, byrow = TRUE)
    }
    
    
    # On obtient l'estimateur de la matrice A.
    A = format_mat(somme_Xi * somme_gam)
    somme_A = somme_A + A
    
    
    # THETA.
    theta_nv = matrix(1,J,K)
    Vits = c()
    for (k in 1:K){
      # On gère les deux cas différents selon la dimension.
      if (dimension == 2){
        model = lm(Y ~ C_temp, weights= c(gam[,k],gam[,k]))
      }
      else {
        model = lm(Z ~ C, weights=gam[,k])
      }
      # On récupère les coefficients.
      Vits = c(Vits, summary(model)$sigma)
      theta_nv[,k] = coef(model)[2:(J+1)]
    }
    # On gère la potentielle moyenne à calculer.
    somme_theta = somme_theta + theta_nv
    print(theta_nv)
    
    # On met à jour la matrice des probabilités des émissions.
    B = proba_emission(Z, C, theta_nv, delta, Vits)
    
    # On met à jour le compteur.
    compteur = compteur + 1
  }
  
  # On gère la moyenne si nécessaire.
  if (moyenne){return(list(A = somme_A/G,Nu = somme_theta/G,
                           Vitesses = sqrt(Vits)))}
  else {return(list(A = A, Nu = theta_nv, Vitesses = sqrt(Vits)))}
}

# E = EM_Langevin( Obs, Lambda, incr, c(0.4,0.4), C, G = 10, moyenne = FALSE)
# print(list(E[[1]],E[[2]],E[[3]]))
# theta




