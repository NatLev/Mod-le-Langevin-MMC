################################################################################
###                     Procédures forward et backward                       ###                               
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
###                   Initialisation et Résultat optimal                     ###                               
################################################################################

Result_opt = function(obs, N_etats, C, Q, J){
  nbr_obs = dim(Obs)[[1]]
  Coef = c()
  Vitesses = numeric(N_etats)
  
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
  for (i in 1:nbr_obs){
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

initialisation = function(obs, N_etats, C, J, methode = 'kmeans', dimension = 2){
  nbr_obs = dim(Obs)[[1]]
  Coef = c()
  Vitesses = numeric(N_etats)
  
  if (methode == 'kmeans'){
    # On gère les différentes dimensions.
    if (dimension == 2){Z = data.frame('Z1' = obs$Z1[1:(nbr_obs-1)],
                                       'Z2' = obs$Z2[1:(nbr_obs-1)])}
    else {Z = obs$Z}
    # On effectue le kmeans sur les observations.
    km = kmeans(Z, N_etats)
    
    # On extrait la suite prédite des états. 
    Q_km = km$cluster
    
    # On découpe Z selon les différents états. Je vais simplement stocker les 
    # indices et non pas directement les valeurs des Z_i. 
    
    # Répartition.
    Z_rep = list()
    for (k in 1:N_etats){
      Z_rep[k] = c(0)
    }
    
    # On ne parcourt qu'une seule fois les données et on les répartit dans les 
    # différentes listes selon l'état prédit pour le déplacement. 
    for (i in 1:(nbr_obs-1)){
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
      vit = summary(model)$sigma
      coef = coef(model)[1:J+1]
      Coef[[i]] = coef
      Vitesses[i] = vit
    }
  }
  
  # On s'occupe maintenant d'estimer la matrice de transition. 
  #print(matrix(Q_km,1,T))
  mcFitMLE <- markovchainFit(data = Q_km)
  A = matrix(mcFitMLE$estimate[1:N_etats],N_etats,N_etats)
  
  return(list('A' = A, 'Beta' = Coef, 'Vitesses' = Vitesses))
}
initialisation(Obs, K, C, J
            )

initialisation_modif = function(obs, N_etats, C, J, methode = 'kmeans', dimension = 2){
  nbr_obs = dim(Obs)[[1]]
  Coef = c()
  Vitesses = numeric(N_etats)
  
  if (methode == 'kmeans'){
    # On gère les différentes dimensions.
    if (dimension == 2){Z = data.frame('Z1' = obs$Z1[1:(nbr_obs-1)],
                                       'Z2' = obs$Z2[1:(nbr_obs-1)])}
    else {Z = obs$Z}
    # On effectue le kmeans sur les observations.
    km = kmeans(Z, N_etats)
    
    # On extrait la suite prédite des états. 
    Q_km = km$cluster
    
    
    # On construit les différentes sous-parties de Z (autant que d'états).  
    for (i in 1:N_etats){
      indices = which(Q_km == i)  # On récupère les indices. 
      positions = as.matrix(obs[indices, c(1,2)]) # On récupère les positions.
      C_modif = as.matrix(C[ c( indices, indices + nbr_obs - 1), c(1:J)])
      
      model = lm(c(positions[,1], positions[,2]) ~ C_modif)
      Coef[[i]] = coef(model)[1:J+1]
      Vitesses[i] = summary(model)$sigma
    }
  }
  
  # On s'occupe maintenant d'estimer la matrice de transition. 
  #print(matrix(Q_km,1,T))
  mcFitMLE <- markovchainFit(data = Q_km)
  A = matrix(mcFitMLE$estimate[1:N_etats],N_etats,N_etats)
  
  return(list('A' = A, 'Beta' = Coef, 'Vitesses' = Vitesses))
}
initialisation_modif(Obs, K, C, J)
initialisation(Obs, K, C, J)



################################################################################
###                               Viterbi                                    ###                               
################################################################################


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
  if (moyenne){return(list(somme_A/G,somme_theta/G,sqrt(vit)))}
  else {return(list(A,theta_nv,sqrt(Vits)))}
}







