library(flexmix)


################################################################################
###                     Probabilités des apparitions                                                       
### Fonction calculant la probabilité des émissions, c'est à dire la probabilité 
### de chaque "partie" du déplacement de l'animal. Il prend en compte la 
### dimension 1 et 2.
###
### increments : dataframe contenant : 
###     $t : les instants où les observations ont lieu.
###     $grad_ci_x/grad_ci_y : gradient de la covariable i selon x ou y. 
###     $etats_caches : la suite des etats caches ayant servi a la simulation.
###     $delta_t : le pas de temps. 
###     $deplacement : deplacement selon une direction (x ou y) qui est precisee
###                    dans la variable suivante.
###     $dimension : variable precisant si le deplacement associe se fait selon
###                  x ou y. 
###     $cov.i : valeur du gradient après la division par deux fois le pas de 
###              temps.
### param : liste de contenant autant de liste qu'il y a d etats caches 
###          differents. Pour chaque etat, on a : 
###     $nu : vecteur contenant les coefficients beta associes a l etat. 
###     $vitesse : valeur de la vitesse (gamma) du processus dans cet etat. 
### 
### Renvoie une matrice avec autant de lignes qu'il y a de déplacements et avec
### autant de colonne qu'il y a d'etats caches differents. 
################################################################################

proba_emission = function(increments, param){
  K = length(param)
  cov_index <- str_detect(colnames(increments), "cov")
  J = sum(cov_index)
  # Création de la matrice.
  B = matrix(0, ncol = K, nrow = nrow(increments)/2)
  for (k in 1:K){
    prov <- matrix(dnorm( increments$deplacement, 
                          mean = as.matrix(increments[, cov_index]) %*%
                            matrix(param[[k]]$nu, nrow=J),
                          sd = param[[k]]$vitesse), ncol= 2) 
    B[,k] <- prov[,1]* prov[,2]    
  }
  return(B)
}

################################################################################
###                       Procedures forward et backward                                
### Procedure Forward.
### A : matrice de transition.
### B : matrice des probabilités d emission sous le modele.
### PI : vecteur des probabilites initiales des etats.
###
### La fonction forward2.0 applique la procedure forward au modele. Elle renvoie
### une matrice contenant les probabilites. La matrice a nbr_obs lignes et 
### autant de colonnes qu il y a d etats.
###
### Procedure Backward.
### A : matrice de transition.
### B : matrice des probabilités d emission sous le modele.
###
### La fonction backward2.0 applique la procedure backward au modele. Elle 
### renvoie une matrice contenant les probabilites. La matrice a nbr_obs lignes  
### et autant de colonnes qu il y a d etats.
###
################################################################################

forward_2.0 = function( A, B, PI){
  
  nbr_obs = dim(B)[1]
  K = dim(B)[2]
  alp = matrix(1, ncol = K, nrow = nbr_obs)
  alp_norm <-  rep(1,nbr_obs)
  # On initialise.
  alp[1,] = PI * B[1,]
  alp_norm[1] = sum(alp[1,])
  alp[1,] = alp[1,] / alp_norm[1]
  
  for (t in 1:(nbr_obs-1)){
    a = (alp[t,] %*% A) %*% diag(B[t+1,])
    alp_norm[t+1] = sum(a)
    alp[t+1,] = a / alp_norm[t+1]
  }
  return(list(alpha = alp, alpha_norm= alp_norm))
}
# alp = forward_2.0( Lambda$A, Lambda$B, Lambda$PI)


backward_2.0 = function( A, B){
  nbr_obs = dim(B)[1]
  K = dim(B)[2]
  beta_norm = rep(1, nbr_obs)
  # Création de la matrice initialisée (dernière ligne égale à 1).
  bet = matrix(1, ncol = K, nrow = nbr_obs)
  
  # On divise par la somme de la première ligne.
  beta_norm[1] = sum(bet[nbr_obs,])
  bet[nbr_obs,] = bet[nbr_obs,]/beta_norm[1]
  
  for (t in (nbr_obs-1):1){
    b = A %*% (B[t+1,] * bet[t+1,])
    beta_norm[t] = sum(b)
    bet[t,] = b / beta_norm[t]
  }
  return(list(beta= bet, neta_norm = beta_norm))
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

################################################################################
###         Estimation des paramètres selon une suite d'etats connue
###
### increments : dataframe contenant : 
###     $t : les instants où les observations ont lieu.
###     $grad_ci_x/grad_ci_y : gradient de la covariable i selon x ou y. 
###     $etats_caches : la suite des etats caches ayant servi a la simulation.
###     $delta_t : le pas de temps. 
###     $deplacement : deplacement selon une direction (x ou y) qui est precisee
###                    dans la variable suivante.
###     $dimension : variable precisant si le deplacement associe se fait selon
###                  x ou y. 
###     $cov.i : valeur du gradient après la division par deux fois le pas de 
###              temps.
### etats : suite des etats supposes. 
###
### Renvoie une liste contenant autant de liste qu'il y a d etats caches 
### differents. La liste associee a chaque etat contient :
###     $nu : vecteur des coefficients associes a l etat.
###     $vitesse : 
###
################################################################################

estim_etatsconnus = function(increments, etats){
  Y=increments$deplacement 
  cov_index <- str_detect(colnames(increments), 'cov')
  Z = increments[, cov_index]
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

################################################################################
###                    Initialisation des parametres 
###
### increments : dataframe contenant : 
###     $t : les instants où les observations ont lieu.
###     $grad_ci_x/grad_ci_y : gradient de la covariable i selon x ou y. 
###     $etats_caches : la suite des etats caches ayant servi a la simulation.
###     $delta_t : le pas de temps. 
###     $deplacement : deplacement selon une direction (x ou y) qui est precisee
###                    dans la variable suivante.
###     $dimension : variable precisant si le deplacement associe se fait selon
###                  x ou y. 
###     $cov.i : valeur du gradient après la division par deux fois le pas de 
###              temps.
### K : nombre d etats differents.
### methode : str representant le nom d une methode. Les choix sont : 
###     - 'kmeans' : methode des kmeans.
### 
### Renvoie une liste contenant : 
###     $A : la matrice de transition prédite. 
###     $param : appel de la fonction estim_etatsconnus avec comme suite d etats 
###              la suite predite par la methode. 
###
################################################################################

initialisation2.0 = function(increments, K){
  ## create formula
  cov_list <- str_extract(colnames(increments), pattern = "cov.[1234567890]?") %>% na.omit()
  rhs <- paste0( c("-1", cov_list), collapse= " + ")
  form <- paste0( "deplacement ~ ", rhs)
  # fit regression mixture
  m1 <- flexmix( formula = as.formula(form), k= K, data= increments)
  coef_est <- parameters(m1)
  
  #format parameter
  param <- lapply(1:K, function(k_){
    list(nu = coef_est[1:J,k_],
         vitesse = coef_est[J+1, k_],
         beta = coef_est[1:J,k_] * coef_est[J+1, k_])
  })
  
  

  A = m1@cluster %>% 
    as_tibble() %>% 
    rename(P2= value) %>% 
    mutate(P1=lag(P2)) %>% 
    na.omit() %>% 
    group_by(P1,P2) %>% 
    summarise(n =n()) %>% 
    mutate(ntot=sum(n)) %>% 
    mutate(p = n/ntot) %>%
    ungroup() %>% 
    dplyr::select(P1,P2,p) %>% 
    pivot_wider( names_from = P2, values_from = p ) %>% as.matrix()

  PI = 1*(c(1:K)==m1@cluster[1])
  return(list(A = A, param = param, PI = PI))
}

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
### increments : dataframe contenant : 
###     $t : les instants où les observations ont lieu.
###     $grad_ci_x/grad_ci_y : gradient de la covariable i selon x ou y. 
###     $etats_caches : la suite des etats caches ayant servi a la simulation.
###     $delta_t : le pas de temps. 
###     $deplacement : deplacement selon une direction (x ou y) qui est precisee
###                    dans la variable suivante.
###     $dimension : variable precisant si le deplacement associe se fait selon
###                  x ou y. 
###     $cov.i : valeur du gradient après la division par deux fois le pas de 
###              temps.
### Lambda : liste à trois elements contenant les paramètres du modèle initial. 
###     $A : la matrice de transition initiale.
###     $B : la matrice de probabilité d'émission initiale. 
###     $PI : Vecteur contenant les probabilités initiales de commencer dans 
###           chaque état.
### Vitesses : liste des vitesses initiales pour chaque etat. 
### G : Nombre d iteration de l'EM (par défaut 20).
### moyenne : booléen (par défaut FALSE) régissant l'utilisation d une moyenne
###           ou non.
### 
### La fonction renvoie une liste : 
###       $A : matrice de transition des etats caches. 
###       $Nu : paramètre Nu des etats caches.
###       $Vitesses : liste des vitesses de chaque etat cache.
################################################################################

EM_Langevin = function(increments, Lambda, G = 10, moyenne = FALSE){
  compteur = 0
  # On gère la dimension du modèle.
  nbr_obs = dim(increments)[1]/2
  
  # Extraction des paramètres du modèle.
  A = Lambda$A
  B = Lambda$B
  PI = Lambda$PI
  K = dim(B)[2]
  
  # On gère l'option moyenne si besoin.
  somme_theta = matrix(0,J,K)
  somme_A = matrix(0,K,K)
  
  # Extraction de la matrice C.
  C = as.matrix(increments[,stringr::str_detect(colnames(increments),"cov")])
  J = dim(C)[2]
  
  while (compteur < G){
    print(paste('Tour',compteur))
    ### EXPECTATION.
    
    # GAMMA.
    alpha_comp = forward_2.0( A, B, PI)
    alp = alpha_comp$alpha
    beta_comp = backward_2.0( A, B)
    bet = beta_comp$beta
    
    gam = alp * bet
    for (t in 1:dim(gam)[1]){gam[t,] = gam[t,]/(sum(gam[t,]))}
    #print(gam)
    
    
    # On gère la potentiel présence de NA dans la matrice gam.
    if (any(is.na(gam))){  # anyNA marche aussi je crois bien.
      warning("Il y a présence d'au moins un NA dans la matrice gam, voici le dernier résultat")
      if (moyenne){return(list(somme_A/G,somme_theta/G,sqrt(vit)))}
      else {return(list(A,theta_nv,sqrt(vit)))}
    }
    
    ## CALCUL DE A.
    
    Xi = array(1,dim = c( K, K, nbr_obs-1),)
    for (t in 1:(nbr_obs-2)){
      Xi[,,t] = diag(alp[t,]) %*% A %*% diag(B[t+1,] * bet[t+1,])
      Xi[,,t] = Xi[,,t]/sum(Xi[,,t])
    }
    
    # Je fais la somme des Xi pour t allant de 1 à T-1.
    somme_Xi = matrix(0,K,K)
    for (t in 1:(nbr_obs-1)){
      somme_Xi = somme_Xi + Xi[,,t]
      }
    
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
    Params = lapply(1:K, function(k){
      model = lm(increments$deplacement ~ -1 + C, weights = c(gam[,k],gam[,k]))
      return(list(nu = coef(model), vitesse = summary(model)$sigma))
    })
    
    # On gère la potentielle moyenne à calculer.
    somme_theta = somme_theta + theta_nv
    
    # On met à jour la matrice des probabilités des émissions.
    B = proba_emission(increments, Params)
    
    Aff = AffParams(Params)
    nu_nv = Aff$nu
    vit_nv = Aff$vitesses
    print(A)
    print(nu_nv)
    print(vit_nv)
    # On met à jour le compteur.
    compteur = compteur + 1
   } 
  
  
  
    
  # On gère la moyenne si nécessaire.
  if (moyenne){return(list(A = somme_A/G,
                           Nu = somme_theta/G,
                           Vitesses = sqrt(Vits)))} else{return(list(A = A, 
                                                                     Nu = nu_nv, 
                                                                     Vitesses = vit_nv))}
}
# E = EM_Langevin(increments_dta, Lambda, Vc(0.4,0.4), G = 10)
# E
# 
# 
# etats_forts = function(B){
#   l = nrow(B)
#   liste_etats = numeric(l)
#   for (t in 1:l){
#     if (B[t,1] > B[t,2]){
#       liste_etats[t] = 1
#     } else if (B[t,1] < B[t,2]){
#         liste_etats[t] = 2
#         }
#   }
#   return(liste_etats)
# }
# Q_test = etats_forts(Lambda$B)
# etats_caches == Q_test
# 


