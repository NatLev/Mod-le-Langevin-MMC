RFoptions(install="no")
source('utils.R')
source('simu.R')
source('estim.R')


nbr_obs = 1000      
K = 2       
J = 3        
dimension = 2  
vit = 0.4            
#PI = c(.5,.3,.2)    
PI = c(.7,.3)

# Paramètre de création des covariables. 
set.seed(1)
lim <- c(-15, 15, -15, 15) # limits of map
resol <- 0.1 # grid resolution
rho <- 4; nu <- 1.5; sigma2 <- 10# Matern covariance parameters
mean_function <- function(z){# mean function
  -log(3 + sum(z^2))}



# Creation de la suite des instants.
tps_final = 1000
ano = 50 # nombre d'anomalies.

instants = seq(1, tps_final, length.out = (nbr_obs + ano))
anomalies = sample(1:(nbr_obs + ano),ano)
tps = instants[-anomalies]
incr = diff(tps)

# Creation de la liste des covariables via Rhabit.

liste_cov = list()
for (i in 1:J){
  liste_cov[[i]] = Rhabit::simSpatialCov(lim, nu, rho, sigma2, resol = resol,
                                 mean_function = mean_function,
                                 raster_like = TRUE)
}

# liste_cov_2 = lapply(1:J, Rhabit::simSpatialCov(lim, nu, rho, sigma2, 
#                                                resol = resol,
#                                                mean_function = mean_function,
#                                                raster_like = TRUE))

# Creation de la suite des etats caches.

# A = matrix(c(.85,.05,.1,.03,.91,.06,.03,.07,.9),
#            ncol = K,
#            nrow = K,
#            byrow = TRUE)
A = matrix(c(.85,.15,.09,.91),
           ncol = K,
           nrow = K,
           byrow = TRUE)

Q = CM_generateur( A, nbr_obs)

# Le parametre de lien. 

theta = Nu(BETA(K,J), vit)
theta

# Simulation des observations en utilisant Rhabit. 

Obs = Generation_observation3.0(liste_theta = matrix_to_list(theta), Q,  liste_cov, 
                                Vits = c(.01,.01), tps)

# On construit le vecteur Y.
Y = c(Obs$Z1[1:nbr_obs-1],Obs$Z2[1:nbr_obs-1])/sqrt(incr)

# On calcule les valeurs du gradient des covariables en les observations et 
# les met sous le bon format. 

# as.matrix(Obs[, c(1,2)])
MatObs = matrix(c(Obs$X1[1:(nbr_obs-1)],Obs$X2[1:(nbr_obs-1)]), (nbr_obs-1), 2, 
                byrow = FALSE)
CovGradLocs = covGradAtLocs(MatObs, liste_cov)
# dire à quoi sert C
C = matrix(NA, 2*(nbr_obs-1), J)
for (t in 1:(nbr_obs-1)){
  for (j in 1:J){
    C[t,j] = CovGradLocs[t, 1, j]
    C[t - 1 + nbr_obs,j] = CovGradLocs[t, 2, j]
  }
}

# Construction de la matrice C avec la division temporelle.
for (i in 1:(nbr_obs-1)){
  C[i] = C[i]/sqrt(incr[i])
  C[nbr_obs - 1 + i] = C[nbr_obs - 1 +i]/sqrt(incr[i])
}


# On initialise les parametres. 

Init = initialisation(Obs, K, C, J)
A_init = Init$A; Beta_init = Init$Beta; Vits_init = Init$Vitesses


Res_opt = estim_etatsconnus(Y, incr, Q[1:nbr_obs-1], C)



theta_initial = BetaToNu(Beta_init, Vits_init)
Lambda = list('A' = A,
              'B' = proba_emission(Obs, C, theta_initial, incr,  Vits_init, 
                                   dimension),
              'PI' = PI)
#Lambda$B
E = EM_Langevin_modif_A( Obs, Lambda, incr, Vits_init, C, G = 10, moyenne = FALSE)
print(list(E[[1]],E[[2]],E[[3]]))
theta



