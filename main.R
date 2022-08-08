RFoptions(install="no")
source('utils.R')
source('simu.R')
source('estim.R')


nbr_obs = 1000      
K = 2       
J = 2        
dimension = 2  
vit = 0.4   ### C est gamma et pas gamma2 !!         
#PI = c(.5,.3,.2)    
PI = c(.5,0.5)

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

# liste_cov_2 = lapply(1:J, function(j){Rhabit::simSpatialCov(lim, nu, rho, sigma2,
#                                                resol = resol,
#                                                mean_function = mean_function,
#                                                raster_like = TRUE)})

# Creation de la suite des etats caches.

# A = matrix(c(.85,.05,.1,.03,.91,.06,.03,.07,.9),
#            ncol = K,
#            nrow = K,
#            byrow = TRUE)
A = matrix(c(0.85,.15,.09,0.91),
           ncol = K,
           nrow = K,
           byrow = TRUE)

Q = CM_generateur( A, nbr_obs)

# Le parametre de lien. 

bets = BETA(K,J)
theta = Nu(bets, vit)
theta

# Simulation des observations en utilisant Rhabit. 

Obs = Generation_observation3.0(liste_theta = matrix_to_list(bets), Q,  liste_cov, 
                                Vits = c(.01,.01), tps)

Obs$Q = Q
Obs$t = tps
colnames(Obs) = c('X1','X2','Z1','Z2',
                  'grad_c1_x','grad_c1_y','grad_c2_x','grad_c2_y','etats','t')

# On construit le vecteur Y.
Y = c(Obs$Z1[1:nbr_obs-1]/sqrt(incr),Obs$Z2[1:nbr_obs-1]/sqrt(incr))

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

# # Construction de la matrice C avec la division temporelle.
# for (i in 1:(nbr_obs-1)){
#   C[i] = C[i]/sqrt(incr[i])
#   C[nbr_obs - 1 + i] = C[nbr_obs - 1 +i]/sqrt(incr[i])
# }


# On initialise les parametres. 

Init = initialisation(Obs, K, C, J)
A_init = Init$A; Beta_init = Init$Beta; Vits_init = Init$Vitesses


Res_opt = estim_etatsconnus(Y, incr, Q[1:nbr_obs-1], C)
Res_opt






################################################################################
#
#   Comparaison de la méthode Rhabit et de notre méthode pour 
#   analyser les données
#
################################################################################

# On utilise Rhabit pour analyser la trajectoire.
fitted_langevin <- fit_langevin_ud(
  cbind(X1,X2) ~  grad_c1 + grad_c2 + grad_c3  , data = Obs)

# Ca ne marche pas avec 3 covariables, on en prend donc que 2.
fitted_langevin <- fit_langevin_ud(
  cbind(X1,X2) ~  grad_c1 + grad_c2 , data = Obs)

# On sort alors les coefficients.
coef(fitted_langevin)
speed_coef(fitted_langevin)

# On utilise notre fonction du coup on dit qu'il n'y a qu'un seul etat. 
Q_1etat = replicate(nbr_obs-1,1)
Res_opt1 = estim_etatsconnus(Y, incr, Q_1etat, C[, 1:2])
Res_opt1







theta_initial = BetaToNu(Beta_init, Vits_init)
Lambda = list('A' = A,
              'B' = proba_emission(Obs, C, theta, incr,  c(0.4,0.4), 
                                   dimension),
              'PI' = PI)
#Lambda$B
E = EM_Langevin_modif_A( Obs, Lambda, incr, Vits_init, C, G = 5, moyenne = FALSE)
print(list(E[[1]],E[[2]],E[[3]]))
theta



