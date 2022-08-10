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
PI = c(.6,0.4)

# Paramètre de création des covariables. 
set.seed(1)
lim <- c(-15, 15, -15, 15) # limits of map
resol <- 0.1 # grid resolution
rho <- 4; nu <- 1.5; sigma2 <- 10# Matern covariance parameters
mean_function <- function(z){# mean function
  -log(3 + sum(z^2))}



# Creation de la suite des instants.
tps_final = 100
ano = round(tps_final/100*5) # nombre d'anomalies.

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

liste_cov_2 = lapply(1:J, function(j){Rhabit::simSpatialCov(lim, nu, rho, sigma2,
                                               resol = resol,
                                               mean_function = mean_function,
                                               raster_like = TRUE)})




# Creation de la suite des etats caches.

# A = matrix(c(.85,.05,.1,.03,.91,.06,.03,.07,.9),
#            ncol = K,
#            nrow = K,
#            byrow = TRUE)
A = matrix(c(0.95,0.05,0.1,0.9),
           ncol = K,
           nrow = K,
           byrow = TRUE)

etats_caches = CM_generateur( A, nbr_obs)

# Le parametre de lien. 


beta_sim <-BETA(K,J) 
theta = Nu(beta_sim, vit)
theta

# Simulation des observations en utilisant Rhabit. 


Obs = Generation_observation3.0(liste_theta = matrix_to_list(beta_sim), etats_caches,  liste_cov, 
                                Vits = c(vit,vit), tps) %>% 
  rowid_to_column()

Obs
# preparation du jeu de données incréments

increments_dta <- Obs %>% 
  mutate(etats_caches = lag(etats_caches)) %>% 
  mutate(delta_t = t- lag(t)) %>% 
  mutate_at(vars(matches("grad")), ~lag(.x))  %>% 
  mutate_at(vars(matches("Z")), ~ .x/sqrt(delta_t)) %>% 
  rename(dep_x = Z1, dep_y = Z2) %>% 
  dplyr::select(-x, -y) %>% 
  na.omit() %>%  
  pivot_longer(matches("dep"), names_to = "dimension", values_to = "deplacement") %>% 
  arrange(dimension) %>% 
  create_covariate_columns()



## Y=increments_dta$deplacement
## Z = increments_dta[, c('cov1', 'cov2' ....)]

Init = initialisation2.0(increments_dta, K,  methode = 'kmeans')
  
A_init = Init$A; param_init = Init$param;

theta_opt = estim_etatsconnus(increments_dta, etats = increments_dta$etats_caches)



Lambda = list('A' = A,
              'B' = proba_emission(increments = increments_dta, param = param_init),
              'PI' = PI)
#Lambda$B

E = EM_Langevin(increments = increments_dta, 
                Lambda = Lambda, 
                Vitesses = c(0.4,0.4),
                G = 10, 
                moyenne = FALSE)
print(list(E[[1]],E[[2]],E[[3]]))
theta



