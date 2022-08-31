RFoptions(install="no")
source('utils.R')
source('simu.R')
source('estim.R')


# # covariate generation ----------------------------------------------------
# 
# 
# liste_cov = lapply(1:J, function(j){Rhabit::simSpatialCov(lim, nu, rho, sigma2,
#                                                           resol = resol,
#                                                           mean_function = mean_function,
#                                                           raster_like = TRUE)})
# save(list = "liste_cov", file = "lise_cov_save.RData")
load("lise_cov_save.RData")


### Options graphiques
ggopts <- theme_light()+theme(text = element_text(size = 10),
                              axis.text = element_blank(),
                              axis.title = element_blank(),
                              legend.key.height= unit(3, "line"),
                              strip.background = element_rect(fill="white",colour = "grey") ,
                              strip.text = element_text(size = 25, color = "black"))
levels <- factor(paste("Covariate", 1:J), levels = paste("Covariate", 1:J))
cov_df <- do.call(rbind.data.frame,
                  lapply(1:J,
                         function(j){
                           Rhabit::rasterToDataFrame(liste_cov[[j]], levels[j])
                         }))


p1 <- ggplot(cov_df, aes(x,y)) + geom_raster(aes(fill = val)) +
  coord_equal() + scale_fill_viridis(name = "Value") + facet_wrap(level~.) +
  ggopts
p1



# parameters simu ---------------------------------------------------------


nbr_obs = 1000      
K = 2       
J = 2        
dimension = 2  
vit = 1 

# Creation de la suite des instants.
pdt = 0.1                  # Pas de temps de base.
ano = round(nbr_obs/100*5) # nombre d'anomalies.
tps = temps(pdt, nbr_obs, ano)

# tps_final = 30
# ano = round(nbr_obs/100*5) # nombre d'anomalies.
# 
# instants = seq(1, tps_final, length.out = (nbr_obs + ano))
# anomalies = sample(1:(nbr_obs + ano),ano)
# tps = instants[-anomalies]

#PI = c(.5,.3,.2)    
#PI = c(.6,0.4)

# Paramètre de création des covariables. 
#set.seed(1)
lim <- c(-30, 30, -30, 30) # limits of map
resol <- 0.1 # grid resolution
rho <- 2; nu <- 1.5; sigma2 <- 30# Matern covariance parameters
mean_function <- function(z){# mean function
  -log(3 + sum(z^2))}




# Creation de la liste des covariables via Rhabit.



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




# increments_dta %>% pivot_longer( matches("cov"),  names_to = "covariate", values_to = "grad") %>% ggplot() + 
#   geom_point(aes(x=grad, y=deplacement, col = as.factor(etats_caches))) +
#   facet_wrap(~covariate) + ggtitle("covariate")

# preparation du jeu de données incréments

# Generation(nbr_obs, pdt, A, liste_cov, theta)




B = proba_emission(increments = increments_dta, param = E$param)
V = Viterbi(E$A, B, E$PI)
V_flexmix = Viterbi(Relab$A, proba_emission(increments_dta, Relab$param), Relab$PI)

print("Réussite de Viterbi pour l'EM")
sum(Obs$etats_caches[1:(nbr_obs-1)] == V)/(nbr_obs-1)
print("Réussite de Viterbi pour Flexmix")
sum(Obs$etats_caches[1:(nbr_obs-1)] == V_flexmix)/(nbr_obs-1)

rel = relabel(A_init, param_init, PI_init)
rel$A; AffParams(rel$params)
E$A; AffParams(E$param)




################################################################################
#
#                              Tests version biaisee
#
################################################################################
A = matrix(c(0.95,0.05,0.1,0.9),
           ncol = K,
           nrow = K,
           byrow = TRUE)


N = 10
liste_theta = list(Nu(matrix(c(5, -5, -5, 5), ncol = K, nrow = J), vit),
               Nu(matrix(c(5, -5, -3, 3), ncol = K, nrow = J), vit),
               Nu(matrix(c(5, -5, -1, 1), ncol = K, nrow = J), vit))

Resultats = data.frame(col_inutile = 1:N)
Resultats_A = data.frame(col_inutile = 1:N)
Resultats_Nu = data.frame(col_inutile = 1:N)

for (i in 1:length(liste_theta)){
  theta = liste_theta[[i]]
  Viterbi_EM = numeric(N)
  Viterbi_Flexmix = numeric(N)
  liste_norm_A_EM = numeric(N)
  liste_norm_A_Flexmix = numeric(N)
  liste_norm_nu_EM = numeric(N)
  liste_norm_nu_Flexmix = numeric(N)
  
  print(theta)
  
  compteur = 0
  while (compteur < N) {
    print(paste0('Tour ',compteur+1, ' pour Nu',i))
    Obs = Generation(nbr_obs, pdt, A, liste_cov, theta)$Obs
    
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
    
    m1 <- flexmix(deplacement ~ -1 + cov.1 + cov.2, k= 2, data= increments_dta)
    table(m1@cluster, increments_dta$etats_caches)
    parameters(m1)
    
    Init = initialisation2.0(increments = increments_dta, K = 2)
    A_init = Init$A[,c(2,3)]; param_init = Init$param; PI_init = Init$PI
    
    Relab = relabel(A_init, param_init, PI_init)
    
    Lambda = list('A' = Relab$A,
                  'param' = Relab$param,
                  'PI' = Relab$PI)
    
    E = EM_Langevin(increments = increments_dta, 
                    Lambda = Lambda, 
                    G = 20, 
                    moyenne = FALSE)
    B = proba_emission(increments = increments_dta, param = E$param)
    V = Viterbi(E$A, B, E$PI)
    V_flexmix = Viterbi(Relab$A, proba_emission(increments_dta, Relab$param), Relab$PI)
    
    Viterbi_EM[compteur+1] = sum(Obs$etats_caches[1:(nbr_obs-1)]==V)/(nbr_obs-1)
    Viterbi_Flexmix[compteur+1] = sum(Obs$etats_caches[1:(nbr_obs-1)]==V_flexmix)/(nbr_obs-1)
    
    liste_norm_A_EM[compteur] = norm(A - E$A)
    liste_norm_A_Flexmix[compteur] = norm(A - Relab$A)
    liste_norm_nu_EM[compteur] = norm(theta - AffParams(E$param)$nu)
    liste_norm_nu_Flexmix[compteur] = norm(theta - AffParams(Relab$param)$nu)
    
   
    compteur = compteur + 1
  }
  # Test de Viterbi.
  Resultats[paste0('Nu',i,' EM')] = Viterbi_EM  
  Resultats[paste0('Nu',i,' Flexmix')] = Viterbi_Flexmix  
  
  # Test des normes.
  # A.
  Resultats_A[paste0('A EM Nu',i)] = liste_norm_A_EM
  Resultats_A[paste0('A Flexmix Nu',i)] = liste_norm_A_Flexmix
  
  # Nu.
  Resultats_Nu[paste0('Nu EM Nu',i)] = liste_norm_nu_EM
  Resultats_Nu[paste0('Nu Flexmix Nu',i)] = liste_norm_nu_Flexmix
  
}
boxplot(Resultats[,-1], 
        col = c('orange','orange','red','red','pink','pink'),
        main = "Test de Viterbi pour l'EM et Flexmix \n 2000 observations, 100 répétitions, 3 nus, \n 20 iterations de l'EM"
        
        )
boxplot(Resultats_A[-1],
        main = "Test sur la prédiction des matrices de transitions \n 
        2000 obs, 100 gen, 20 tours d'EM trois Nu, pdt = 0.1",
        col = c('red','red','orange','orange','yellow','yellow'),
        ylab = 'Erreurs')

boxplot(Resultats_Nu[-1],
        main = "Test sur la prédiction des Nu \n 
        2000 obs, 100 gen, 20 tours d'EM trois Nu, pdt = 0.1",
        col = c('red','red','orange','orange','yellow','yellow'),
        ylab = 'Erreurs')

################################################################################
#
#                              Tests version non biaisee
#
################################################################################
nbr_obs = 1000      
K = 2       
J = 2        
dimension = 2  
vit = 1.2 
pdt = 0.1     

A = matrix(c(0.95,0.05,0.1,0.9),
           ncol = K,
           nrow = K,
           byrow = TRUE)


N = 75
liste_theta = list(Nu(matrix(c(5, -5, 1, -.9)/vit**2, ncol = K, nrow = J), vit), 
                   Nu(matrix(c(5, -5, -5, 5)/vit**2, ncol = K, nrow = J), vit),
                   Nu(matrix(c(5, -5, -1.2, 0.6)/vit**2, ncol = K, nrow = J), vit))
                   

# liste_theta = list(Nu(matrix(c(5,-5,2,-0.5),ncol =K, nrow = J),vit))
Resultats = data.frame(col_inutile = 1:N)
Resultats_A = data.frame(col_inutile = 1:N)
Resultats_Nu = data.frame(col_inutile = 1:N)
nbr_erreurs = 0

for (i in 1:length(liste_theta)){
  theta = liste_theta[[i]]
  Viterbi_EM = numeric(N)
  Viterbi_Flexmix = numeric(N)
  liste_norm_A_EM = numeric(N)
  liste_norm_A_Flexmix = numeric(N)
  liste_norm_nu_EM = numeric(N)
  liste_norm_nu_Flexmix = numeric(N)
  print(theta)
  
  compteur = 0
  while (compteur < N) {
    print(paste0('Tour ',compteur+1, ' pour Nu',i))
    Obs = Generation(nbr_obs, pdt, A, liste_cov, theta)$Obs
    
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
    
    m1 <- flexmix(deplacement ~ -1 + cov.1 + cov.2, k= 2, data= increments_dta)
    table(m1@cluster, increments_dta$etats_caches)
    parameters(m1)
    
    Init = initialisation2.0(increments = increments_dta, K = 2)
    A_init = Init$A[,c(2,3)]; param_init = Init$param; PI_init = Init$PI
    
    Relab = relabel(A_init, param_init, PI_init)
    
    Lambda = list('A' = Relab$A,
                  'param' = Relab$param,
                  'PI' = Relab$PI)
    
    E = EM_Langevin(increments = increments_dta, 
                    Lambda = Lambda, 
                    G = 20, 
                    moyenne = FALSE)
    
    if (!E$erreur){
      print("Pas d erreur")
      B = proba_emission(increments = increments_dta, param = E$param)
      V = Viterbi(E$A, B, E$PI)
      V_flexmix = Viterbi(Relab$A, proba_emission(increments_dta, Relab$param), Relab$PI)
      
      Viterbi_EM[compteur+1] = sum(Obs$etats_caches[1:(nbr_obs-1)]==V)/(nbr_obs-1)
      Viterbi_Flexmix[compteur+1] = sum(Obs$etats_caches[1:(nbr_obs-1)]==V_flexmix)/(nbr_obs-1)
      
      liste_norm_A_EM[compteur] = norm(A - E$A)
      liste_norm_A_Flexmix[compteur] = norm(A - Relab$A)
      liste_norm_nu_EM[compteur] = norm(theta - AffParams(E$param)$nu)
      liste_norm_nu_Flexmix[compteur] = norm(theta - AffParams(Relab$param)$nu)
      compteur = compteur + 1
      }else{nbr_erreurs = nbr_erreurs + 1}
  }
  # Test de Viterbi.
  Resultats[paste0('Nu',i,' EM')] = Viterbi_EM  
  Resultats[paste0('Nu',i,' Flexmix')] = Viterbi_Flexmix  
  
  # Test des normes.
  # A.
  Resultats_A[paste0('A EM Nu',i)] = liste_norm_A_EM
  Resultats_A[paste0('A Flexmix Nu',i)] = liste_norm_A_Flexmix
  
  # Nu.
  Resultats_Nu[paste0('Nu EM Nu',i)] = liste_norm_nu_EM
  Resultats_Nu[paste0('Nu Flexmix Nu',i)] = liste_norm_nu_Flexmix
  
}
# Je ne prends pas les résultats dans l'ordre car les thetas ne le sont pas.
boxplot(Resultats[,c(4,5,6,7,2,3)]*100, 
        col = c('red','orange','red','orange','red','orange'),
        main = paste0("Test de Viterbi (",nbr_obs," observations, ",N,
                      " répétitions)\n vitesses = (",vit,', ',vit,")"),
        at=c(1,2,4,5,7,8),
        names = c("Em Nu1",'Fl Nu1',"Em Nu2",'Fl Nu2',"Em Nu3",'Fl Nu3'),
        ylab = '% de réussite',
        outline = FALSE)
        
boxplot(Resultats_A[,c(4,5,6,7,2,3)],
        main = paste0("Test sur la prédiction de la matrice de transition (",nbr_obs,
                      " obs, ",N," gen)\n vitesses = (",vit,', ',vit,")"),
        col = c('red','orange','red','orange','red','orange'),
        at=c(1,2,4,5,7,8),
        names = c("Em Nu1",'Fl Nu1',"Em Nu2",'Fl Nu2',"Em Nu3",'Fl Nu3'),
        ylab = 'Erreurs',
        outline = FALSE)

boxplot(Resultats_Nu[,c(4,5,6,7,2,3)],
        main = paste0("Test sur la prédiction des Nu (",nbr_obs," obs, ",N
                      ," gen) \n vitesses = (",vit,', ',vit,")"),
        col = c('red','orange','red','orange','red','orange'),
        at=c(1,2,4,5,7,8),
        names = c("Em Nu1",'Fl Nu1',"Em Nu2",'Fl Nu2',"Em Nu3",'Fl Nu3'),
        ylab = 'Erreurs', 
        outline = FALSE)





################################################################################
wResultats = data.frame(col_inutile = 1:N)

N = 100 
compteur = 0 

boxplot(Resultats[-1],
        main = "Test des normes \n sur la prédiction de Flexmix et de l'EM \n 100 générations, 20 tours d'EM, Nu 1, pdt = 0.1" ,
        col = c('orange','red','orange','red'),
        ylab = 'erreur')









# install.packages("ggalluvial")
# library(ggalluvial)




# df = data.frame(Vit_EM = V, Vit_Flexmix = V_flexmix, Reel = Obs$etats_caches[1:(nbr_obs-1)] )


# ggplot(data = df,
#        aes(x = vit_em, stratum = stratum, alluvium = alluvium,
#            y = Freq, label = stratum)) +
#   geom_alluvium(aes(fill = Survived)) +
#   geom_stratum() + geom_text(stat = "stratum") +
#   theme_minimal() +
#   ggtitle("passengers on the maiden voyage of the Titanic",
#           "stratified by demographics and survival")
# 




Aff = AffParams(param_init)
Aff
nu_flexmix = Aff$nu
nu_EM = E$Nu


Mat_EM = theta - (nu_EM)
Mat_flexmix = theta - (nu_flexmix)
Mat_EM_r = theta - ret(nu_EM)
Mat_flexmix_r = theta - ret(nu_flexmix)

norme(Mat_EM)
norme(Mat_flexmix)
norme(Mat_EM_r)
norme(Mat_flexmix_r)

#Pour x :
## Pour k=1 :
k = 2
flexmix_y = m1@posterior[1]$scaled[500:998,k]
m1@cluster[1:499]
Gam = E$gam
plot(flexmix_x, Gam)
