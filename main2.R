RFoptions(install="no")
source('utils.R')
source('simu.R')
source('estim.R')
load("lise_cov_save.RData")

################################################################################
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














################################################################################

nbr_obs = 1000      
K = 2       
J = 2     
pdt = 0.1

lim <- c(-30, 30, -30, 30) # limits of map
resol <- 0.1 # grid resolution
rho <- 2; nu <- 1.5; sigma2 <- 30# Matern covariance parameters
mean_function <- function(z){# mean function
  -log(3 + sum(z^2))}


A = matrix(c(0.95,0.05,0.1,0.9),
           ncol = K,
           nrow = K,
           byrow = TRUE)

v1 = 5
v2 = -5
beta_sim <- matrix(c(v1, -v1, v2, -v2), ncol = K, nrow = J)
theta = Nu(beta_sim, vit)
theta

N = 50
c = 0
CPTR = 0

ano = round(nbr_obs/100*5)   # nombre d'anomalies.
tps = temps(pdt, nbr_obs, ano)


liste_resultats = as.list(numeric(N))

while (c < N){
  print(c)
  # GEN = Generation_observation3.0(matrix_to_list(theta), 
  #                                 etats_caches = CM_generateur( A, nbr_obs),
  #                                 liste_cov = liste_cov,
  #                                 Vits = c(0.4,0.4),
  #                                 tps = tps   
  #                                 )
  GEN = Generation(nbr_obs, pdt, A, liste_cov, theta)
  Obs = GEN$Obs
  increments_dta = GEN$increments_dta
  
  m1 <- flexmix(deplacement ~ -1 + cov.1 + cov.2, k= 2, data= increments_dta)
  table(m1@cluster, increments_dta$etats_caches)
  parameters(m1)
  
  Init = initialisation2.0(increments = increments_dta, K = 2)
  A_init = Init$A[,c(2,3)]; param_init = Init$param; PI_init = Init$PI
  
  Relab = relabel(A_init, param_init, PI_init)
  
  Lambda = list(A = Relab$A,
                param = Relab$param,
                PI = Relab$PI)
  
  E = EM_Langevin(increments = increments_dta, 
                  Lambda = Lambda, 
                  G = 7, 
                  moyenne = FALSE)
  
  
  Aff = AffParams(param_init)
  
  N_EM = norme(theta - (E$Nu))
  N_flexmix = norme(theta - (Aff$nu))
  N_flexmix_r = norme(theta - ret(Aff$nu))
  
  if (N_EM < N_flexmix && N_EM < N_flexmix_r){
    CPTR = CPTR + 1
    Aff$test = T
  }
  liste_resultats[[c+1]] = Aff
  c = c + 1
}
CPTR/N


