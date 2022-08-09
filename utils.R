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
# 
# increments = function(liste){
#   # Ca n'est pas très efficace de faire grossir un vecteur au fur et à mesure, 
#   # ca fait bcp d'allocation memoire inutile et de copie dan sla zone mémoire
#   # n <-  length(liste) - 1
#   # l <-  numeric(length = n)
#   # for (i in 1:(length(liste)-1)){
#   #   l[i] <- liste[i+1] - liste[i]
#   # }
#   # Une autre option encore plus efficace
#   # l <- diff(liste)  
#   
#   l = c()
#   for (i in 1:(length(liste)-1)){
#     l = c(l,liste[i+1] - liste[i])
#   }
#   return(l)
# }

matrix_to_list = function(mat){
  # meme remarque que dans la fonction précédente, autant que possible quand tu connais la taille
  #   d'un objet il faut l'allouer des le début
  l = list()
  N = ncol(mat)
  for (i in 1:N){
    l[[i]] = c(mat[,i])
  }
  return(l)
}

list_to_matrix = function(liste){
  K = length(liste)
  J = length(liste[[1]])
  m = c()
  for (i in 1:K){
    m = c(m, liste[[i]])
  }
  return(matrix(m,J,K))
}
A = matrix(c(1,2,3,4,5,6),3,2)
list_to_matrix(matrix_to_list(A))


retourner = function(liste){
  l = length(liste)
  liste_retournee = c()
  for (i in l:1){
    liste_retournee = c(liste_retournee, liste[i])
  }
  return(liste_retournee)
}



################################################################################
###                   Algorithme EM
###
### Paramètres : 
### increments : dataframe contenant :
###   
### Lambda : liste contenant le modèle initial. 
### Vitesses : liste des vitesses initiales selon les etats. 
### G : nombre d iterations de l'EM. (Par défaut, 20).
### moyenne :
### 






# THETA.
theta_nv = matrix(1,J,K)
Vits = numeric(K)


Params = lapply(1:K, function(k){
  model = lm(increments$deplacement ~ -1 + C, weights = c(gam[,k],gam[,k]))
  return(list(nu = coef(model), vitesse = summary(model)$sigma))
})


for (k in 1:K){
  # On gere les deux cas differents selon la dimension.
  model = lm(increments_dta$deplacement ~ as.matrix(C), weights= c(gam[,k],gam[,k]))
  
  # On recupere les coefficients.
  Vits[k] = summary(model)$sigma
  theta_nv[,k] = coef(model)[2:(J+1)]
}
print(theta_nv)
# On range les nouveaux parametres dans une liste.
Params = list(Nu = theta_nv, vitesse = Vits)


















