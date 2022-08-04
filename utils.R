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
  # Ca n'est pas très efficace de faire grossir un vecteur au fur et à mesure, 
  # ca fait bcp d'allocation memoire inutile et de copie dan sla zone mémoire
  # n <-  length(liste) - 1
  # l <-  numeric(length = n)
  # for (i in 1:(length(liste)-1)){
  #   l[i] <- liste[i+1] - liste[i]
  # }
  # Une autre option encore plus efficace
  # l <- diff(liste)  
  
  l = c()
  for (i in 1:(length(liste)-1)){
    l = c(l,liste[i+1] - liste[i])
  }
  return(l)
}

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

retourner = function(liste){
  l = length(liste)
  liste_retournee = c()
  for (i in l:1){
    liste_retournee = c(liste_retournee, liste[i])
  }
  return(liste_retournee)
}

