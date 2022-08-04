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

# Transforme une matrice en une liste des colonnes. 
matrix_to_list = function(mat){
  N = ncol(mat)
  return( lapply(1:N, function(d){mat[,d]}))
}

# Fonction qui retourner une liste. Elle peut être utile après l'utilisation 
# de l'algorithme de Viterbi.
reverse = function(liste){
  l = length(liste)
  liste_retournee = numeric(length = l)
  for (i in l:1){
    liste_retournee[l-i] = liste[i]
  }
  return(liste_retournee)
}

