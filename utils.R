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



##create covariable
create_covariate_columns <- function(df){
  ## extarct covariates names
  cov_names <- str_extract(colnames(df), "grad_c[0123456789]+?_x") %>% 
    na.omit() %>% as_tibble() %>% 
    mutate(cov = str_remove(value, "grad_c")) %>% 
    mutate(cov = str_remove(cov, "_x")) %>% 
    mutate(cov = as.numeric(cov)) %>% 
    pull(cov)
  
  ## function to select grad_dx or gradient_dy 
  test <- function(name, df ){
    nomx <- paste0("grad_c",name,"_x")
    nomy <- paste0("grad_c",name,"_y")
    covname <- paste0("cov",name)
    mutate(df, !!covname:= ifelse(dimension=="dep_x", !!sym(nomx), !!sym(nomy))) %>% 
      pull(!!covname)
  }
  
  cov_to_add <- cov_names %>% 
    map_dfc(.f = test, df= df) %>% 
    rename_with(~str_c("cov", str_remove(., ".."), sep = ""), matches( "..[0123456789]+?")) 
  
  
  df %>% bind_cols(cov_to_add) %>% 
    mutate(across(matches("cov.[1234567890]+?"),.fns =  ~ 0.5*.x * sqrt(delta_t))) 
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

ret = function(nu){return(cbind(nu[,2], nu[,1]))}
norme = function(mat){return(tr(mat %*% aperm(mat)))}
AffParams = function(param){
  K = length(param)
  J = length(param[[1]]$nu)
  nu = c()
  vitesses = numeric(K)
  for (k in 1:K){
    nu = c(nu, as.numeric(param[[k]]$nu))
    vitesses[k] = param[[k]]$vitesse
  }
  return(list(nu = matrix(nu, nrow = J, ncol = K),
              vitesses = vitesses))
}

