
get_df_varname <- function(allnames){
  dfname <- data.frame(fullname = character(),
                       name = character(),
                       year = numeric(),
                       dim1 = numeric(),
                       dim2 = numeric(),
                       dim3 = numeric())
  
  
  # dim0 var ----------------------------------------------------------------
  
  listvar_dim0 <- allnames %>% 
    str_subset(fixed("["), negate = TRUE)
  if(length(listvar_dim0) > 0) {
    listdim0 <- get_dim0_var(listvar_dim0)
    dfname <- rbind(
      dfname,
      data.frame(fullname = listdim0$name,
                 name = listdim0$name,
                 dim1 = NA,
                 dim2 = NA,
                 dim3 = NA)
    )
  }
  
  
  # dim1 var ---------------------
  # var[year]
  
  listvar_dim1 <- allnames %>% 
    str_subset(fixed("[")) %>% 
    str_subset(fixed(","), negate = TRUE)
  
  if(length(listvar_dim1) > 0) {
    listdim1 <- get_dim1_var(listvar_dim1)
    
    dfname <- rbind(
      dfname,
      data.frame(fullname = listvar_dim1,
                 name = listdim1$name,
                 dim1 = listdim1$dim1,
                 dim2 = NA,
                 dim3 = NA)
    )
  }
  
  # dim 2 var --------------------------
  
  listvar_dim2 <-  allnames[str_count(allnames, pattern = fixed(",")) == 1]
  
  if(length(listvar_dim2) > 0) {
    
    listdim2 <- get_dim2_var(listvar_dim2)
    dfname <- rbind(
      dfname,
      data.frame(fullname = listdim2$fullname,
                 name = listdim2$name,
                 dim1 = listdim2$dim1,
                 dim2 = as.numeric(listdim2$dim2),
                 dim3 = NA)
    )
  }
  return(dfname)
  
  
}


# dim 0 var ----------------------------
get_dim0_var <- function(allnames){
  #  listvar_dim0 <-  str_subset(allnames,fixed("["),negate = TRUE)
  return(list("name" = allnames))
}

# dim 1 var ----------------------------
get_dim1_var <- function(allnames){
  # listvar_dim1 <- allnames %>% 
  #   str_subset(fixed("[")) %>% 
  #   str_subset(fixed(","), negate = TRUE)
  
  splitvar <- str_split(allnames, pattern = fixed("["))
  lapply(splitvar, function(mystr){
    mystr[1]
  }) %>%
    do.call("c",.) -> listname_dim1
  
  
  splitvar <- str_split(allnames, pattern = fixed("["))
  lapply(splitvar, function(mystr){
    mystr[2]
  }) %>% 
    do.call("c",.) %>% 
    str_remove(., fixed("]")) %>% 
    as.numeric() -> dim1
  
  return(list("fullname" = allnames,
              "name" = listname_dim1,
              "dim1" = dim1))
}

# dim 2 var ----------------------------
get_dim2_var <- function(allnames){
  splitvar <- str_split(allnames, pattern = fixed("["))
  lapply(splitvar, function(mystr){
    mystr[1]
  }) %>% do.call("c",.) -> listname_dim2
  
  lapply(splitvar, function(mystr){
    mystr[2]
  }) %>% 
    do.call("c",.) %>% 
    str_split(., pattern = fixed(",")) %>% 
    lapply(., function(mystr){
      mystr[1]
    }) %>%
    do.call("c",.) %>% 
    as.numeric() -> dim1
  
  lapply(splitvar, function(mystr){
    mystr[2]
  }) %>% 
    do.call("c",.) %>% 
    str_split(., pattern = fixed(",")) %>% 
    lapply(., function(mystr){
      mystr[2]
    }) %>%
    do.call("c",.) %>% 
    str_remove(., fixed("]")) %>% 
    str_remove(., fixed(" ")) -> dim2
  
  return(list("fullname" = allnames,
              "name" = listname_dim2,
              "dim1" = dim1,
              "dim2" = dim2))
}


# dim 3 var ----------------------------
get_dim3_var <- function(allnames){
  splitvar <- str_split(allnames, pattern = fixed("["))
  lapply(splitvar, function(mystr){
    mystr[1]
  }) %>% do.call("c",.) -> listname_dim3
  
  lapply(splitvar, function(mystr){
    mystr[2]
  }) %>% 
    do.call("c",.) %>% 
    str_split(., pattern = fixed(",")) %>% 
    lapply(., function(mystr){
      mystr[1]
    }) %>%
    do.call("c",.) %>% 
    as.numeric() -> dim1
  
  lapply(splitvar, function(mystr){
    mystr[2]
  }) %>% 
    do.call("c",.) %>% 
    str_split(., pattern = fixed(",")) %>% 
    lapply(., function(mystr){
      mystr[2]
    }) %>%
    do.call("c",.) %>% 
    as.numeric()-> dim2
  
  lapply(splitvar, function(mystr){
    mystr[2]
  }) %>% 
    do.call("c",.) %>% 
    str_split(., pattern = fixed(",")) %>% 
    lapply(., function(mystr){
      mystr[3]
    }) %>%
    do.call("c",.) %>% 
    str_remove(., fixed("]")) %>% 
    str_remove(., fixed(" ")) -> dim3
  
  return(list("fullname" = allnames,
              "name" = listname_dim3,
              "dim1" = dim1,
              "dim2" = dim2,
              "dim3" = dim3))
}

