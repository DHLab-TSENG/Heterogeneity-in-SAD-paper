cramerV_matrix <- 
  function(dataset = dataset, vars = vars) {
    cramerV_result <- 
      matrix(data     = NA, 
             nrow     = length(vars),
             ncol     = length(vars),
             dimnames = list(vars,vars))
    for (i in seq_along(vars)) {
      for (j in seq_along(vars)) {
        
        cramerV_result[i,j] <- 
          dataset[,table(.SD), .SDcols = c(vars[i],vars[j])] %>% 
          assocstats(.) %>% 
          .["cramer"] %>% 
          unname(.) %>% 
          unlist(.) %>% 
          round(.,5)
      }
    }
    data.frame(cramerV_result)
  }