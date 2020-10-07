## Functions for estimating c.hat based on bootstrapping
## #############################################################################

## Source:
##  https://sites.google.com/site/cmrsoftware/lecture-lab-schedule/6---model-fit-and-mmi/methods-for-assessing-fit/dipper.gof.R?attredirects=0
## Modified February 2020, Meinolf Ottensmann

## function to get number of new releases for each group * occasion
## data = data frame
## n.ocassion = number of capture-resight events
## groups = group levels
## group.col = column name of group
Marked.burnham <- function(data = NULL,
                           n.occasions = NULL,
                           groups = c("Female", "Male"), group.col = "sex") {
  
  #convert ldld to ll ch history as in CJS model
  data$ch <- sapply(1:nrow(data), function(x) {
    temp <- stringr::str_split(data$ch[x], "")[[1]]
    paste0(temp[seq(1, length(temp),2)], collapse = "")
  })
  
  group <- data[, group.col]
  marked <- matrix(nrow = length(groups), ncol = n.occasions)
  for (g in 1:length(groups)) {
    data_ <- subset(data, group == groups[g])
    data_
    ch <- data_$ch
    for (i in 1:n.occasions) {
      ch1 <- ch[(as.numeric(substr(ch,1,i))) == 1]
      marked[g,i] <- length(ch1)
    }
  }
  return(marked)
}


## simulate processes S, p, F and r
simul.burnham <- function(S, p, F, r, marked) {
  n.occasions <- length(p) + 1 
  S <- matrix(S, n.occasions - 1, nrow = sum(marked), byrow = T)
  P <- matrix(p, n.occasions - 1, nrow = sum(marked), byrow = T)
  F <- matrix(F, n.occasions - 1, nrow = sum(marked), byrow = T)
  r <- matrix(r, n.occasions - 1, nrow = sum(marked), byrow = T)
  
  CH <- matrix(0, ncol = 2 * n.occasions, nrow = sum(marked))
  # define a vector with marking occasion 
  mark.occ <- rep(seq(1, length(marked), 1), marked[1:length(marked)])
  
  # Corresponding column in Capture history follows
  alive.slot = function(i) i +  (-1 + i)
  dead.slot = function(i) i +  (-2 + i)
  
  # fill in CH
  for (i in 1:sum(marked)) { # i = individual
    CH[i, alive.slot(mark.occ[i])] <- 1 # mark the capture 
    if (mark.occ[i] == n.occasions) next # if last occasion reached
    for (t in (mark.occ[i] + 1):n.occasions) { ## all possible occasions + 1
      #survive?
      sur <- rbinom(1, 1 , S[i, t - 1])
      if (sur == 0) {
        # noticed or not
        death.rec <- rbinom(1, 1, r)
        if (death.rec == 1) CH[i, dead.slot(t)] <- 1
        if (sur == 0) break # check if expression is correct!
      } else if (sur == 1) {
        ## stays or leaves permanently ?
        stays <- rbinom(1,1, F[i, t - 1])  
        if (stays == 0) break # check if expression is correct!
        ## resighted ?
        #recaptured?
        rp <- rbinom(1,1,P[i,t - 1])
        if (rp == 1) CH[i, alive.slot(t)] <- 1
      }
    } #t
  } #i
  return(CH)
}

#function to create capture history character strings
pasty <- function(x) {
  k <- ncol(x)
  n <- nrow(x)
  out <- array(dim = n)
  for (i in 1:n) {
    out[i] <- paste(x[i,  ],  collapse = "")
  }
  return(out)
}


#estimates go into simulation within loop
sims <- function(Phi.fem, Phi.mal, p.fem, p.mal, marked,reps, group.col = "sex", model.parameters = NULL) {
  deviance <- dim(reps)
  for (i in 1:reps) {
    cat("iteration = ", i, "\n")
    sim.fem <- simul.cjs(Phi.fem,p.fem,marked[1,])
    sim.mal <- simul.cjs(Phi.mal,p.mal,marked[2,])
    fem.hist <- data.frame(ch = pasty(sim.fem), sex = "Female")
    mal.hist <- data.frame(ch = pasty(sim.mal), sex = "Male")
    sim.data <- rbind(fem.hist, mal.hist)
    sim.processed = process.data(sim.data, model = "CJS", groups = "sex")
    sim.ddl = make.design.data(sim.processed)
    global.sim <- mark(sim.processed,sim.ddl, 
                       model.parameters = model.parameters, output = F, silent = T)
    
    
    deviance[i] <- global.sim$results$deviance 
  }
  out <- list(deviance.mean = mean(deviance),
              deviance.025 = quantile(deviance, 0.025), 
              deviance.975 = quantile(deviance, 0.975))
}


sims.burnham <- function(S.fem, S.mal, p.fem, p.mal, F.fem, F.mal, r.mal, r.fem,
                         marked, reps, group.col = "sex", model.parameters = NULL) {
  deviance <- dim(reps)
  for (i in 1:reps) {
    cat("iteration = ", i, "\n")
    sim.fem <- simul.burnham(S = S.fem, p = p.fem, F = F.fem, r = r.fem, marked[1,])
    sim.mal <- simul.burnham(S = S.mal, p = p.mal, F = F.mal, r = r.mal, marked[2,])
    
    fem.hist <- data.frame(ch = pasty(sim.fem), sex = "Female")
    mal.hist <- data.frame(ch = pasty(sim.mal), sex = "Male")
    sim.data <- rbind(fem.hist, mal.hist)
    sim.processed = process.data(sim.data, model = "Burnham", groups = "sex")
    sim.ddl = make.design.data(sim.processed)
    global.sim <- mark(sim.processed,sim.ddl, 
                       model.parameters = model.parameters, output = F, silent = T)
    
    
    deviance[i] <- global.sim$results$deviance 
  }
  out <- list(deviance.mean = mean(deviance),
              deviance.025 = quantile(deviance, 0.025), 
              deviance.975 = quantile(deviance, 0.975))
}



# Wraper to Rmar::covariate.predictions with further data wrangling
# models: List of CJS models
# data: data frame with covariates to predict
# index.meta data frame with model.index information
covar.pred <- function(models = NULL, df = NULL, index.meta = index.meta, ddl = ddl) {
  
  ## obtain predictions from model(s)
  mod.pred <- covariate.predictions(  
    model = models, 
    data = df,
    #drop = T,
    indices = index.meta$model.index)
  
  ## fix model.index and subset 
  ests <- mod.pred$estimates %>% 
    mutate(model.index = par.index) %>%  # renaming see above
    subset(., select = -c(vcv.index, par.index)) %>% 
    unique.data.frame()
  
  ## add meta information
  ests <- 
    left_join(ests, ddl$Phi[,c("model.index", "age")], by = c("model.index")) %>% 
    mutate(age = age.y) %>% 
    select(., subset = -c(age.x, age.y, fixed, se)) %>% 
    mutate(dummy = paste0(model.index, Sex, age)) %>% 
    filter(., dummy %in% index.meta$dummy) %>% 
    select(., subset = -c(dummy)) %>%
    unique.data.frame()
  
  return(ests)
  
}

# Create cohort function
create.cohort = function(x) {
  # split the capture histories into a list with each list element
  # being a vector of the occasion values (0/1). 1001 becomes
  # "1","0","0","1"
  split.ch = sapply(x$ch, strsplit, split = "")
  # combine these all into a matrix representation for the ch
  # rows are animals and columns are occasions
  chmat = do.call("rbind", split.ch)
  # use the defined function on the rows (apply(chmat,1...) of the
  # matrix. The defined function figures out the column containing
  # the first 1 (its initial release column).
  return(factor(apply(chmat, 1, function(x) min(which(x != "0"))))) }
