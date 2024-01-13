## Custom functions used for analyses of common buzzard survival:
## Functions are sorted in alphabetical order.
## 
## Implementation of median c-hat simulation is derived from:
## https://sites.google.com/site/cmrsoftware/lecture-lab-schedule/6---model-fit-and-mmi/methods-for-assessing-fit/dipper.gof.R?attredirects=0 
## 
## -----------------------------------------------------------------------------

#' @description Function that creates live-dead encounter histories
#' @param start first year of sampling period
#' @param last last year of sampling period
#' @param ringing_data data frame with ringing data of individuals
#' @param resights data frame with re-sighting data of individuals
create.ld <- function(start, end, ringing_data, resights) {
  ## Identify all individuals 
  individuals <- ringing_data$Ring
  
  ## define Study periods based on start and end years
  years <- c(start, end)
  periods <- rep(seq(start, end), each = 1)
  
  ## create matrix 
  ## rows: Individuals
  ## cols: 2 Periods per year with Live and Dead encounters respectively
  mat <- matrix(0, nrow = length(individuals), ncol = length(periods)*2) %>% 
    set_colnames(paste0(rep(periods, each = 2), c("L", "D"))) %>% 
    set_rownames(individuals)
  
  ## code year of capture as 1,
  ## Code sightings in subsequent intervals by 1 too
  for (x in individuals) {
    year.ringed <- filter(ringing_data, Ring == x)[["Year"]][1]
    ## code birth year in living encounters
    mat[rownames(mat) == x, colnames(mat) == paste0(year.ringed, "L")] <- 1  
    
    ## Sampling intervals follow the breeding cycle (May - April).
    ## Hence, spring month count for the previous year, i.e. March 2021 -> 2020
    year.resighted <- 
      filter(resights, Ring == x)[,c("Year", "Month")] %>% 
      mutate(period = case_when(
        Month == "01" ~ paste0(as.numeric(Year) - 1 , "L"),
        Month == "02" ~ paste0(as.numeric(Year) - 1 , "L"),
        Month == "03" ~ paste0(as.numeric(Year) - 1 , "L"),
        Month == "04" ~ paste0(as.numeric(Year) - 1 , "L"),
        Month == "05" ~ paste0(Year, "L"),
        Month == "06" ~ paste0(Year, "L"),
        Month == "07" ~ paste0(Year, "L"),
        Month == "08" ~ paste0(Year, "L"),
        Month == "09" ~ paste0(Year, "L"),
        Month == "10" ~ paste0(Year, "L"),
        Month == "11" ~ paste0(Year, "L"),
        Month == "12" ~ paste0(Year, "L")
      ))
    
    ## write re-sights to matrix 
    for (col in year.resighted$period) {
      mat[rownames(mat) == x, colnames(mat) == col] <- 1
    }
    
    ## check if reported as dead
    death <- filter(ringing_data, Ring == x)[,c("Dead", "DateDeath")][1,]
    if (!is.na(death$Dead)) {
      death <- death %>% 
        mutate(Month = lubridate::month(death$DateDeath)) %>% 
        mutate(Year = lubridate::year(death$DateDeath)) %>% 
        mutate(period = case_when(
          Month == "1" ~ paste0(as.numeric(Year) - 1 , "D"),
          Month == "2" ~ paste0(as.numeric(Year) - 1 , "D"),
          Month == "3" ~ paste0(as.numeric(Year) - 1 , "D"),
          Month == "4" ~ paste0(as.numeric(Year) - 1 , "D"),
          Month == "5" ~ paste0(Year, "D"),
          Month == "6" ~ paste0(Year, "D"),
          Month == "7" ~ paste0(Year, "D"),
          Month == "8" ~ paste0(Year, "D"),
          Month == "9" ~ paste0(Year, "D"),
          Month == "10" ~ paste0(Year, "D"),
          Month == "11" ~ paste0(Year, "D"),
          Month == "12" ~ paste0(Year, "D")
        ))
      col <- death[["period"]]
      mat[rownames(mat) == x, colnames(mat) == col] <- 1
    }
    
  }
  
  ## coerce to character vector
  return(apply(mat[,], 1, function(x) paste0(x, collapse = "")) %>% 
           as.character())
  
}


#' @description Translate LD capture histories to release corhorts
#' @param data data frame with columns ch and sex
#' @param n.occassion number of occasions described by ch  
#' @param groups vector of group levels
#' @param group.col grouping variable
Marked <- function(data, n.occasions, groups, group.col) {
  
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



#' @description  Condense vector to string
#' @param x vector
pasty <- function(x) {
  k <- ncol(x)
  n <- nrow(x)
  out <- array(dim = n)
  for (i in 1:n) {
    out[i] <- paste(x[i,  ],  collapse = "")
  }
  return(out)
}

#' @description Simulate encounter histories (ch) based on original marking events
#' @details Uses estimates from the empirical model to simulate LD encouter sequences
#' @param S.fem data frame with S parameters for females
#' @param S.mal data frame with S parameters for males
#' @param p.fem data frame with p parameters for females
#' @param p.mal data frame with S parameters for males
#' @param r.fem data frame with r parameters for females
#' @param r.mal data frame with S parameters for males
#' @param F.fem F parameter value
#' @param F.mal F parameter value
#' @param marked matrix of releases (row 1 = females, row 2 = males)
#' @param model.parameters list of model parameters
#' 
sims.burnham <- function(S.fem, S.mal, p.fem, p.mal, F.fem, F.mal, r.mal, r.fem,
                         marked, reps, model.parameters) {
  
  ## set up vector to save deviance estimates
  deviance <- dim(reps)
  for (i in 1:reps) {
    cat("iteration ", i, "out of", reps, "\n")
    
    ## Simulate capture histories, separately for groups
    ## ----------------------------------------------------------------
    sim.femDInf <- simul.burnham(S = S.fem, p = p.fem, F = F.fem, r = r.fem, marked[1,])
    # 1
    fem.histDInf <- data.frame(ch = pasty(sim.femDInf), Sex = "Female", Morph = "D", Lbinom = "Infected")
    sim.femDUnInf <- simul.burnham(S = S.fem, p = p.fem, F = F.fem, r = r.fem, marked[7,])
    #2
    fem.histDUnInf <- data.frame(ch = pasty(sim.femDUnInf), Sex = "Female", Morph = "D", Lbinom = "Uninfected")
    
    sim.femIInf <- simul.burnham(S = S.fem, p = p.fem, F = F.fem, r = r.fem, marked[3,])
    #3
    fem.histIInf <- data.frame(ch = pasty(sim.femIInf), Sex = "Female", Morph = "I", Lbinom = "Infected")
    sim.femIUnInf <- simul.burnham(S = S.fem, p = p.fem, F = F.fem, r = r.fem, marked[9,])
    #4
    fem.histIUnInf <- data.frame(ch = pasty(sim.femIUnInf), Sex = "Female", Morph = "I", Lbinom = "Uninfected")
    
    sim.femLInf <- simul.burnham(S = S.fem, p = p.fem, F = F.fem, r = r.fem, marked[5,])
    #5
    fem.histLInf <- data.frame(ch = pasty(sim.femLInf), Sex = "Female", Morph = "L", Lbinom = "Infected")
    sim.femLUnInf <- simul.burnham(S = S.fem, p = p.fem, F = F.fem, r = r.fem, marked[11,])
    #6
    fem.histLUnInf <- data.frame(ch = pasty(sim.femLUnInf), Sex = "Female", Morph = "L", Lbinom = "Uninfected")
    
    sim.malDInf <- simul.burnham(S = S.mal, p = p.mal, F = F.mal, r = r.mal, marked[2,])
    #7
    mal.histDInf <- data.frame(ch = pasty(sim.malDInf), Sex = "Male", Morph = "D", Lbinom = "Infected")
    sim.malDUnInf <- simul.burnham(S = S.mal, p = p.mal, F = F.mal, r = r.mal, marked[8,])
    #8
    mal.histDUnInf <- data.frame(ch = pasty(sim.malDUnInf), Sex = "Male", Morph = "D", Lbinom = "Uninfected")
    
    sim.malIInf <- simul.burnham(S = S.mal, p = p.mal, F = F.mal, r = r.mal, marked[4,])
    #9
    mal.histIInf <- data.frame(ch = pasty(sim.malIInf), Sex = "Male", Morph = "I", Lbinom = "Infected")
    sim.malIUnInf <- simul.burnham(S = S.mal, p = p.mal, F = F.mal, r = r.mal, marked[10,])
    #10
    mal.histIUnInf <- data.frame(ch = pasty(sim.malIUnInf), Sex = "Male", Morph = "I", Lbinom = "Uninfected")
    
    sim.malLInf <- simul.burnham(S = S.mal, p = p.mal, F = F.mal, r = r.mal, marked[6,])
    #11
    mal.histLInf <- data.frame(ch = pasty(sim.malLInf), Sex = "Male", Morph = "L", Lbinom = "Infected")
    sim.malLUnInf <- simul.burnham(S = S.mal, p = p.mal, F = F.mal, r = r.mal, marked[12,])
    #12
    mal.histLUnInf <- data.frame(ch = pasty(sim.malLUnInf), Sex = "Male", Morph = "L", Lbinom = "Uninfected")
    
    
    ## Merge datasets and process for burnham model
    ## ----------------------------------------------------------------
    sim.data <- rbind(fem.histDInf, fem.histDUnInf, fem.histIInf, fem.histIUnInf,
                      fem.histLInf, fem.histLUnInf, mal.histDInf, mal.histDUnInf,
                      mal.histIInf, mal.histIUnInf, mal.histLInf, mal.histLUnInf) %>% 
      mutate(Lbinom = factor(Lbinom)) %>% 
      mutate(Sex = factor(Sex)) %>% 
      mutate(Morph = factor(Morph))
    sim.processed = process.data(sim.data, model = "Burnham",
                                 groups = c("Sex", "Morph", "Lbinom"),
                                 begin.time = 2007)
    
    ## Add design data. Coded the same way as for the empirical data    
    ## ----------------------------------------------------------------
    sim.ddl = make.design.data(sim.processed, remove.unused = FALSE)
    sim.ddl <- add.design.data(data = sim.processed,
                               ddl = sim.ddl,
                               parameter="S", 
                               type = "age",
                               bins = c(0, 1, 2, 3, 13),
                               right = FALSE,
                               name = "age",
                               replace = TRUE)
    sim.ddl <- add.design.data(data = sim.processed,
                               ddl = sim.ddl,
                               parameter="p", 
                               type = "age",
                               bins = c(0, 1, 2, 3, 13),
                               right = FALSE,
                               name = "age",
                               replace = TRUE)
    sim.ddl <- add.design.data(data = sim.processed,
                               ddl = sim.ddl,
                               parameter="r", 
                               type = "age",
                               bins = c(0, 1, 2, 3, 13),
                               right = FALSE,
                               name = "age",
                               replace = TRUE)
    
    ## Add second age-class argument for r parameter
    sim.ddl$r <- sim.ddl$r %>% 
      mutate(ageclass2 = case_when(
        age == "[0,1)" ~ "CY1",
        age != "[0,1)" ~ "CY2plus"
      ))
    ## Run the model on the simulated dataset
    ## ----------------------------------------------------------------
    global.sim <- mark(
      sim.processed, sim.ddl, model.parameters = model.parameters,
      output = FALSE, silent = TRUE, invisible = TRUE)
    
    ## extract and save model deviance
    ## ----------------------------------------------------------------
    deviance[i] <- global.sim$results$deviance 
  }
  ## summarise bootstrapping results
  out <- list(
    deviance.raw = deviance,
    deviance.median = median(deviance),
    deviance.025 = quantile(deviance, 0.025), 
    deviance.975 = quantile(deviance, 0.975))
  return(out)
}

#' @description Simulate encounter histories (ch) based on original marking events
#' @details Uses estimates from the empirical model to simulate LD encounter sequences
#' @param S.fem data frame with S parameters for females
#' @param S.mal data frame with S parameters for males
#' @param R.fem data frame with R parameters for females
#' @param R.mal data frame with R parameters for males
#' @param RPrime.fem data frame with R' parameters for females
#' @param RPrime.mal data frame with R' parameters for males
#' @param r.fem data frame with r parameters for females
#' @param r.mal data frame with S parameters for males
#' @param marked matrix of releases (row 1 = females, row 2 = males)
#' @param model.parameters list of model parameters
#' 
sims.barker <- function(S.fem, S.mal, R.fem, R.mal, RPrime.fem, RPrime.mal, r.mal, r.fem,
                        marked, reps, model.parameters) {
  
  ## set up vector to save deviance estimates
  deviance <- dim(reps)
  for (i in 1:reps) {
    cat("iteration ", i, "out of", reps, "\n")
    
    ## Simulate capture histories, separately for groups
    ## ----------------------------------------------------------------
    sim.femDInf <- simul.barker(S = S.fem, R = R.fem, RPrime = RPrime.fem, r = r.fem, marked = marked[1,])
    fem.histDInf <- data.frame(ch = pasty(sim.femDInf), Sex = "Female", Morph = "D", Lbinom = "Infected")
    sim.femDUnInf <- simul.barker(S = S.fem, R = R.fem, RPrime = RPrime.fem, r = r.fem, marked = marked[7,])
    fem.histDUnInf <- data.frame(ch = pasty(sim.femDUnInf), Sex = "Female", Morph = "D", Lbinom = "Uninfected")
    sim.femIInf <- simul.barker(S = S.fem, R = R.fem, RPrime = RPrime.fem, r = r.fem, marked = marked[3,])
    fem.histIInf <- data.frame(ch = pasty(sim.femIInf), Sex = "Female", Morph = "I", Lbinom = "Infected")
    sim.femIUnInf <- simul.barker(S = S.fem, R = R.fem, RPrime = RPrime.fem, r = r.fem, marked = marked[9,])
    fem.histIUnInf <- data.frame(ch = pasty(sim.femIUnInf), Sex = "Female", Morph = "I", Lbinom = "Uninfected")
    sim.femLInf <- simul.barker(S = S.fem, R = R.fem, RPrime = RPrime.fem, r = r.fem, marked = marked[5,])
    fem.histLInf <- data.frame(ch = pasty(sim.femLInf), Sex = "Female", Morph = "L", Lbinom = "Infected")
    sim.femLUnInf <- simul.barker(S = S.fem, R = R.fem, RPrime = RPrime.fem, r = r.fem, marked = marked[11,])
    fem.histLUnInf <- data.frame(ch = pasty(sim.femLUnInf), Sex = "Female", Morph = "L", Lbinom = "Uninfected")
    
    sim.malDInf <- simul.barker(S = S.mal, R = R.mal, RPrime = RPrime.mal, r = r.mal, marked = marked[2,])
    mal.histDInf <- data.frame(ch = pasty(sim.malDInf), Sex = "Male", Morph = "D", Lbinom = "Infected")
    sim.malDUnInf <- simul.barker(S = S.mal, R = R.mal, RPrime = RPrime.mal, r = r.mal, marked = marked[8,])
    mal.histDUnInf <- data.frame(ch = pasty(sim.malDUnInf), Sex = "Male", Morph = "D", Lbinom = "Uninfected")
    sim.malIInf <- simul.barker(S = S.mal, R = R.mal, RPrime = RPrime.mal, r = r.mal, marked = marked[4,])
    mal.histIInf <- data.frame(ch = pasty(sim.malIInf), Sex = "Male", Morph = "I", Lbinom = "Infected")
    sim.malIUnInf <- simul.barker(S = S.mal, R = R.mal, RPrime = RPrime.mal, r = r.mal, marked = marked[10,])
    mal.histIUnInf <- data.frame(ch = pasty(sim.malIUnInf), Sex = "Male", Morph = "I", Lbinom = "Uninfected")
    sim.malLInf <- simul.barker(S = S.mal, R = R.mal, RPrime = RPrime.mal, r = r.mal, marked = marked[6,])
    mal.histLInf <- data.frame(ch = pasty(sim.malLInf), Sex = "Male", Morph = "L", Lbinom = "Infected")
    sim.malLUnInf <- simul.barker(S = S.mal, R = R.mal, RPrime = RPrime.mal, r = r.mal, marked = marked[12,])
    mal.histLUnInf <- data.frame(ch = pasty(sim.malLUnInf), Sex = "Male", Morph = "L", Lbinom = "Uninfected")
    
    ## Merge datasets and process for burnham model
    ## ----------------------------------------------------------------
    sim.data <- rbind(fem.histDInf, fem.histDUnInf, fem.histIInf, fem.histIUnInf,
                      fem.histLInf, fem.histLUnInf, mal.histDInf, mal.histDUnInf,
                      mal.histIInf, mal.histIUnInf, mal.histLInf, mal.histLUnInf) %>% 
      mutate(Lbinom = factor(Lbinom)) %>% 
      mutate(Sex = factor(Sex)) %>% 
      mutate(Morph = factor(Morph))
   
    dp.barker  <- process.data(data = sim.data, model = "Barker", begin.time = 2007,
                               groups = c("Sex", "Morph", "Lbinom"))
    
    ## create design data
    ## ------------------------------------------------------
    ddl.barker <- make.design.data(dp.barker) 
    
    ## Bin age-classes
    ## Survival S:
    ## ------------------------------------------------------ 
    ddl.barker <- add.design.data(
      data = dp.barker, ddl = ddl.barker, parameter="S", type = "age",
      name = "age", bins = c(0, 1, 2, 13),
      right = FALSE, replace = TRUE)
    
    ## Further parameters
    ## -----------------------------------------------------------------------------
    ddl.barker <- add.design.data(
      data = dp.barker, ddl = ddl.barker, parameter="R", type = "age",
      name = "age", bins = c(0, 1, 2, 13),
      right = FALSE, replace = TRUE)
    
    ddl.barker <- add.design.data(
      data = dp.barker, ddl = ddl.barker, parameter="r", type = "age",
      name = "age", bins = c(0, 1, 2, 13),
      right = FALSE, replace = TRUE)
    
    ddl.barker <- add.design.data(
      data = dp.barker, ddl = ddl.barker, parameter="RPrime", type = "age",
      name = "age", bins = c(0, 1, 2, 13),
      right = FALSE, replace = TRUE)
    
    ddl.barker$r <- ddl.barker$r %>% 
      mutate(ageclass2 = case_when(
        age == "[0,1)" ~ "Juvenile",
        age != "[0,1)" ~ "Adult"
      ))
    
    ddl.barker$R <- ddl.barker$R %>% 
      mutate(ageclass2 = case_when(
        age == "[0,1)" ~ "Juvenile",
        age != "[0,1)" ~ "Adult"
      ))
    
    ddl.barker$RPrime <- ddl.barker$RPrime %>% 
      mutate(ageclass2 = case_when(
        age == "[0,1)" ~ "Juvenile",
        age != "[0,1)" ~ "Adult"
      ))
   
     ## Run the model on the simulated dataset
    ## ----------------------------------------------------------------
    global.sim <- mark(
      dp.barker, ddl.barker, model.parameters = model.parameters,
      output = FALSE, silent = TRUE, invisible = TRUE)
    
    ## extract and save model deviance
    ## ----------------------------------------------------------------
    deviance[i] <- global.sim$results$deviance 
  }
  ## summarise bootstrapping results
  out <- list(
    deviance.raw = deviance,
    deviance.median = median(deviance),
    deviance.025 = quantile(deviance, 0.025), 
    deviance.975 = quantile(deviance, 0.975))
  return(out)
}



#' Simulate capture histories using a set of parameters
#' @param S data frame 
#' @param p data frame
#' @param F numeric value
#' @param r data frame
#' @param marked vector of releases per cohort
#' 
## simulate processes S, p, F and r
simul.burnham <- function(S, p, F, r, marked) {
  
  ## 1.) Set-up parameters
  ## -----------------------------------------------------------------------------
  int <- min(S$time):max(S$time) 
  oc <- length(int) 
  
  ## 2.) Initialise capture history matrix
  ## -----------------------------------------------------------------------------
  CH <- matrix(data = 0, ncol = 2 * oc, nrow = sum(marked) * oc) %>% 
    set_colnames(paste0(c("L", "D"), rep(int, each = 2)))
  
  # define a vector with marking occasion 
  mark.occ <- rep(seq(1, length(marked), 1), marked[1:length(marked)])
  
  # Functions to link to ch locations
  alive.slot <- function(i) i +  (-1 + i)
  dead.slot <- function(i) i +  (-2 + i)
  
  ## which columns are live encounters
  live <- sapply(1:oc, alive.slot)
  ## which columns are dead encounter
  dead <- live + 1
  
  ## Code releases in CH
  CH <- lapply(1:length(mark.occ), function(i) {
    CH[i, alive.slot(mark.occ[i])] <- 1
    return(CH[i,])
  }) %>% 
    do.call("rbind",.)
  
  ## 3.) Run simulation
  ## -----------------------------------------------------------------------------
  CH <- lapply(1:nrow(CH), function(ind) {
    CHx <- CH[ind,]
    ## column of the release
    col <- which(CHx == 1)
    ## corresponding mark.occ
    release <- which(live == col)
    
    ## Simulate fate of individual
    for (t in release:oc) { 
      ## get covariates to pick parameters 
      ## focal year
      time <- int[t]
      
      # 1.) Survival: Function of time and age-class 
      # ------------------------------------------------------------------------
      # Check if individual is within the young age groups
      # change if more than 3 age groups are used!
      if ((t - release) < 2) {
        Sx <- filter(S, time == (int[t]), age == (t - release))[["estimate"]]  
      } else {
        Sx <- filter(S, time == (int[t]), age >= (t - release))[["estimate"]]  
      }
      if (length(Sx) != 1) stop("Check Sx")
      
      sur <- rbinom(1, 1 , Sx)
      if (sur == 0) {
        # 2.) If dead, recorded or not: Functions of two age classes
        # ----------------------------------------------------------------------
        if ((t - release) == 0) {
          rx <- filter(r, age == 0)[["estimate"]]
        } else if ((t - release) >= 1) {
          rx <- filter(r, age == 1)[["estimate"]]
        } 
        if (length(rx) != 1) stop("Check rx")
        death.rec <- rbinom(1, 1, rx)
        if (death.rec == 1) CHx[dead[t]] <- 1
        if (sur == 0) break 
      } else if (sur == 1) {
        # 3.) Is ind site faithful?: Function of age
        # ----------------------------------------------------------------------
        ## stays or leaves permanently ? --> Fixed to 1 for now
        stays <- rbinom(1,1, F)  
        if (stays == 0) break # check if expression is correct!
        ## 4.) re-sighted ?
        # ----------------------------------------------------------------------
        if (t > release) {
          px <- filter(p, fixed != "Fixed", time == (int[t]))[["estimate"]]
          if (length(px) != 1) stop("Check px")
          rp <- rbinom(1, 1, px)
          if (rp == 1) CHx[live[t]] <- 1  
        }
      }
    }  
    return(CHx)
  }) %>% 
    do.call("rbind",.)
  
  return(CH)
}

#' Simulate capture histories using a set of parameters
#' @param S data frame 
#' @param p data frame
#' @param F numeric value
#' @param r data frame
#' @param marked vector of releases per cohort
#' 
## simulate processes S, R, R' and r
simul.barker <- function(S, R, RPrime, r, marked, F = 1) {
  
  ## 1.) Set-up parameters
  ## -----------------------------------------------------------------------------
  int <- min(S$time):max(S$time) 
  oc <- length(int) 
  
  ## 2.) Initialise capture history matrix
  ## -----------------------------------------------------------------------------
  CH <- matrix(data = 0, ncol = 2 * oc, nrow = sum(marked) * oc) %>% 
    set_colnames(paste0(c("L", "D"), rep(int, each = 2)))
  
  # define a vector with marking occasion 
  mark.occ <- rep(seq(1, length(marked), 1), marked[1:length(marked)])
  
  # Functions to link to ch locations
  alive.slot <- function(i) i +  (-1 + i)
  dead.slot <- function(i) i +  (-2 + i)
  
  ## which columns are live encounters
  live <- sapply(1:oc, alive.slot)
  ## which columns are dead encounter
  dead <- live + 1
  
  ## Code releases in CH
  CH <- lapply(1:length(mark.occ), function(i) {
    CH[i, alive.slot(mark.occ[i])] <- 1
    return(CH[i,])
  }) %>% 
    do.call("rbind",.)
  
  ## 3.) Run simulation
  ## -----------------------------------------------------------------------------
  CH <- lapply(1:nrow(CH), function(ind) {
    CHx <- CH[ind,]
    ## column of the release
    col <- which(CHx == 1)
    ## corresponding mark.occ
    release <- which(live == col)
    
    ## Simulate fate of individual
    for (t in release:oc) { 
      ## get covariates to pick parameters 
      ## focal year
      time <- int[t]
      
      # 1.) Survival: Function of time and age-class 
      # ------------------------------------------------------------------------
      # Check if individual is within the young age groups
      # change if more than 3 age groups are used!
      if ((t - release) < 2) {
        Sx <- filter(S, time == (int[t]), age == (t - release))[["estimate"]]  
      } else {
        Sx <- filter(S, time == (int[t]), age >= (t - release))[["estimate"]]  
      }
      if (length(Sx) != 1) stop("Check Sx")
      
      sur <- rbinom(1, 1 , Sx)
      if (sur == 0) {
        # 2.) If dead, recorded or not: Functions of two age classes
        # ----------------------------------------------------------------------
        if ((t - release) == 0) {
          rx <- filter(r, age == 0)[["estimate"]]
        } else if ((t - release) >= 1) {
          rx <- filter(r, age == 1)[["estimate"]]
        } 
        if (length(rx) != 1) stop("Check rx")
        death.rec <- rbinom(1, 1, rx)
        if (death.rec == 1) CHx[dead[t]] <- 1
        if (sur == 0) break 
      } else if (sur == 1) {
        # 3.) Is ind site faithful?: Function of age
        # ----------------------------------------------------------------------
        ## stays or leaves permanently ? --> Fixed to 1 for now
        stays <- rbinom(1,1, F)  
        if (stays == 0) break # check if expression is correct!
        ## 4.) re-sighted ?
        # ----------------------------------------------------------------------
        if (t > release) {
          px <- filter(R, time == (int[t]))[["estimate"]]
          if (length(px) != 1) stop("Check px")
          rp <- rbinom(1, 1, px)
          if (rp == 1) CHx[dead[t]] <- 2  
        }
      }
    }  
    return(CHx)
  }) %>% 
    do.call("rbind",.)
  
  return(CH)
}

#' Convert LD to L 
#' @param ch ld capture history
#' @return character vector with life encounters
ld2l <- function(ch.ld = df$ch) {
  output <- sapply(ch.ld, function(ch.ld.x) {
    ch <- stringr::str_split(ch.ld.x, "")[[1]]
    ch <- paste0(ch[seq(1, length(ch),2)], collapse = "")
    return(ch)
  })
}


# start = 2007
# end = 2020
# resights = resighting_data
# x <- 3407451

#' @description Function that creates live-dead encounter histories
#' @param start first year of sampling period
#' @param last last year of sampling period
#' @param ringing_data data frame with ringing data of individuals
#' @param resights data frame with re-sighting data of individuals
create.ld.barker <- function(start = 2007, end = 2020, ringing_data = ringing_data, resights = resighting_data) {
  ## Identify all individuals 
  individuals <- ringing_data$Ring
  
  ## define Study periods based on start and end years
  years <- c(start, end)
  periods <- rep(seq(start, end), each = 1)
  
  ## create matrix 
  ## rows: Individuals
  ## cols: 2 Periods per year with Live and Dead encounters respectively
  mat <- matrix(0, nrow = length(individuals), ncol = length(periods)*2) %>% 
    set_colnames(paste0(rep(periods, each = 2), c("L", "D"))) %>% 
    set_rownames(individuals)
  
  ## code year of capture as 1,
  ## Code sightings in subsequent intervals as 1 if within study_area too
  for (x in individuals) {
    year.ringed <- filter(ringing_data, Ring == x)[["Year"]][1]
    ## code birth year in living encounters
    mat[rownames(mat) == x, colnames(mat) == paste0(year.ringed, "L")] <- 1  
    
    ## Distinguish 'recaptures' (study area and breeding period) from
    ## additional resightings outside study area and/or breeding period
    
    ## Sampling intervals follow the breeding cycle (May - April).
    ## Hence, spring month count for the previous year, i.e. March 2021 -> 2020
    ## For Barker model write encounter in D bins as a 2
    year.resighted <- 
      filter(resights, Ring == x)[,c("Year", "Month", "study_area")] %>% 
      mutate(period = case_when(
        Month == "01" ~ paste0(as.numeric(Year) - 1 , "D"),
        Month == "02" ~ paste0(as.numeric(Year) - 1 , "D"),
        Month == "03" ~ paste0(as.numeric(Year) - 1 , "D"),
        Month == "04" & study_area == FALSE ~ paste0(as.numeric(Year) - 1 , "D"),
        Month == "04" & study_area == TRUE ~ paste0(as.numeric(Year) - 1, "D"),
        Month == "05" & study_area == FALSE ~ paste0(Year, "D"),
        Month == "05" & study_area == TRUE ~ paste0(Year, "D"),
        Month == "06" & study_area == FALSE ~ paste0(Year, "D"),
        Month == "06" & study_area == TRUE ~ paste0(Year, "D"),
        Month == "06" ~ paste0(Year, "D"),
        Month == "07" ~ paste0(Year, "D"),
        Month == "08" ~ paste0(Year, "D"),
        Month == "09" ~ paste0(Year, "D"),
        Month == "10" ~ paste0(Year, "D"),
        Month == "11" ~ paste0(Year, "D"),
        Month == "12" ~ paste0(Year, "D")
      ))
    
    ## write live captures or re-sights to matrix 
    for (col in year.resighted$period) {
      mat[rownames(mat) == x, colnames(mat) == col] <- 
        ifelse(stringr::str_detect(col, "L"), 1, 2)
    }
    
    ## check if reported as dead
    death <- filter(ringing_data, Ring == x)[,c("Dead", "DateDeath")][1,]
    if (!is.na(death$Dead)) {
      death <- death %>% 
        mutate(Month = lubridate::month(death$DateDeath)) %>% 
        mutate(Year = lubridate::year(death$DateDeath)) %>% 
        mutate(period = case_when(
          Month == "1" ~ paste0(as.numeric(Year) - 1 , "D"),
          Month == "2" ~ paste0(as.numeric(Year) - 1 , "D"),
          Month == "3" ~ paste0(as.numeric(Year) - 1 , "D"),
          Month == "4" ~ paste0(as.numeric(Year) - 1 , "D"),
          Month == "5" ~ paste0(Year, "D"),
          Month == "6" ~ paste0(Year, "D"),
          Month == "7" ~ paste0(Year, "D"),
          Month == "8" ~ paste0(Year, "D"),
          Month == "9" ~ paste0(Year, "D"),
          Month == "10" ~ paste0(Year, "D"),
          Month == "11" ~ paste0(Year, "D"),
          Month == "12" ~ paste0(Year, "D")
        ))
      col <- death[["period"]]
      mat[rownames(mat) == x, colnames(mat) == col] <- 1
    }
    
  }
  
  ## coerce to character vector
  return(apply(mat[,], 1, function(x) paste0(x, collapse = "")) %>% 
           as.character())
  
}

#' Create corhors
#' @source Section 12 in RMark Workshop notes (modified)
#' @param x data frame
create.cohort <- function(x) {
  ## split ch string and convert to a matrix
  chmat <- sapply(x$ch,strsplit,split="") %>% 
    do.call("rbind",.)
  ## obtain releases
  releases <- chmat[,seq(1, ncol(chmat), 2)]
  ## return a factor based on each release ocassion
  return(factor(apply(releases, 1, function(x) min(which(x!="0"))))) 
}
