## Objective Identification: Station data
# Kriging

#### Settings ####
# Load package
library(dplyr)
library(lubridate)
library(fields)
library(viridis)
library(sf)
library(bestNormalize)
library(mvtnorm)
library(sp)
library(maps)
library(maptools)

# Path of R-Data
data_path <- "/.../"

# Path for Kriging data
res_path <- "/.../"

#### Load data ####
# Name of file
file_name <- paste0(data_path, "station_preds_for_kriging.csv") 

# Names of the winter storms
mydata <- read.csv(file = file_name,
                   sep = ";",
                   header = T)

# Storm names as characters
mydata <- mydata %>% 
  mutate(storm = as.character(storm))

#### Initiation ####
# Storms
storm_vec <- unique(mydata[["storm"]])

# Labels
lab_vec <- sort(unique(mydata[["lab"]]))

# Vector of countries to plot
regions_vec <- c("austria", "belgium", "czech republic", "denmark", "france", 
                 "germany", "luxembourg",  "netherlands", "poland", "switzerland", "uk")

# Colors of labels
lab_cols <- c("lightgrey",
              "#808080", # gray
              "#FF0000", # red
              "#008000", # green
              "#0000FF", # blue
              "#FFD700") # gold

#### Functions #####
# We use the Matern covariance model for kriging
# my.matern is the univariate Matern covariance model
my.matern<-function(u.h, u.a, u.sigma, u.nu)
{
  u.h[u.h == 0] <- 1e-10
  num1 <- (u.sigma^2)*(2^(1-u.nu))/gamma(u.nu)
  num2 <- (u.h*u.a)^u.nu
  num3 <- besselK(x = (u.h*u.a), nu = u.nu)
  return(num1*num2*num3)
}

## Now creating a function for a 5-variate Matern model
# rho_par is a function to create a 5x5 non-negative definite matrix as a function of any random values
rho_par<-function(k)
{ 
  L1 <- matrix(c(1, k[1], k[2], k[3], k[4],
                 0, 1, k[5], k[6], k[7], 
                 0, 0, 1, k[8], k[9], 
                 0, 0, 0, 1, k[10], 
                 0, 0, 0, 0, 1),byrow = T,nrow = 5)
  un_norm <- L1%*%t(L1)
  dterm <- 1/sqrt(diag(un_norm))
  dmat <- matrix(0,nrow = 5,ncol = 5)
  diag(dmat) <- dterm
  return(dmat%*%un_norm%*%dmat)
}

# my_5v.matern is a 5-variate Matern model, for details, refer to 2012 JASA paper by Tatiyana V. Apanasovich, Marc G. Genton and Ying Sun
my_5v.matern<-function(h,a,nu,sigma,beta)
{
  dmat <- h
  
  a1 <- a[1]
  a2 <- a[2]
  a3 <- a[3]
  a4 <- a[4]
  a5 <- a[5]
  
  nu1 <- nu[1]
  nu2 <- nu[2]
  nu3 <- nu[3]
  nu4 <- nu[4]
  nu5 <- nu[5]
  
  a12 <- sqrt(((a1^2)+(a2^2))/2)
  a13 <- sqrt(((a1^2)+(a3^2))/2)
  a14 <- sqrt(((a1^2)+(a4^2))/2)
  a15 <- sqrt(((a1^2)+(a5^2))/2)
  a23 <- sqrt(((a2^2)+(a3^2))/2)
  a24 <- sqrt(((a2^2)+(a4^2))/2)
  a25 <- sqrt(((a2^2)+(a5^2))/2)
  a34 <- sqrt(((a3^2)+(a4^2))/2)
  a35 <- sqrt(((a3^2)+(a5^2))/2)
  a45 <- sqrt(((a4^2)+(a5^2))/2)
  
  
  nu12 <- mean(c(nu1,nu2))
  nu13 <- mean(c(nu1,nu3))
  nu14 <- mean(c(nu1,nu4))
  nu15 <- mean(c(nu1,nu5))
  nu23 <- mean(c(nu2,nu3))
  nu24 <- mean(c(nu2,nu4))
  nu25 <- mean(c(nu2,nu5))
  nu34 <- mean(c(nu3,nu4))
  nu35 <- mean(c(nu3,nu5))
  nu45 <- mean(c(nu4,nu5))
  
  sigma1 <- sigma[1]
  sigma2 <- sigma[2]
  sigma3 <- sigma[3]
  sigma4 <- sigma[4]
  sigma5 <- sigma[5]
  
  Rvmat <- rho_par(beta)
  
  sigma12 <- sigma1*sigma2*Rvmat[1,2]*gamma(nu12)*(a1^nu1)*(a2^nu2)/(sqrt(gamma(nu1)*gamma(nu2))*(a12^(2*nu12)))
  sigma13 <- sigma1*sigma3*Rvmat[1,3]*gamma(nu13)*(a1^nu1)*(a3^nu3)/(sqrt(gamma(nu1)*gamma(nu3))*(a13^(2*nu13)))
  sigma14 <- sigma1*sigma4*Rvmat[1,4]*gamma(nu14)*(a1^nu1)*(a4^nu4)/(sqrt(gamma(nu1)*gamma(nu4))*(a14^(2*nu14)))
  sigma15 <- sigma1*sigma5*Rvmat[1,5]*gamma(nu15)*(a1^nu1)*(a5^nu5)/(sqrt(gamma(nu1)*gamma(nu5))*(a15^(2*nu15)))
  sigma23 <- sigma2*sigma3*Rvmat[2,3]*gamma(nu23)*(a2^nu2)*(a3^nu3)/(sqrt(gamma(nu2)*gamma(nu3))*(a23^(2*nu23)))
  sigma24 <- sigma2*sigma4*Rvmat[2,4]*gamma(nu24)*(a2^nu2)*(a4^nu4)/(sqrt(gamma(nu2)*gamma(nu4))*(a24^(2*nu24)))
  sigma25 <- sigma2*sigma5*Rvmat[2,5]*gamma(nu25)*(a2^nu2)*(a5^nu5)/(sqrt(gamma(nu2)*gamma(nu5))*(a25^(2*nu25)))
  sigma34 <- sigma3*sigma4*Rvmat[3,4]*gamma(nu34)*(a3^nu3)*(a4^nu4)/(sqrt(gamma(nu3)*gamma(nu4))*(a34^(2*nu34)))
  sigma35 <- sigma3*sigma5*Rvmat[3,5]*gamma(nu35)*(a3^nu3)*(a5^nu5)/(sqrt(gamma(nu3)*gamma(nu5))*(a35^(2*nu35)))
  sigma45 <- sigma4*sigma5*Rvmat[4,5]*gamma(nu45)*(a4^nu4)*(a5^nu5)/(sqrt(gamma(nu4)*gamma(nu5))*(a45^(2*nu45)))
  
  
  C11 <- my.matern(u.h = dmat,u.a = a1,u.sigma = sigma1,u.nu = nu1)
  C22 <- my.matern(u.h = dmat,u.a = a2,u.sigma = sigma2,u.nu = nu2)
  C33 <- my.matern(u.h = dmat,u.a = a3,u.sigma = sigma3,u.nu = nu3)
  C44 <- my.matern(u.h = dmat,u.a = a4,u.sigma = sigma4,u.nu = nu4)
  C55 <- my.matern(u.h = dmat,u.a = a5,u.sigma = sigma5,u.nu = nu5)
  
  
  C12 <- C21 <- sigma12*my.matern(u.h = dmat,u.a = a12,u.sigma = 1,u.nu = nu12)
  C13 <- C31 <- sigma13*my.matern(u.h = dmat,u.a = a13,u.sigma = 1,u.nu = nu13)
  C14 <- C41 <- sigma14*my.matern(u.h = dmat,u.a = a14,u.sigma = 1,u.nu = nu14)
  C15 <- C51 <- sigma15*my.matern(u.h = dmat,u.a = a15,u.sigma = 1,u.nu = nu15)
  C23 <- C32 <- sigma23*my.matern(u.h = dmat,u.a = a23,u.sigma = 1,u.nu = nu23)
  C24 <- C42 <- sigma24*my.matern(u.h = dmat,u.a = a24,u.sigma = 1,u.nu = nu24)
  C25 <- C52 <- sigma25*my.matern(u.h = dmat,u.a = a25,u.sigma = 1,u.nu = nu25)
  C34 <- C43 <- sigma34*my.matern(u.h = dmat,u.a = a34,u.sigma = 1,u.nu = nu34)
  C35 <- C53 <- sigma35*my.matern(u.h = dmat,u.a = a35,u.sigma = 1,u.nu = nu35)
  C45 <- C54 <- sigma45*my.matern(u.h = dmat,u.a = a45,u.sigma = 1,u.nu = nu45)
  
  C <- rbind(cbind(C11,C12,C13,C14,C15),
             cbind(C21,C22,C23,C24,C25),
             cbind(C31,C32,C33,C34,C35),
             cbind(C41,C42,C43,C44,C45),
             cbind(C51,C52,C53,C54,C55))
  
  return(C)
  
}

#### For-Loop over storms ####
# For-Loop over storms
for(temp_storm in storm_vec){
  #### Initiation ####
  # Console
  print(paste0("Storm: ", temp_storm))
  
  # Get subset of storm
  mydata_storm <- subset(mydata, storm == temp_storm)
  
  # Get time steps
  time_vec <- unique(mydata_storm[["time"]])
  
  #### For-Loop over time steps ####
  for(temp_step in time_vec){
    #### Initiation ####
    # Console
    print(paste0("Storm: ", temp_storm, "; Time: ", temp_step))
    
    # Get subset of time step
    l1 <- subset(mydata_storm, time == temp_step)
    
    #### Functions and preparation ####
    # Data
    t.data <- data.frame(lon = l1$lon,
                         lat = l1$lat,
                         v1 = l1$p0,
                         v2 = l1$p1,
                         v3 = l1$p2,
                         v4 = l1$p3,
                         v5 = l1$p5)
    
    # Training function
    my.uni.train <- function(t.uni.data, myseed = 123)
    {
      lon <- t.uni.data$lon
      lat <- t.uni.data$lat
      mydist <- rdist(cbind(lon,lat))
      
      
      ### Taking bestnormalize Transform ####
      colnames(t.uni.data)[3] <- 'v1'
      bnormalize <- bestNormalize(c(t.uni.data$v1))
      
      tp <- predict(bnormalize)
      
      train.data <- matrix(tp,ncol = 1,byrow = F)
      
      ####### Computing componentwise means ########
      l.v1.m <- mean(train.data[,1])
      
      mean.mat <- matrix(rep(c(l.v1.m),each = nrow(train.data)),ncol = 1,byrow = F)
      
      s.train.data <- train.data-mean.mat
      
      myz <- c(s.train.data[,1])
      
      mle.comp.all <- function(spdist,z,p)
      {
        
        
        
        
        if(sum(c(p[1:3] <= 0)) != 0)
        {
          return(list(mlv = Inf))
        }
        else
        {
          
          
          C <- my.matern(u.h = spdist,u.a = p[1],u.nu = p[2],u.sigma = p[3])
          nlogl <-- dmvnorm(x = z,mean = (rep(0,times = length(z))),sigma = C,log = T)
          
          return(list(mlv = nlogl))
        }
      }
      
      mle.comp_only_mlv <- function(par)
      {
        return(mle.comp.all(spdist = mydist,z = myz,p = par)$mlv
        )
      }
      optim_marg_comp.loglik <- function(par){
        optim(par = par,
              fn = mle.comp_only_mlv,
              hessian = FALSE,
              control = list(trace = 0, # No output, before 6
                             pgtol = 0,
                             parscale = rep(0.1,length(par)),
                             maxit = 5000))
      }
      set.seed(myseed)
      init.v <- c(runif(3,min = 0.1,max = 2))
      matern.mle <- optim_marg_comp.loglik(init.v)
      
      return(list(original.data = t.uni.data,
                  transform = bnormalize,
                  mean.mat = mean.mat,
                  matern.mle = matern.mle,
                  train.points = train.data,
                  s.train.data = s.train.data))
    }
    
    # Making a function that assigns a state, country to a given coordinate
    lonlat_to_state_sp <- function(pointsDF) {
      # Prepare SpatialPolygons object with one SpatialPolygon
      # per state (plus DC, minus HI & AK)
      states <- map('world',
                    regions = regions_vec, 
                    fill = TRUE, col = "transparent", 
                    plot = FALSE)
      IDs <- sapply(strsplit(states$names, ":"), function(x) x[1])
      states_sp <- map2SpatialPolygons(states, 
                                       IDs = IDs,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))
      
      # Convert pointsDF to a SpatialPoints object 
      pointsSP <- SpatialPoints(pointsDF, 
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
      
      # Use 'over' to get _indices_ of the Polygons object containing each point 
      indices <- over(pointsSP, states_sp)
      
      # Return the state names of the Polygons object containing each point
      stateNames <- sapply(states_sp@polygons, function(x) x@ID)
      stateNames[indices]
    }
    
    ## Now preparing prediction grid ##
    xmin <- min(mydata$lon) - 10
    xmax <- max(mydata$lon) + 10
    
    ymin <- min(mydata$lat) - 10
    ymax <- max(mydata$lat) + 20
    
    x.seq <- seq(xmin, xmax, length = 150)
    y.seq <- seq(ymin, ymax, length = 150)
    
    pred.grid <- expand.grid(x.seq, y.seq)
    colnames(pred.grid) <- c('x', 'y')
    check.points <- lonlat_to_state_sp(pred.grid)
    fin.pgrid <- pred.grid[!is.na(check.points),]
    
    # Throw out points too far South
    fin.pgrid <- fin.pgrid[(fin.pgrid[["x"]] <= 20) & 
                             (fin.pgrid[["x"]] >= -5) & 
                             (fin.pgrid[["y"]] <= 60) & 
                             (fin.pgrid[["y"]] >= 43),]

    # Making a function for univariate prediction
    uni.pred <- function(pred.points = fin.pgrid, mytrain_object)
    {
      train.locs <- cbind(mytrain_object$original.data$lon,mytrain_object$original.data$lat)
      pred.locs <- cbind(pred.points$x,pred.points$y)
      
      full.locs <- rbind(pred.locs,train.locs)
      prdl <- length(pred.locs[,1])
      
      full.dist <- rdist(full.locs)
      
      full.cov <- my.matern(u.h = full.dist,u.a  = mytrain_object$matern.mle$par[1],
                            u.nu = mytrain_object$matern.mle$par[2],
                            u.sigma = mytrain_object$matern.mle$par[3])
      
      cross.cov <- full.cov[1:prdl,-(1:prdl)]
      train.cov <- full.cov[-(1:prdl),-(1:prdl)]
      
      wts <- cross.cov%*%solve(train.cov)
      
      pred.vals <- wts%*%matrix(mytrain_object$s.train.data,ncol = 1)
      
      pred.vals2 <- pred.vals+mytrain_object$mean.mat[1]
      
      pred.vals3 <- predict(mytrain_object$transform,newdata = pred.vals2,inverse = T)
      
      # Values in (0, 1)
      pred.vals3 <- pmax(0, pred.vals3)
      
      return(pred.vals3)
    }
    
    #### Kriging ######
    
    t.uni.data.p0 <- t.data[,c(1,2,3)]
    t.uni.data.p1 <- t.data[,c(1,2,4)]
    t.uni.data.p2 <- t.data[,c(1,2,5)]
    t.uni.data.p3 <- t.data[,c(1,2,6)]
    t.uni.data.p5 <- t.data[,c(1,2,7)]
    
    v1.mle <- my.uni.train(t.uni.data.p0)
    v2.mle <- my.uni.train(t.uni.data.p1, myseed = 111)
    v3.mle <- my.uni.train(t.uni.data.p2)
    v4.mle <- my.uni.train(t.uni.data.p3)
    v5.mle <- my.uni.train(t.uni.data.p5)
    
    v1.pred <- uni.pred(mytrain_object = v1.mle)
    v2.pred <- uni.pred(mytrain_object = v2.mle)
    v3.pred <- uni.pred(mytrain_object = v3.mle)
    v4.pred <- uni.pred(mytrain_object = v4.mle)
    v5.pred <- uni.pred(mytrain_object = v5.mle)
    
    #### Normalization ####
    # Normalizing factor
    norm.factor <- v1.pred + v2.pred + v3.pred + v4.pred + v5.pred
    
    # No normalization for low accumulated probabilities (maybe incl. in max. prob.)
    norm.factor[norm.factor < 0.2] <- 1 
    
    # No normalization for low maximum probability
    norm.factor[pmax(v1.pred, v2.pred, v3.pred, v4.pred, v5.pred) < 0.2] <- 1 
    
    # Normalize
    v1.pred.n <- v1.pred/norm.factor
    v2.pred.n <- v2.pred/norm.factor
    v3.pred.n <- v3.pred/norm.factor
    v4.pred.n <- v4.pred/norm.factor
    v5.pred.n <- v5.pred/norm.factor
    
    # Only grid
    full.normal <- data.frame(lon = fin.pgrid[,1],
                              lat = fin.pgrid[,2],
                              P0 = v1.pred.n,
                              P1 = v2.pred.n,
                              P2 = v3.pred.n,
                              P3 = v4.pred.n,
                              P5 = v5.pred.n)
    
    # Get label of maximum probability
    norm.pred.lable <- apply(full.normal[,-c(1,2)], 1, "which.max")
    
    # Get maximum probabilities
    norm.adj.lable <- numeric()
    for(i in 1:length(full.normal$lon))
    {
      norm.adj.lable[i] <- full.normal[i,norm.pred.lable[i]+2]
    }
    
    # If corresponding probability smaller than 0.2, assign no feature 
    norm.lab.adf.fac <- which(norm.adj.lable < 0.2)
    norm.pred.lable[norm.lab.adf.fac] <-- 0
    
    # Assign to data frame
    full.normal$label <- norm.pred.lable
    
    # Transform labels
    full.normal[["label"]] <- full.normal[["label"]] - (full.normal[["label"]] < 5)
    
    # Write csv
    write.csv(full.normal,paste0(res_path, temp_storm, "_", temp_step, ".csv"))
  }
}