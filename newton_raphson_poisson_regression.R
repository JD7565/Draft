#Author: Juan Diego Mejia
#e-main: judmejiabe@unal.edu.co
#The following code provides an analysis for the number of born livestock       

rm(list=ls())
setwd("C:/Users/Usuario/Documents/DatosGanado")

library(ggplot2)

NewtonRaphsonForPoissonRegression <- function(y, X, initial.beta){
#This function estimates the parameters in a Poisson regression by means of the
#Newton Raphson algorithm
#@param: y="Response variable"
#@param: X="Design matrix"
#@param: initial.beta="initial value of the algorithm"
#@return: estimated.beta="Estimated parameters"
#@return: var.beta="Approximated variance matrix of the estimated parameters"
#@return: deviance"Calculated deviance"
  eta <- function(y, X, beta){
  #This function calculates the systematic component of the model
  #@return: systematic component of the model
    return (X %*% beta)
  }
  mu <- function(y, X, beta){
  #This function calculates the vector of means
    return (exp(eta(y=y, X=X, beta=beta)))
  }
  dmu <- function(y, X, beta){
  #This function may be changed in order to change the link function
  #If changing link function, please change d2mu
    return (exp(eta(y=y, X=X, beta=beta))) 
  }
  d2mu <- function(y, X, beta){
    return (exp(eta(y=y, X=X, beta=beta))) 
  }
  var.y <- function(y, X, beta){
    return (exp(eta(y=y, X=X, beta=beta))) 
  }
  D <- function(y, X, beta){
  #This function returns a matrix useful for reducing calculations in the NR
  #algorithm
    return (diag(as.numeric(dmu(y, X, beta)/mu(y, X, beta))))
  }
  u <- function(y, X, beta){
    return (t(X) %*% D(y, X, beta) %*% (y - mu(y, X, beta)))
  }
  V <- function(y, X, beta){
    return (diag(as.numeric((y-mu(y, X, beta))/mu(y, X, beta)*d2mu(y, X, beta)-
      y/(mu(y, X, beta)^2)*dmu(y, X, beta)^2)))
  }
  H <- function(y, X, beta){
    return (t(X) %*% V(y, X, beta) %*% X)
  }
  W <- function(y, X, beta){
    return (diag(as.numeric(dmu(y, X, beta)^2/var.y(y, X, beta))))
  }
  var.beta <- function(y, X, beta){
    return(solve(t(X) %*% W(y, X, beta) %*% X))
  }
  VarMatrixY <- function(y, X, beta){
    return(diag(as.numeric(mu(y, X, beta))))
  }
  HatMatrix <- function(y, X, beta){
    return(sqrt(W(y, X, beta)) %*% X %*% var.beta(y, X, beta) %*% 
      t(X) %*% sqrt(W(y, X, beta)))
  }
  VarOfResiduals <- function(y, X, beta){
    return(sqrt(VarMatrixY(y, X, beta)) %*% 
      (diag(1, length(y)) - HatMatrix(y, X, beta)) %*% 
      sqrt(VarMatrixY(y, X, beta)))
  }
  StandarizedResiduals <- function(y, X, beta){
    diagonal.var.residuals <- diag(VarOfResiduals(y, X, beta))
    return((y-mu(y, X, beta))/sqrt(diagonal.var.residuals))
  }
  tol<-1
  previous.beta <- initial.beta
  while (tol>10^(-6)){
    beta <- previous.beta - solve(H(y, X, previous.beta)) %*% 
      u(y, X, previous.beta)
    tol <- max(abs((beta-previous.beta)/previous.beta))
    previous.beta <- beta
  }
  Deviance <- function(y, X, beta){
    d <- rep(0, length(y))
    for (i in 1:length(y)){
      if (y[i] == 0){
        d[i] <- 2 * mu(y, X, beta)[i]
      } else{
        d[i] <- 2* ( y[i] * log(y[i]/mu(y, X, beta)[i]) + 
          mu(y, X, beta)[i] - y[i])
      }
    }
    return(sum(d))
  }
  colnames(beta) <- c("estimated.value")
  return (list(beta = beta, var.beta = var.beta(y, X, beta), deviance =
    Deviance(y, X, beta), prediction = mu(y, X, beta),
    standarized.residuals=StandarizedResiduals(y, X, beta)))
}

TestNewtonRaphsonForPoissonRegression <- function(){
  data <- read.csv(file="data.csv", sep=";")
  y <- as.matrix(data$born)
  ones <- rep(1, 54)
  X <- cbind(ones, data$precipitation, data$rainy.days, 
    data$livestock.at.reproductive.age)
  X <- as.matrix(X)
  colnames(X) <- c("constant.term","precipitation","rainy.days","livestock")
  initial.beta <- as.matrix(c(0, 0, 0, 0))
  print(NewtonRaphsonForPoissonRegression(y, X, initial.beta))
}

ModelSelection <- function(y, X, variable.names){
  #variable.names must only contain the names of the variables except the 
  #intercept
  #@return: important variables
  output <- NewtonRaphsonForPoissonRegression(y, X, rep(0, dim(X)[2]) )
  estimated.sd <- sqrt(diag(output$var.beta)[2:dim(X)[2]])
  t.scores <- output$beta[2:dim(X)[2]]/estimated.sd
  selected.variables <- variable.names
  #on.vector is the vector of indexes of variables in the model
  on.vector <- seq(1, dim(X)[2])

  while (TRUE){
    #1 must be added to proposal.off since the constant term of the design
    #matrix X occupies the first place
    proposal.off <- which.min(t.scores[(on.vector[-1]-1)])+1
    deviance.null <- NewtonRaphsonForPoissonRegression(y, X[ ,
      on.vector[-proposal.off]], 
      as.matrix(rep(0, length(on.vector[-proposal.off]))))$deviance
    deviance.alt <- NewtonRaphsonForPoissonRegression(y, X[ , 
      on.vector], as.matrix(rep(0, length(on.vector))))$deviance
    if (deviance.null-deviance.alt > qchisq(0.95, 1)){
      break
    }else{
      on.vector <- on.vector[-proposal.off]
      if ( length(on.vector) == 1 ){
        print("No variables seem to be important")
        break
      }
    }
  }
  return (selected.variables[on.vector[-1]])
}

TestModelSelection <- function(){
  data <- read.csv(file="data.csv", sep=";")
  y <- as.matrix(data$born)
  ones <- rep(1, 54)
  X <- cbind(ones, data$precipitation, data$rainy.days, 
    data$livestock.at.reproductive.age)
  X <- as.matrix(X)
  colnames(X) <- variable.names <- c("constant.term",
    "precipitation","rainy.days","livestock")
  print(ModelSelection(y, X, variable.names))
}

OutlierDetection <- function(y, X, confidence.level = 0.95){
#Detects outlier using likelihood ratio tests
  n <- length(y)
  p <- dim(X)[2]
  outliers <- NULL
  #Current is a model with all data
  current <- NewtonRaphsonForPoissonRegression(y, X, 
      initial.beta = rep(0, p))
  deviance.current <- current$deviance
  beta.current <- current$beta
  for ( i in 1:n ){
    new.model <- NewtonRaphsonForPoissonRegression(y[-i], X[-i, ], 
      initial.beta = rep(0, p))
    if (deviance.current - new.model$deviance > qchisq(confidence.level, 1)){
      outliers <- c(outliers, i)
    }
  }
  if ( length(outliers) == 0 ){
    print("No outliers were found")
  } else{
    return(outliers)
  }
}

TestOutlierDetection <- function(){
  data <- read.csv(file="data.csv", sep=";")
  y <- as.matrix(data$born)
  ones <- rep(1, 54)
  X <- cbind(ones, data$precipitation, data$rainy.days, 
    data$livestock.at.reproductive.age)
  X <- as.matrix(X)
  colnames(X) <- variable.names <- c("constant.term",
    "precipitation","rainy.days","livestock")
  OutlierDetection(y, X)
}

##############################Analysis##############################

data <- read.csv(file="data.csv", sep=";")
n <- dim(data)[1]
p <- dim(data)[2]
y <- as.matrix(data$born)
ones <- rep(1, 54)
X <- cbind(ones, data$precipitation, data$rainy.days, 
  data$livestock.at.reproductive.age)
X <- as.matrix(X)
colnames(X) <- variable.names <- c("constant.term",
  "precipitation","rainy.days","livestock")
#livestock is the only important variable
#ModelSelection(y, X, variable.names)

model <- NewtonRaphsonForPoissonRegression(y, X[ ,c(1,4)], c(0,0))
outliers <- OutlierDetection(y, X[ ,c(1,4)])
outlier.column <- rep(0, n)
outlier.column[outliers] <- 1

data.for.plotting <- as.data.frame(cbind(data$born, X[ ,-1], model$prediction,
  model$standarized.residuals, outlier.column))
names (data.for.plotting)[1] <- "born"
names (data.for.plotting)[5] <- "estimated.mean"
names (data.for.plotting)[6] <- "standarized.residuals"
names (data.for.plotting)[7] <- "outlier"


ModelMean <- function(x){
  return ( exp(-0.617 + 0.0305*x) )
}

VarOfMean <- function(x){
  mu <- ModelMean(x)
  m <- length(x)
  variance.vector <- rep(0, m)
  model$var.beta
  for ( i in 1:m ){
    x.vector <- t(as.matrix(c(1,x[i])))
    variance.vector[i] <- mu[i]^2 * x.vector %*% model$var.beta %*% 
      t(x.vector)
  }
  return (variance.vector)
}


livestock.sequence <- seq(16, 46, length = 100)
mean.model <- ModelMean(livestock.sequence)
var.of.estimated.mean <- VarOfMean(livestock.sequence)

livestock.sequence <- as.data.frame(livestock.sequence)
names(livestock.sequence) <- "livestock"
lower.confidence.band <- mean.model - 1.95 * sqrt(var.of.estimated.mean)
upper.confidence.band <- mean.model + 1.95 * sqrt(var.of.estimated.mean)

ggplot(livestock.sequence, aes(livestock)) +
  geom_line(aes(y = mean.model), colour = "blue") + 
  geom_ribbon(aes(ymin = lower.confidence.band, ymax = upper.confidence.band),
    alpha = 0.2)+
  geom_point(data = data.for.plotting, aes(x = livestock, y = born, 
    colour = factor (outlier))) +
  scale_colour_manual(values = c("black", "#FF3333"), breaks = c("1", "0"), 
    labels = c("Si", "No"), name = "Outlier") +
  ggtitle("Predicción de Nacimientos \n") + 
  xlab("Población de Ganado en Edad Reproductiva") + 
  ylab("Número de Nacimientos") + 
  annotate("text", x = 46, y = 0.2, label = "Mayo 2016", colour = "red", 
    size = 2) + 
  annotate("text", x = 27.5, y = 2, 
    label = "mu(x)==italic(e)^{-0.617 + 0.0305*x}", 
    parse = TRUE, size = 4) + 
  annotate("text", x = 20, y = 4.5, 
    label = "Distribución: Poisson(mu(x))", parse = TRUE, size = 5 )

ggplot(data.for.plotting, aes(x = factor(0), standarized.residuals)) +
  geom_boxplot()

ggplot(data.for.plotting, aes(x = estimated.mean, y = standarized.residuals,
  colour = factor(outlier))) + 
  scale_colour_manual(values = c("black", "#FF3333"), breaks = c("1", "0"), 
    labels = c("Si", "No"), name = "Outlier") +
  geom_point() +
  geom_smooth(colour = "blue")

ggplot(data.for.plotting, aes(x = livestock, y = standarized.residuals,
  colour = factor(outlier))) + 
  scale_colour_manual(values = c("black", "#FF3333"), breaks = c("1", "0"), 
    labels = c("Si", "No"), name = "Outlier") +
  geom_point() +
  geom_smooth(colour = "blue")

ggplot(data.for.plotting, aes(x = precipitation, y = standarized.residuals,
  colour = factor(outlier))) + 
  scale_colour_manual(values = c("black", "#FF3333"), breaks = c("1", "0"), 
    labels = c("Si", "No"), name = "Outlier") +
  geom_point() +
  geom_smooth(colour = "blue")

ggplot(data.for.plotting, aes(x = rainy.days, y = standarized.residuals,
  colour = factor(outlier))) + 
  scale_colour_manual(values = c("black", "#FF3333"), breaks = c("1", "0"), 
    labels = c("Si", "No"), name = "Outlier") +
  geom_point() +
  geom_smooth(colour = "blue")

standarized.residuals.ts <- as.data.frame(cbind(seq(1, n), 
  data.for.plotting$standarized.residuals))
names(standarized.residuals.ts) = c("observation", "standarized.residuals")
qplot(observation, standarized.residuals, data = standarized.residuals.ts, 
  geom = "line")

acf(data.for.plotting$standarized.residuals)
