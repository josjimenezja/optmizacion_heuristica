
f_rastrigin2d<- function(x1,x2){
  y <- 20+x1^2-10*cos(2*pi*x1)+x2^2-10*cos(2*pi*x2)
  return(y)
}

n<-50
x1 <- seq(-5.12,5.12,length.out=n)
x2 <- seq(-5.12,5.12,length.out=n)

X <- expand.grid(x1,x2)
z <- f_rastrigin2d(X[,1],X[,2])
Z <- matrix(z,ncol=n,nrow=n)

#Curvas de nivel:
contour(x1,x2,Z)
filled.contour(x1, x2, Z, color.palette = bl2gr.colors)

#persp(x1,x2,Z,theta = 50, d=2)
persp3D(x1,x2,Z,theta = 50, phi = 20, col.palette = bl2gr.colors)

library(plotly)
fig <- plot_ly(z = as.matrix(Z), type = "surface")
fig <- fig %>% add_surface()

#Definamos la funci�n vectorizada:
f_rastrigin2d_vec <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  y <- 20+x1^2-10*cos(2*pi*x1)+x2^2-10*cos(2*pi*x2)
  return(y)
} 

#optimicemos la funci�n con optim:

opt_ras01 <- optim(par=c(0.01,-0.01), fn = f_rastrigin2d_vec)
#"par", que representa al valor inicial desde el cual queremos comenzar la b�squeda
# counts - function es la cantidad de iteraciones del procedimiento

#veamos c�mo se optimiza esta funci�n con evoluci�n diferencial:

#install.packages('DEoptim')
library(DEoptim)

opt_ras01 <- DEoptim(fn=f_rastrigin2d_vec, lower=c(-5.12,-5.12), upper = c(5.12,5.12))


#algoritmos gen�ticos
#install.packages('GA')
library(GA)
library(foreach)
library(iterators)

#Modifiquemos la funci�n de Rastrigin para minimizarla usando un m�todo de maximizaci�n:

f_rastrigin2d_vec_inv <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  y <- 20+x1^2-10*cos(2*pi*x1)+x2^2-10*cos(2*pi*x2)
  return(-y)
}

opt_ras03 <- ga(type="real-valued", fitness = f_rastrigin2d_vec_inv, lower=c(-5.12,-5.12), upper = c(5.12,5.12), seed=5)
#maxiter=100 y popSize=50 por default

contour(x1,x2,Z)
lines(opt_ras03@population,type="p",pch=2,col="red",lwd=3)

summary(opt_ras03)
mejor_poblacion <- cbind(opt_ras03@population,opt_ras03@fitness)

#veamos c�mo el promedio, la mediana y el mejor valor de la funci�n objetivo (el fitness) 
#evoluciona en cada iteraci�n repecto a la poblaci�n vigente en la iteraci�n
plot(opt_ras03)

optimArgs = list(control = list(trace = 1))

# Almacenar cada poblaci�n y c�mo evoluciona en cada iteraci�n
#bibliograf�a
# https://cran.r-project.org/web/packages/GA/vignettes/GA.html


