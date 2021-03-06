#install.packages('GA')
#install.packages('DEoptim')
#install.packages("pso")
#install.packages("rgl")
library(ggplot2)
library(DEoptim)
library(GA)
library(pso)
library(rgl)
library(plotly)

# ----------------PUNTO 1: Escoja dos funciones de prueba

#Escogemos Rastrigin y Schwefel:


#determinamos la funci�n de rastrigin a utilizar para graficar:

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

persp(x1,x2,Z,theta = 50, d=2)
persp3d(x1,x2,Z,theta = 50, phi = 20, col.palette = bl2gr.colors)

fig <- plot_ly(z = as.matrix(Z), type = "surface")
fig <- fig %>% add_surface()
fig

#rastrigin vectorizada: para optimizar
 
f_rastrigin2d_vec <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  y <- 20+x1^2-10*cos(2*pi*x1)+x2^2-10*cos(2*pi*x2)
  return(y)
} 

# --------HACEMOS LO MISMO PERO PARA SCHWEFEL:

#determinamos la funci�n de schwefel a utilizar para graficar:

f_schwefel2d<- function(x1,x2){
  y <- (-((x1*sin(sqrt(abs(x1)))) + (x2*sin(sqrt(abs(x2)))) ))
  return(y)
}

n<-50
x1 <- seq(-500,500,length.out=n)
x2 <- seq(-500,500,length.out=n)

X <- expand.grid(x1,x2)
z <- f_schwefel2d(X[,1],X[,2])
Z <- matrix(z,ncol=n,nrow=n)

#Curvas de nivel:

contour(x1,x2,Z)
filled.contour(x1, x2, Z, color.palette = bl2gr.colors)

persp(x1,x2,Z,theta = 50, d=2)

#la persp3d Esta poniendo problema para graficar con colores
persp3d(x1,x2,Z,theta = 50, phi = 20, col.palette = bl2gr.colors)

fig <- plot_ly(z = as.matrix(Z), type = "surface")
fig <- fig %>% add_surface()
fig

#schwefel  vectorizada: para optimizar

f_schwefel2d_vec <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  y <- (-((x1*sin(sqrt(abs(x1)))) + (x2*sin(sqrt(abs(x2))))))
  return(y)
} 

# ----------- PUNTO 2: DESCENSO POR GRADIENTE

#calculamos la derivada parcial
partial_dev <- function(x,i,fun,h=0.01){
  e <- x*0 # crea un vector de ceros de la misma longitud de x
  e[i] <- h
  y <- (fun(x+e)-fun(x-e))/(2*h)
  return(y)
}

#obtenemos cada una de las derivadas parciales de f en x:

num_grad <- function(x,fun,h=0.01){
  # x: punto del espacio donde se debe evaluar el gradiente
  # fun: funci�n para la que se desea calcular el gradiente en x
  # h: es el tama�o de ventana para el c�lculo de la derivada num�rica
  d <- length(x)
  y <- mapply(FUN=partial_dev,i=1:d,MoreArgs=list(x=x,h=h,fun=fun))
  return(y)
}

# SIN PRECONDICIONAMIENTO:

#determinamos la funci�n para evaluar el gradiente de la funci�n Rastrigin:

deriv_grad <- function(x,fun,i=1,h=0.01){
  # x: punto en el que se eval�a el gradiente
  # fun: funci�n para la cual se calcula la derivada del gradiente respecto a la �esima componente
  # i: i-�sima componente del vector x con respecto a la que se deriva
  e <- x*0 # crea un vector de ceros de la misma longitud de x
  e[i] <- h
  y <- (num_grad(x+e,fun=fun,h=h)-num_grad(x-e,fun=fun,h=h))/(2*h)
  return(y)
}

# Calculamos la matriz hessiana:

matriz_hessiana <- function(x,fun,h=0.01){
  # x: punto en el que se eval�a la matriz hessiana
  # fun: funci�n a la que se le calcula la matriz hessiana en x
  # h: es el tama�o de ventana para el c�lculo de la derivada num�rica
  d <- length(x)
  y <- mapply(FUN=deriv_grad,i=1:d,MoreArgs=list(x=x,h=h,fun=fun),SIMPLIFY = TRUE)
  return(y)
}

# Determinamos la funci�n de optimizador multivariado

optimizador_mult_numdev <- function(x0,fun,max_eval=100,h=0.01,eta=0.01){
  x <- matrix(NA,ncol =length(x0), nrow = max_eval)
  x[1,] <- x0
  for (i in 2:max_eval){
    num_grad_fun <- num_grad(x[i-1,],fun,h)
    H <- matriz_hessiana(x[i-1,],fun,h)
    cambio <- - eta*solve(H)%*%num_grad_fun
    x[i,] <- x[i-1,] + cambio
    cambio_opt <- sqrt(sum((x[i-1,]-x[i,])^2))
    if (cambio_opt<0.00001){
      break
    }
  }
  return(x[1:i,])
}

#pasamos rastrigin al optimizador iniciando desde (0.8,4.8)

sol_ras <- optimizador_mult_numdev(f_rastrigin2d_vec,x0=c(0.8,4.8),eta=1)

contour(x1, x2, Z, las=1,
        xlab = expression(x[1]), ylab = expression(x[2]),
        main = "Funci�n de Rastrigin",
        sub = "Curvas de nivel de la funci�n")
lines(sol_ras, type="b",cex=1.5,col="red")


#CON PRECONDICIONAMIENTO

#hacemos el precondicionamiento para la cantidad de dimensiones

deriv_segunda <- function(x,fun,i=1,h=0.01){
  e <- x*0 # crea un vector de ceros de la misma longitud de x
  e[i] <- h
  y <- (fun(x+e)-2*fun(x)+fun(x-e))/(h^2)
  return(y)
}

# Determinamos la funci�n de optimizador multivariado

optimizador_mult_precond <- function(x0,fun,max_eval=100,h=0.01,eta=0.01){
  x <- matrix(NA,ncol =length(x0), nrow = max_eval)
  d <- length(x0)
  x[1,] <- x0
  for (i in 2:max_eval){
    num_grad_fun <- num_grad(x[i-1,],fun,h)
    diag_H <- mapply(FUN=deriv_segunda,i=1:d,MoreArgs=list(x=x[i-1,],h=h,fun=fun))
    # print(diag_H)
    # H <- matriz_hessiana(x[i-1,],fun,h)
    # print(H)
    H_precond <- diag(1/diag_H)
    cambio <- - eta*H_precond%*%num_grad_fun
    x[i,] <- x[i-1,] + cambio
    cambio_opt <- sqrt(sum((x[i-1,]-x[i,])^2))
    if (cambio_opt<0.0000001){
      break
    }
  }
  return(x[1:i,])
}

#pasamos rastrigin al optimizador iniciando desde (0.8,4.8)

sol_ras_precon <- optimizador_mult_precond(f_rastrigin2d_vec,x0=c(0.8,4.8),eta=1)

contour(x1, x2, Z, las=1,
        xlab = expression(x[1]), ylab = expression(x[2]),
        main = "Funci�n de Rastrigin",
        sub = "Curvas de nivel de la funci�n")
lines(sol_ras_precon, type="b",cex=1.5,col="red")
lines(sol_ras, type="b",cex=1.5,col="blue")
legend("topright",col=c("red","blue"),legend = c("Con precondicionamiento","Sin precondicionamiento"),lty=1)


#----------------------- AHORA CON FUNCION SCHWEFEL


#----SIN PRECONDICIONAMIENTO

#pasamos schwefel al optimizador SIN precondicionamiento iniciando desde (380,390)

sol_schwefel <- optimizador_mult_numdev(f_schwefel2d_vec,x0=c(380,390),eta=1)

contour(x1, x2, Z, las=1,
        xlab = expression(x[1]), ylab = expression(x[2]),
        main = "Funci�n de Schwefel",
        sub = "Curvas de nivel de la funci�n")
lines(sol_schwefel, type="b",cex=1.5,col="red")


#----CON PRECONDICIONAMIENTO

#pasamos schwefel al optimizador CON precondicionamiento iniciando desde (380,390)

sol_schwefel_precon <- optimizador_mult_precond(f_schwefel2d_vec,x0=c(380,390),eta=1)

contour(x1, x2, Z, las=1,
        xlab = expression(x[1]), ylab = expression(x[2]),
        main = "Funci�n de Schwefel",
        sub = "Curvas de nivel de la funci�n")
lines(sol_schwefel_precon, type="b",cex=1.5,col="red")
lines(sol_schwefel, type="b",cex=1.5,col="blue")
legend("topleft",col=c("red","blue"),legend = c("Con precondicionamiento","Sin precondicionamiento"),lty=1)




# -------PUNTO 3: algoritmos evolutivos, optimizaci�n de part�culas y evoluci�n diferencial

# RASTRIGIN

#------------------------ #optimizacion con evoluci�n diferencial:

opt_ras01 <- DEoptim(fn=f_rastrigin2d_vec, lower=c(-5.12,-5.12), upper = c(5.12,5.12))

#-------------------- optimizacion con algoritmos gen�ticos

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


#-------------------- optimizacion con particulas

opt_ras_pso <- psoptim(par=c(NA,NA), fn=f_rastrigin2d_vec, lower=c(-5.12,-5.12),
                             upper = c(5.12,5.12), control = list(trace.stats=TRUE,trace=1))

contour(x1,x2,Z,main="Optimizaci�n con part�culas")
lines(t(opt_ras_pso$stats$x[[7]]), type="p", pch=2,col="red",lwd=3)
#Ir cambiando el x[[1]] para ver c�mo se comporta



# ------SCHWEFEL


#------------------------ #optimizacion con evoluci�n diferencial:


opt_schwefel_DE <- DEoptim(fn=f_schwefel2d_vec, lower=c(-500,-500), upper = c(500,500))

#-------------------- optimizacion con algoritmos gen�ticos


#Modificamos la funci�n de Schwefel para minimizarla usando un m�todo de maximizaci�n:

f_schwefel2d_vec_inv <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  y <- (-((x1*sin(sqrt(abs(x1)))) + (x2*sin(sqrt(abs(x2)))) ))
  return(-y)
}

opt_schwefel_GA <- ga(type="real-valued", fitness = f_schwefel2d_vec_inv, lower=c(-500,-500), upper = c(500,500), seed=5)
#maxiter=100 y popSize=50 por default

contour(x1,x2,Z)
lines(opt_schwefel_GA@population,type="p",pch=2,col="red",lwd=3)

summary(opt_schwefel_GA)
mejor_poblacion <- cbind(opt_schwefel_GA@population,opt_schwefel_GA@fitness)

#veamos c�mo el promedio, la mediana y el mejor valor de la funci�n objetivo (el fitness) 
#evoluciona en cada iteraci�n repecto a la poblaci�n vigente en la iteraci�n
plot(opt_schwefel_GA)

optimArgs = list(control = list(trace = 1))


#-------------------- optimizacion con particulas


opt_schwefel_pso <- psoptim(par=c(NA,NA), fn=f_schwefel2d_vec, lower=c(-500,-500),
                            upper = c(500,500), control = list(trace.stats=TRUE,trace=1))

contour(x1,x2,Z)
lines(t(opt_schwefel_pso$stats$x[[1]]), type="p", pch=2,col="red",lwd=3)
#Ir cambiando el x[[1]] para ver c�mo se comporta


#bibliograf�a
# https://cran.r-project.org/web/packages/GA/vignettes/GA.html



