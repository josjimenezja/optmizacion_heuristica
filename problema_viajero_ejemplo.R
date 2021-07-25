## Usango algoritmos geneticos
library(GA)
library(tidyverse)
library(animation)
library(ggrepel)

# Creamos la matrix de coordenadas
ciudades <- c("Palmira",	"Pasto",	"Tuluá",	"Bogota",	"Pereira",	"Armenia",	"Manizales",	"Valledupar",	"Montería",	"Soledad",
              "Cartagena",	"Barranquilla",	"Medellín",	"Bucaramanga",	"Cúcuta")
coordenadas <- data.frame( long = c(-76.30361, -77.28111, -76.19536, -74.08175, -75.69611, -75.68111, -75.51738, -73.25322, -75.88143, -74.76459,
                                    -75.51444, -74.78132, -75.56359, -73.1198, -72.50782),
                           lat= c(3.53944, 1.21361, 4.08466, 4.60971, 4.81333, 4.53389, 5.06889, 10.46314, 8.74798, 10.91843,
                                  10.39972, 10.96854, 6.25184, 7.12539, 7.89391),
                           stringsAsFactors = F)
# Leemos la matrix de costos
D <- read.csv("costos_ciudades2.csv", sep = ";", row.name = 1)

#Function to calculate tour length

tourLength <- function(tour, distMatrix) {
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}

#Firness function to be maximized

tspFitness <- function(tour, ...) 1/tourLength(tour, ...)
# Iteracciones maximas
Maxter <- c(25, 50, 75, 100, 200, 400, 600, 800, 1000,  1200)

# Ciudades
x <- coordenadas[,1]
y <- coordenadas[,2]

# Utilizamos la funcion saveGIF para crear el gif durante la iteraccion de la funcion
# ga, ademas guardamos cada imagen creada
saveGIF(
  # For para cambiar el numero maximo de iteraciones en la funcion ga
  for (i in 1:length(Maxter)) {
    set.seed(124)

    GAiter <- ga(type = "permutation", fitness = tspFitness, distMatrix = D,
                 lower = 1, upper = 15, popSize = 50, maxiter = Maxter[i],
                 run = 500, pmutation = 0.2)

    tour <- GAiter@solution[1, ]
    tour <- c(tour, tour[1])
    n <- length(tour)

    ## Grafica
    map <- ggplot() +
      borders("world", "colombia") +
      geom_point(data = coordenadas,
                 aes(x = long,
                     y = lat), colour = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))

    map_vectors <- map +
      ggtitle(paste( paste("Max. numero de iteraciones:", Maxter[i]), paste(" - Iteracciones realizadas: ", GAiter@iter))) +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_segment(aes(x = x[tour[-n]], y = y[tour[-n]], xend = x[tour[-1]], yend = y[tour[-1]]),
                   arrow = arrow(length = unit(0.25, "cm")), col = "steelblue", lwd= 0.6) +
      geom_text_repel(aes(x, y, label = ciudades),
                      size = 3.5)
    plot(map_vectors)

    # Guardamos las imagenes para el GIF
    ggsave(map_vectors, filename = gsub(" ", "", paste("Ga_", i,".png")),  width = 5.78 , height = 4.66)

  }

  , movie.name = "TSP_GA.gif")



## Usango colonia de hormigas

# Creamos la matrix de coordenadas
#ciudades <- c("Palmira",	"Pasto",	"Tuluá",	"Bogota",	"Pereira",	"Armenia",	"Manizales",	"Valledupar",	"Montería",	"Soledad",
#              "Cartagena",	"Barranquilla",	"Medellín",	"Bucaramanga",	"Cúcuta")
#coordenadas <- data.frame( long = c(-76.30361, -77.28111, -76.19536, -74.08175, -75.69611, -75.68111, -75.51738, -73.25322, -75.88143, -74.76459,
#                                    -75.51444, -74.78132, -75.56359, -73.1198, -72.50782),
#                           lat= c(3.53944, 1.21361, 4.08466, 4.60971, 4.81333, 4.53389, 5.06889, 10.46314, 8.74798, 10.91843,
#                                  10.39972, 10.96854, 6.25184, 7.12539, 7.89391),
#                           stringsAsFactors = F)

# Leemos el archivo resultante luego de ejecutar el metodo de colonia de
# hormigas en python
Resultado_Ant <- read.csv("tours_results.csv", sep = ",", row.name = 1)
iter <- c(4,6,8,10,12,14,16,18,20,22)

# Ciudades
#x <- coordenadas[,1]
#y <- coordenadas[,2]

saveGIF(
  # For para cambiar los datos que vamos a graficar
  for (i in c(1,2,3,4,5,6,7,8,9,10)) {

    tour <- as.matrix(Resultado_Ant[i])
    tour <- tour + 1
    n <- length(tour)
    
    ## Grafica
    map <- ggplot() +
      borders("world", "colombia") +
      geom_point(data = coordenadas,
                 aes(x = long,
                     y = lat), colour = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
    
    map_vectors <- map +
      ggtitle(paste( paste("Numero de iteraciones:", iter[i]), "Hormigas: 100")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_segment(aes(x = x[tour[-n]], y = y[tour[-n]], xend = x[tour[-1]], yend = y[tour[-1]]),
                   arrow = arrow(length = unit(0.25, "cm")), col = "steelblue", lwd= 0.6) +
      geom_text_repel(aes(x, y, label = ciudades),
                      size = 3.5)
    plot(map_vectors)
    
    # Guardamos las imagenes para el GIF
    ggsave(map_vectors, filename = gsub(" ", "", paste("Ants_", i,".png")),  width = 5.78 , height = 4.66)
    
  }
  , movie.name = "TSP_ANTS.gif")


#Mejor resultado GA
MGA_1 <- c(11,  9, 13,  7,  5,  6,  3,  1,  2, 4, 14,  15,   8,  10,  12)
MGA_2 <- c(11, 12, 10,  8, 15, 14,  4,  2,  1,   3,   6,   5,   7,  13,   9)