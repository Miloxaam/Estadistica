
#                                   Actividad 4
#                                   EJERCICIO 1
#                                                          13 de septiembre 2019


# A continuación se encuentra la obtención numérica del jacobiano y se compara el resultado
# con el determinante de Vandermonde

# Datos
H <- replicate(n=10,{sample(rnorm(n=10, 0, 0.2))})    # Matriz con entradas iid N(0,1)
Hs <- (H+t(H))/2                                    # SIMETRÍA
Lambda.Hs <- eigen(Hs)
eigen.values <- Lambda.Hs$values                    # Matriz de eigenvalues
eigen.vectores <- Lambda.Hs$vectors                 # Matriz de eigenvectors
Jacob.matriz <- matrix(0,55,55)
epsilon <- exp(-7)
count <- 1

# Creación de Jacobbiano  
for (i in 1:10) {
  for (j in i:10) {
    
    E <- matrix(0,10,10)
    E[i,j] <- 1
    E[j,i] <- 1
    Hs.per <- Hs + epsilon*E
    
    
    Lambda.Hs.per <- eigen(Hs.per)
    eigen.values.p <- Lambda.Hs.per$values
    eigen.vectores.p <- Lambda.Hs.per$vectors
    dl <- (eigen.values.p - eigen.values)/epsilon
    QdQ <- (t(eigen.vectores)%*%(eigen.vectores.p-eigen.vectores))/epsilon
    QdQ <-QdQ[lower.tri(QdQ, diag = FALSE)]
    dl <- as.matrix(dl, 10,1)
    QdQ <- as.matrix(QdQ, 45, 1)
    jac.col <- rbind(dl,QdQ)
    
    Jacob.matriz[ ,count] <- jac.col
    count <- count+1
  }
}

# Determinante del Jacobiano
det.jac <- abs(det(Jacob.matriz))
det.jac
# Determinante de Vandermonde
library(matrixcalc)
Van.matriz <- vandermonde.matrix(eigen.values, 10)
Det.Vander <- abs(det(Van.matriz))
Det.Vander <- 1/Det.Vander
Det.Vander  

# EJERCICIO 2
library(elliptic)
# Se tiene como objetivo comprobar la consistencia del círculo de Wigner en
# la expresión del revolvente.

# Para ello se realizaron dos casos

f<-function(x, z){                          # Convertir la función de 2 variables a una
  function(x){
    sqrt(2-x^2)/z-x  
  }
}

integral.compleja1 <- function(z){            # Función que integra caso 1 lambda en (0,sqrt(2))
  fx <- f(x, z)
  myintegrate(fx,0,sqrt(2))  
}

integral.compleja2 <- function(z){            # Función que integra caso 2 lambda en (-sqrt(2), 0)
  fx <- f(x, z)
  myintegrate(fx,-sqrt(2), 0)  
}

# Caso 1 lambda en (0,sqrt(2))
z1 <- complex(length.out = 10, real = sample(runif(n = 10, min = 0, max = sqrt(2))), 
              imaginary = 0.02)             # Generando los complejos

int.solve1 <- sapply(z1, integral.compleja1)

G1 <- function(z){
  z-sqrt(z^2-2)
}

G1.solve <- sapply(z1, G1)


# Caso 2 lambda en (-sqrt(2), 0)
z2 <- complex(length.out = 10, real = sample(runif(n = 10, min = -sqrt(2), max = 0)), 
              imaginary = 0.02)   # Generando los complejos

int.solve2 <- sapply(z2, integral.compleja2)

G2 <- function(z){
  z+sqrt(z^2-2)
}



