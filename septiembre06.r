# EJERCICIO 1

# PARAMETROS
m <- 10000
p1 <- 10
p2 <- 100

# Semic�culo Wigner
x <- seq(from = -sqrt(2), to = sqrt(2), by = 0.1)
Y <- function(x){
  y <- sqrt(2-x^2)
  y <- y/pi
  y
}

plot(x, Y(x), type = 'l', main = 'Semicírculo Wigner', col='blue')

# GOE
# Esta funci�n genera las lambdas escaladas
Lambda.GOE <- function(p1){
  H <- replicate(n=p1,{sample(rnorm(n=p1, 0, 1))})
  Hs <- (H+t(H))/2 # GOE
  Lambda.Hs <- eigen(Hs)
  Lambda.Hs <- Lambda.Hs$values
  Lambda.Hs <- (1/sqrt(p1))*Lambda.Hs
  Lambda.Hs
}

# P=10
x.GOE <- replicate(n = m, Lambda.GOE(10))
hist(x.GOE, freq = FALSE)

#P=100
x.GOE <- replicate(n = m, Lambda.GOE(100))
hist(x.GOE, freq = FALSE)


# GUE
# Esta función genera las lambdas escaladas
Lambda.GUE <- function(p1){
  H2 <- replicate(n=p1, complex(length.out = p1, real = sample(rnorm(n = p1, mean = 0, sd = 1)), imaginary = sample(rnorm(n = p1, mean = 0, sd = 1))))
  H2.conjug <- Conj(H2)
  H.her <- (H2+t(H2.conjug))/2
  Lambda.H.her <- eigen(H.her)
  Lambda.H.her <- Lambda.H.her$values
  Lambda.H.her <- (1/sqrt(20))*Lambda.H.her
  Lambda.H.her
}

# P=10
x.GUE <- replicate(n = m, Lambda.GUE(10))
hist(x.GUE, freq = FALSE)

#P=100
x.GUE <- replicate(n = m, Lambda.GUE(100))
hist(x.GUE, freq = FALSE)

# GSE
# Esta función genera las lambdas escalonadas

Lambda.GSE <- function(p1){
  X <- replicate(n=p1, complex(length.out = p1, real = sample(rnorm(n = p1, mean = 0, sd = 1)), imaginary = sample(rnorm(n = p1, mean = 0, sd = 1))))
  Y <- replicate(n=p1, complex(length.out = p1, real = sample(rnorm(n = p1, mean = 0, sd = 1)), imaginary = sample(rnorm(n = p1, mean = 0, sd = 1))))
  X.conjug <- Conj(X)
  Y.conjug <- Conj(Y)
  a <- cbind(X, Y)
  b <- cbind(-Y.conjug, X.conjug)
  H3 <- rbind(a,b)
  H.sc <- (H3+t(Conj(H3)))/2
  Lambda.H.sc <- eigen(H.sc)
  Lambda.H.sc <- Lambda.H.sc$values
  Lambda.H.sc <- (1/sqrt(4*p1))*Lambda.H.sc
  Lambda.H.sc
}

# P=10
x.GSE <- replicate(n = m, Lambda.GSE(10))
hist(x.GSE, freq = FALSE)
lines(density(x.GSE ))

#P=100
x.GSE <- replicate(n = m, Lambda.GSE(100))
hist(x.GSE, freq = FALSE, density = TRUE)

# EJERCICIO 2

# Función que calcula las diferecias

diferencial.lambda <- function(p, p1=100){
  Li <- Lambda.GOE(p1) # eigenvalores de matriz simétrica de pxp 
  dif.LI.Lj <- diff(Li)
  dif.LI.Lj <- abs(dif.LI.Lj)
}

m <- 1000
GOE.100 <- replicate(n=m, diferencial.lambda(100))
hist(GOE.100, freq = FALSE)
lines(density(GOE.100))


GOE.2 <- replicate(n=m, diferencial.lambda(2))
hist(GOE.2, freq = FALSE)
lines(density(GOE.2))



