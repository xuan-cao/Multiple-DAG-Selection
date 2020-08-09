
rm(list=ls())
# install.packages("doMC")
# install.packages("lars")

# # to use data set from Ryan's code
# X = cbind(X0, as.matrix(y))
# j = ncol(X)

# read source files
setwd("~/OneDrive - University of Cincinnati/multiple_DAG/sim/")
source("functions.R")

########################################################################
# EXAMPLE
########################################################################
p = 150
n = 100
K = 3

# generation of precision matrix (based on MCD) for k=1
set.seed(123)
A01.vec = rep(0, p*(p-1)/2)
A01.len = length(A01.vec)
s01.rowwise = rep(0, p-1)
s01.tot = floor(A01.len*0.02) # number of nonzero signals on A01 sparsity setting!!!
S01.tot = sample(1:A01.len, size = s01.tot, replace = FALSE) # locations of nonzero signals on A01
sign.tot = sample(c(1,-1), size = s01.tot, replace = TRUE) # signs of nonzero signals on A01
A01.vec[ S01.tot ] = runif(n = s01.tot, min = 0.3, max = 0.7)*sign.tot
A01 = matrix(0, p,p)
for(i in 2:p){
  for(j in 1:(i-1)){
    A01[i,j] = A01.vec[(i-1)*(i-2)/2 + j]
  }
  s01.rowwise[i] = length(which(A01[i,] != 0))
}
D01 = diag(runif(n = p, min = 2, max = 5))
Omega01 = t(diag(p) - A01)%*%diag(1/diag(D01))%*%(diag(p) - A01)
eig01 = eigen(Omega01, symmetric = TRUE, only.values = TRUE)$val

# generation of precision matrix (based on MCD) for k=2, remove 5 edges and add 5 edges from DAG1 under 2%, 5 edges under 4% 
A02.vec = rep(0, p*(p-1)/2)
s02.rowwise = rep(0, p-1)
A02.vec[S01.tot] = 1
add1 = sample(which(A02.vec == 0), 5, replace = FALSE)
remove1 = sample(which(A02.vec == 1), 5, replace = FALSE)
A02.vec[add1] = 1
A02.vec[remove1] = 0
S02.tot = which(A02.vec != 0)
A02.vec[ S02.tot ] = runif(n = s01.tot, min = 0.3, max = 0.7)*sign.tot
A02 = matrix(0, p,p)
for(i in 2:p){
  for(j in 1:(i-1)){
    A02[i,j] = A02.vec[(i-1)*(i-2)/2 + j]
  }
  s02.rowwise[i] = length(which(A02[i,] != 0))
}
D02 = diag(runif(n = p, min = 2, max = 5))
Omega02 = t(diag(p) - A02)%*%diag(1/diag(D02))%*%(diag(p) - A02)
eig02 = eigen(Omega02, symmetric = TRUE, only.values = TRUE)$val


# generation of precision matrix (based on MCD) for k=3, remove 20 edges and add 20 edges from DAG2 under 2%, 20 edges under 4%
A03.vec = rep(0, p*(p-1)/2)
s03.rowwise = rep(0, p-1)
A03.vec[S02.tot] = 1
add2 = sample(which(A03.vec == 0), 5, replace = FALSE)
remove2 = sample(which(A03.vec == 1), 5, replace = FALSE)
A03.vec[add2] = 1
A03.vec[remove2] = 0
S03.tot = which(A03.vec != 0)
A03.vec[ S03.tot ] = runif(n = s01.tot, min = 0.3, max = 0.7)*sign.tot
A03 = matrix(0, p,p)
for(i in 2:p){
  for(j in 1:(i-1)){
    A03[i,j] = A03.vec[(i-1)*(i-2)/2 + j]
  }
  s03.rowwise[i] = length(which(A03[i,] != 0))
}
D03 = diag(runif(n = p, min = 2, max = 5))
Omega03 = t(diag(p) - A03)%*%diag(1/diag(D03))%*%(diag(p) - A03)
eig03 = eigen(Omega03, symmetric = TRUE, only.values = TRUE)$val

# hyperparameters
alpha = 0.999
gamma = 1
nu0 = 0.1
c1 = 1 #c1 in ESC paper, doesn't matter when set to 1
c2 = 2
b = 1/(p*(K - 1)) #c2 in the paper, when comparing, make sure under the same data!!
#b = 10
niter = 1000
nburn = 0.2*niter
#nburn = 100
nadap = 0


########################################################################
# single simulation study
########################################################################
set.seed(123)
X1 = matrix(0, nrow=n, ncol=p)
X1[,1] = rnorm(n, sd = sqrt(D01[1,1]))
for(j in 2:p){
  mean.vec.j = X1[, 1:(j-1)]%*%as.matrix(A01[j, 1:(j-1)])
  X1[,j] = rnorm(n, mean = mean.vec.j, sd = sqrt(D01[j,j]))
}

X2 = matrix(0, nrow=n, ncol=p)
X2[,1] = rnorm(n, sd = sqrt(D02[1,1]))
for(j in 2:p){
  mean.vec.j = X2[, 1:(j-1)]%*%as.matrix(A02[j, 1:(j-1)])
  X2[,j] = rnorm(n, mean = mean.vec.j, sd = sqrt(D02[j,j]))
}


X3 = matrix(0, nrow=n, ncol=p)
X3[,1] = rnorm(n, sd = sqrt(D03[1,1]))
for(j in 2:p){
  mean.vec.j = X3[, 1:(j-1)]%*%as.matrix(A03[j, 1:(j-1)])
  X3[,j] = rnorm(n, mean = mean.vec.j, sd = sqrt(D03[j,j]))
}





# lambdas <- c(1,2*qnorm(p=(200/(p*(1:(p-1)))),lower.tail=FALSE))/sqrt(n) #for n 100, p 300
# LassoAdjseq <- DAGLassoseq(Y=scale(X1),lambda.seq=lambdas)
# LassoAdjseq[which(LassoAdjseq!=0)] <- 1
# evaluation.dag((A01 != 0), LassoAdjseq)

X = array(0, dim = c(K, dim(X1)))
X[1,,] = X1
X[2,,] = X2
X[3,,] = X3

########################################################################
# joint DAG posterior inference with ESC prior
########################################################################
time = Sys.time()
library(parallel)
source("functions.R")
res = mdag(X, alpha, gamma, nu0, c1, c2, c3=NULL, b, niter, nburn)
elapsed = Sys.time() - time
elapsed

#sort return list
l1 = list()
for (k in 1:K) {
  l2 = list(NULL)
  for (j in 1:(p-1)) {
    m <- matrix(0, niter, j)
    for (t in (nburn + 1):(nburn + niter)) {
      m[t-nburn, ] = res[[j]][[t]][k, ]
    }
    l2[[j]] = m
  }
  l1[[k]] = l2
}



########################################################################
# results
########################################################################
#DAG1 evaluation
incl.pr1 <- matrix(0,p,p)

for(j in 2:p){
  Sj.mat = l1[[1]][[j-1]]
  incl.pr1[j,seq(1:(j-1))] <- apply(Sj.mat, 2, mean)}
print(evaluation.dag((A01[lower.tri(A01)] != 0),1*(incl.pr1[lower.tri(incl.pr1)]>0.5)))
print(Joint_auc(as.vector(1*(A01[lower.tri(A01)] != 0)), as.vector(1*(incl.pr1[lower.tri(incl.pr1)]))))

#DAG2 evaluation
incl.pr2 <- matrix(0,p,p)

for(j in 2:p){
  Sj.mat = l1[[2]][[j-1]]
  incl.pr2[j,seq(1:(j-1))] <- apply(Sj.mat, 2, mean)}
print(evaluation.dag((A02[lower.tri(A02)] != 0),1*(incl.pr2[lower.tri(incl.pr2)]>0.5)))
print(Joint_auc(as.vector(1*(A02[lower.tri(A02)] != 0)), as.vector(1*(incl.pr2[lower.tri(incl.pr2)]))))

#DAG3 evaluation
incl.pr3 <- matrix(0,p,p)

for(j in 2:p){
  Sj.mat = l1[[3]][[j-1]]
  incl.pr3[j,seq(1:(j-1))] <- apply(Sj.mat, 2, mean)}
print(evaluation.dag((A03[lower.tri(A03)] != 0),1*(incl.pr3[lower.tri(incl.pr3)]>0.5)))
print(Joint_auc(as.vector(1*(A03[lower.tri(A03)] != 0)), as.vector(1*(incl.pr3[lower.tri(incl.pr3)]))))


true.vec1 = 1*(as.vector(A01[lower.tri(A01)] != 0))
true.vec2 = 1*(as.vector(A02[lower.tri(A02)] != 0))
true.vec3 = 1*(as.vector(A03[lower.tri(A03)] != 0))

est.vec11 = 1*(as.vector(incl.pr1[lower.tri(incl.pr1)]>0.5))
est.vec21 = 1*(as.vector(incl.pr2[lower.tri(incl.pr2)]>0.5))
est.vec31 = 1*(as.vector(incl.pr3[lower.tri(incl.pr3)]>0.5))

est.vec1 = 1*(as.vector(incl.pr1[lower.tri(incl.pr1)]>0.5))
est.vec2 = 1*(as.vector(incl.pr2[lower.tri(incl.pr2)]>0.5))
est.vec3 = 1*(as.vector(incl.pr3[lower.tri(incl.pr3)]>0.5))


print(evaluation.dag(c(true.vec1, true.vec2, true.vec3), c(est.vec1, est.vec2, est.vec3)))
print(Joint_auc(c(true.vec1, true.vec2, true.vec3), c(est.vec11, est.vec21, est.vec31)))


print(evaluation.dag(c(true.vec2[which(true.vec1 - true.vec2 == -1)], true.vec2[which(true.vec1 - true.vec2 == 1)], true.vec3[which(true.vec2 - true.vec3 == -1)], true.vec3[which(true.vec2 - true.vec3 == 1)]), c(est.vec2[which(true.vec1 - true.vec2 == -1)], est.vec2[which(true.vec1 - true.vec2 == 1)], est.vec3[which(true.vec2 - true.vec3 == -1)], est.vec3[which(true.vec2 - true.vec3 == 1)])))
print(Joint_auc(c(true.vec2[which(true.vec1 - true.vec2 == -1)], true.vec2[which(true.vec1 - true.vec2 == 1)], true.vec3[which(true.vec2 - true.vec3 == -1)], true.vec3[which(true.vec2 - true.vec3 == 1)]), c(est.vec21[which(true.vec1 - true.vec2 == -1)], est.vec21[which(true.vec1 - true.vec2 == 1)], est.vec31[which(true.vec2 - true.vec3 == -1)], est.vec31[which(true.vec2 - true.vec3 == 1)])))







