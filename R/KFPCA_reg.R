#' @title Kendall Functional Principal Component Analysis (KFPCA) for dense and regular design
#'
#' @description KFPCA for non-Gaussian functional data with dense and regular design.
#'
#' @param Lt A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the observation time in ascending order for each subject. The observation times are the same for each subject.
#' @param Ly A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the measurements of each subject at the observation time correspond to \code{Lt}.
#' @param nGrid An integer denoting the number of observation time for each subject.
#' @param nK An integer denoting the number of FPCs.
#' @param fdParobj A functional parameter object for the smoothing of mean function and eigenfunctions. For more detail, see \code{\link[fda]{smooth.basis}}.
#'
#' @return A \code{list} containing the following components:
#' \item{meanfd}{A functional data object for the mean function estimates.}
#' \item{FPC_list}{A \code{list} containing \code{nK} functional data objects, which are the eigenfunction estimates.}
#' \item{score}{A \emph{n} by \code{nK} \code{matrix} containing the estimates of the FPC scores, where \emph{n} is the sample size.}
#' \item{CompTime}{A scalar denoting the computation time.}
#' @export
#'
#' @references \cite{Rou Zhong, Shishi Liu, Haocheng Li, Jingxiao Zhang (2021). "Functional principal component analysis estimator for non-Gaussian data." <arXiv: https://arxiv.org/abs/2102.01286>.}
#'
#' @examples
#' # Generate data
#' n <- 100
#' interval <- c(0, 10)
#' lambda_1 <- 16 #the first eigenvalue
#' lambda_2 <- 9 #the second eigenvalue
#' eigfun <- list()
#' eigfun[[1]] <- function(x){cos(pi * x/10)/sqrt(5)}
#' eigfun[[2]] <- function(x){sin(pi * x/10)/sqrt(5)}
#' score <- cbind(rnorm(n, 0, sqrt(lambda_1)), rnorm(n, 0, sqrt(lambda_2)))
#' DataNew <- GenDataKL(n, interval = interval, sparse = 51, regular = TRUE,
#'                      meanfun = function(x){0}, score = score,
#'                      eigfun = eigfun, sd = sqrt(0.25))
#' basis <- fda::create.bspline.basis(interval, nbasis = 13, norder = 4,
#'                               breaks = seq(0, 10, length.out = 11))
#' #KFPCA
#' Klist <- KFPCA_reg(DataNew$Lt, DataNew$Ly, nGrid = 51, nK = 2, fdParobj = basis)
#' plot(Klist$FPC_list[[1]])
#' plot(Klist$FPC_list[[2]])
#'
KFPCA_reg <- function(Lt, Ly, nGrid, nK, fdParobj){

  startFPCA <- Sys.time() #Starting time for KFPCA

  # n: sample size
  n <- length(Lt)
  gridequal <- Lt[[1]]

  X <- matrix(0, nrow = n, ncol = nGrid)
  for(i in 1:n){
    X[i,] <- Ly[[i]]
  }

  # mean function estimate
  Xfd <- fda::smooth.basis(gridequal, t(X), fdParobj)$fd
  muest_fd <- fda::mean.fd(Xfd)

  # Kendall's tau function
  K <- matrix(0, ncol = nGrid, nrow = nGrid)
  for(i in 2:n){
    for(j in 1:(i-1)){
      XX <- X[i,] - X[j,]
      K <- K + (XX %*% t(XX))/sum(XX^2)
    }
  }
  K <- 2 * K/(n * (n-1))

  # eigenfunction estimates
  eig_K <- eigen(K)
  eigK_fd <- list()
  for(k in 1:nK){
    eigK_fd_k <- fda::smooth.basis(gridequal, eig_K$vectors[,k], fdParobj)$fd
    eigK_fd[[k]] <- eigK_fd_k * (as.numeric(fda::inprod(eigK_fd_k, eigK_fd_k)))^(-1/2)
  }

  # score estimates
  score_K <- matrix(0, nrow = n, ncol = nK)
  for(i in 1:n){
    for(j in 1:nK){
      score_K[i, j] <- fda::inprod(Xfd[i] - muest_fd, eigK_fd[[j]])
    }
  }

  # Computation time
  endFPCA <- Sys.time()
  CompTime <- endFPCA - startFPCA

  # The return object
  ret <- list()
  ret$meanfd <- muest_fd
  ret$FPC_list <- eigK_fd
  ret$score <- score_K
  ret$CompTime <- CompTime

  return(ret)

}
