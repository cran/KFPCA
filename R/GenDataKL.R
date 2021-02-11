#' @title Generate functional/longitudinal data via KL expansion
#'
#' @description Generate functional/longitudinal data via Karhunen–Loève expansion.
#'
#' @param n number of sample size.
#' @param interval A \code{vector} of length two denoting the supporting interval.
#' @param sparse A \code{vector} denoting the possible numbers of observation size. The elements are chosen with equal chance.
#' @param meanfun A function for the mean.
#' @param score A \emph{n} by \code{nK} \code{matrix} containing the estimates of the FPC scores, where \code{nK} is the number of FPCs.
#' @param eigfun A \code{list} containing the eigenfunctions.
#' @param sd A scalar denoting the standard deviation of measurement errors.
#'
#' @return A \code{list} containing the following components:
#' \item{Lt}{A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the observation time in ascending order for each subject.}
#' \item{Ly}{A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the measurements of each subject at the observation time correspond to \code{Lt}.}
#' @export
#'
#' @examples
#' n <- 100
#' interval <- c(0, 10)
#' lambda_1 <- 9 #the first eigenvalue
#' lambda_2 <- 1.5 #the second eigenvalue
#' eigfun <- list()
#' eigfun[[1]] <- function(x){cos(pi * x/10)/sqrt(5)}
#' eigfun[[2]] <- function(x){sin(pi * x/10)/sqrt(5)}
#' score <- cbind(rnorm(n, 0, sqrt(lambda_1)), rnorm(n, 0, sqrt(lambda_2)))
#' DataNew <- GenDataKL(n, interval = interval, sparse = 6:8,
#'                      meanfun = function(x){0}, score = score,
#'                      eigfun = eigfun, sd = sqrt(0.1))
GenDataKL <- function(n, interval, sparse, meanfun, score, eigfun, sd){

  # nK: the number of FPCs
  nK <- ncol(score)

  Lt <- list()
  Ly <- list()
  for(i in 1:n){
    num <- sample(sparse, 1) #observation size for the i-th subject
    Lt[[i]] <- sort(unique(stats::runif(num, min = interval[1], max = interval[2])))
    scorei <- score[i,]
    y <- meanfun(Lt[[i]])
    for(j in 1:nK){
      y <- y + scorei[j] * eigfun[[j]](Lt[[i]])
    }
    Ly[[i]] <- y + stats::rnorm(length(Lt[[i]]), 0, sd = sd)
  }

  ret <- list()
  ret$Lt <- Lt
  ret$Ly <- Ly

  return(ret)
}
