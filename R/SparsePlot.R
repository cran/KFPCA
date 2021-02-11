#' @title Sparse plot
#'
#' @description Create sparse plot to see the sparsity of the data.
#'
#' @param Lt A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the observation time in ascending order for each subject.
#' @param interval A \code{vector} of length two denoting the supporting interval.
#' @param ... Other arguments passed into \code{\link[graphics]{plot}}.
#'
#' @return Create the corresponding sparse plot.
#'
#' @details For the sparse plot, x-axis is the observation time while y-axis represents various subjects.
#' @export
#'
#' @examples
#' # Generate data
#' n <- 100
#' interval <- c(0, 10)
#' lambda_1 <- 9 #the first eigenvalue
#' lambda_2 <- 1.5 #the second eigenvalue
#' eigfun <- list()
#' eigfun[[1]] <- function(x){cos(pi * x/10)/sqrt(5)}
#' eigfun[[2]] <- function(x){sin(pi * x/10)/sqrt(5)}
#' score <- cbind(rnorm(n, 0, sqrt(lambda_1)), rnorm(n, 0, sqrt(lambda_2)))
#' # DataNew1 and DataNew2 have different sparsity
#' DataNew1 <- GenDataKL(n, interval = interval, sparse = 6:8,
#'                       meanfun = function(x){0}, score = score,
#'                       eigfun = eigfun, sd = sqrt(0.1))
#' DataNew2 <- GenDataKL(n, interval = interval, sparse = 2:4,
#'                       meanfun = function(x){0}, score = score,
#'                       eigfun = eigfun, sd = sqrt(0.1))
#' # Create sparse plots
#' par(mfrow = c(1, 2))
#' SparsePlot(DataNew1$Lt, interval = interval)
#' SparsePlot(DataNew2$Lt, interval = interval)
#' par(mfrow = c(1, 1))
SparsePlot <- function(Lt, interval, ...){
  # n: sample size
  n <- length(Lt)

  # para: parameters for plot()
  para <- list(main = "Sparse Plot", xlab = "time", ylab = "ID", tck = -0.025,
               mgp = c(1.5, 0.5, 0))
  parauser <- list(...)
  para[names(parauser)] <- parauser

  pin_init <- graphics::par(no.readonly = T)$pin
  on.exit(graphics::par(pin = pin_init))

  if(is.null(para$pin)){
    graphics::par(pin = c(3,3))
  }else{
    graphics::par(pin = para$pin)
  }
  do.call(graphics::plot, c(para, list(x = interval, y = c(0, n), type = "n")))
  for(i in 1:n){
    graphics::points(Lt[[i]], rep(i, length(Lt[[i]])))
  }
}
