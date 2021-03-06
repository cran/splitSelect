#' 
#' @import foreach
#' 
#' @title Best Split Selection Modeling for Low-Dimensional Data
#'
#' @description \code{splitSelect} performs the best split selection algorithm.
#' 
#' @param x Design matrix.
#' @param y Response vector.
#' @param intercept Boolean variable to determine if there is intercept (default is TRUE) or not.
#' @param G Number of groups into which the variables are split. Can have more than one value.
#' @param group.model Model used for the groups. Must be one of "glmnet" or "LS".
#' @param family Description of the error distribution and link function to be used for the model. Must be one of "gaussian" or "binomial".
#' @param lambdas The shinkrage parameters for the "glmnet" regularization. If NULL (default), optimal values are chosen.
#' @param alphas Elastic net mixing parameter. Should be between 0 (default) and 1.
#' @param nsample Number of sample splits for each value of G. If NULL, then all splits will be considered (unless there is overflow).
#' @param use.all Boolean variable to determine if all variables must be used (default is TRUE).
#' @param fix.partition Optional list with G elements indicating the partitions (in each row) to be considered for the splits.
#' @param fix.split Optional matrix with p columns indicating the groups (in each row) to be considered for the splits.
#' @param parallel Boolean variable to determine if parallelization of the function. Default is FALSE.
#' @param cores Number of cores for the parallelization for the function.
#' @param verbose Boolean variable to determine if console output for cross-validation progress is printed (default is TRUE).
#' 
#' @return An object of class splitSelect.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @seealso \code{\link{coef.splitSelect}}, \code{\link{predict.splitSelect}}
#' 
#' @examples
#' # Setting the parameters
#' p <- 4
#' n <- 30
#' n.test <- 5000
#' beta <- rep(5,4)
#' rho <- 0.1
#' r <- 0.9
#' SNR <- 3
#' # Creating the target matrix with "kernel" set to rho
#' target_cor <- function(r, p){
#'   Gamma <- diag(p)
#'   for(i in 1:(p-1)){
#'     for(j in (i+1):p){
#'       Gamma[i,j] <- Gamma[j,i] <- r^(abs(i-j))
#'     }
#'   }
#'   return(Gamma)
#' }
#' # AR Correlation Structure
#' Sigma.r <- target_cor(r, p)
#' Sigma.rho <- target_cor(rho, p)
#' sigma.epsilon <- as.numeric(sqrt((t(beta) %*% Sigma.rho %*% beta)/SNR))
#' # Simulate some data
#' x.train <- mvnfast::rmvn(30, mu=rep(0,p), sigma=Sigma.r)
#' y.train <- 1 + x.train %*% beta + rnorm(n=n, mean=0, sd=sigma.epsilon)
#' 
#' # Generating the coefficients for a fixed partition of the variables
#' \donttest{
#' split.out <- splitSelect(x.train, y.train, G=2, use.all=TRUE,
#'                          fix.partition=list(matrix(c(2,2), 
#'                                              ncol=2, byrow=TRUE)), 
#'                          fix.split=NULL,
#'                          intercept=TRUE, group.model="glmnet", alphas=0)
#' }
#'
splitSelect <- function(x, y, intercept = TRUE,
                        G, use.all = TRUE,
                        family=c("gaussian", "binomial")[1],
                        group.model=c("glmnet", "LS", "Logistic")[1], lambdas=NULL, alphas = 0,
                        nsample = NULL, fix.partition = NULL, fix.split = NULL,
                        parallel=FALSE, cores=getOption('mc.cores', 2L),
                        verbose=TRUE){

  # Encure numeric entry for response
  y <- as.numeric(y)
  # Check input data
  if (all(!inherits(x, "matrix"), !inherits(x, "data.frame"))) {
    stop("x should belong to one of the following classes: matrix, data.frame.")
  } else if (all(!inherits(y, "matrix"), all(!inherits(y, "numeric")))) {
    stop("y should belong to one of the following classes: matrix, numeric.")
  } else if (any(anyNA(x), any(is.nan(x)), any(is.infinite(x)))) {
    stop("x should not have missing, infinite or nan values.")
  } else if (any(anyNA(y), any(is.nan(y)), any(is.infinite(y)))) {
    stop("y should not have missing, infinite or nan values.")
  } else {
    if(inherits(y, "matrix")) {
      if (ncol(y)>1){
        stop("y should be a vector")
      }
      y <- as.numeric(y)
    }
    len_y <- length(y)
    if (len_y != nrow(x)) {
      stop("y and x should have the same number of rows.")
    }
  }
  if(!(family %in% c("gaussian", "binomial"))){
    stop("group.model should be one of \"gaussian\" or \"binomial\".")
  }
  if(!is.null(alphas)){
    if (!inherits(alphas, "numeric")) {
      stop("alphas should be numeric")
    } else if (any(any(alphas < 0), any(alphas > 1))) {
      stop("alphas should be a numeric value between 0 and 1.")
    }
  }
  if(!is.null(lambdas)){
    if (!inherits(lambdas, "numeric")) {
      stop("lambdas should be numeric")
    } else if (any(lambdas < 0)) {
      stop("lambdas should be a numeric non-negative vector of length G.")
    }
  }
  p <- ncol(x) # Storing the number of variables
  if(any(c(!is.numeric(p), length(p)!=1, p<1, !(p%%1==0))))
    stop("p should be a positive interger.")
  if(any(c(!is.numeric(G), any(G<1), !any(G%%1==0), any(G>p))))
    stop("G should be a positive interger less or equal to p.")
  if(!is.null(nsample)){
    if(any(c(!is.numeric(nsample), length(nsample)!=1, nsample<1, !(nsample%%1==0))))
      stop("nsample should be a positive interger.")
  }
  if(!is.null(fix.partition))
    if(any(c(!is.list(fix.partition), length(fix.partition)!=length(G))))
      stop("fix.partition should be a list of size the number of elements of G.")
  if(!is.null(fix.split)){
    if(any(!is.matrix(fix.split), ncol(fix.split)!=p))
      stop("fix.split should be a matrix with p columns.")
    if(use.all)
      if(!is.null(fix.split) && anyNA(fix.split))
        stop("fix.split contains NA values. Set use.all to FALSE if needed.")
  }
  if(!(group.model %in% c("glmnet", "LS", "Logistic"))){
    stop("group.model should be one of \"glmnet\" or \"LS\" or \"Logistic\".")
  }

  # Determining if random splits must be generated (total splits may be too large)
  if(is.null(nsample)){
    
    total.splits <- numeric(1)
    for(G.ind in 1:length(G))
      if(is.null(fix.partition))
        total.splits <- total.splits + nsplit(p=p, G=G[G.ind], use.all=use.all, fix.partition=fix.partition) else
          total.splits <- total.splits + nsplit(p=p, G=G[G.ind], use.all=use.all, fix.partition=fix.partition[[G.ind]])

    if(total.splits > 1e5){
      if(verbose)
        cat("There are over", 1e5, "possible splits. Random splits will be considered instead.\n")
      splits.candidates <- "Sample"
      n <- 1e3
    }
  }
  
  # Generating the splits
  final.splits <- matrix(nrow=0, ncol=p)
  if(!is.null(fix.split))
    final.splits <- fix.split else{
      for(G.ind in 1:length(G)){
        if(is.null(fix.partition)){
          if(is.null(nsample))
            final.splits <- rbind(final.splits, generate_splits(p=p, G=G[G.ind], use.all=use.all, fix.partition=fix.partition, verbose=verbose)) else
              final.splits <- rbind(final.splits, rsplit(n=nsample, p=p, G=G[G.ind], use.all=use.all, fix.partition=fix.partition, verbose=verbose))
        } else{
          if(is.null(nsample))
            final.splits <- rbind(final.splits, generate_splits(p=p, G=G[G.ind], use.all=use.all, fix.partition=fix.partition[[G.ind]], verbose=verbose)) else
              final.splits <- rbind(final.splits, rsplit(n=nsample, p=p, G=G[G.ind], use.all=use.all, fix.partition=fix.partition[[G.ind]], verbose=verbose))
        }
      }
    }
  # Setting NA values as 0
  final.splits[is.na(final.splits)] <- 0
  
  # Consider random splits if there are too many
  if(nrow(final.splits) > 5e3){
    warning("There are over 5,000 splits. A random sample of 5,000 splits will be taken.")
    final.splits <- final.splits[sample(1:nrow(final.splits), 5000, replace=FALSE),]
  }
  
  if(parallel){

    # Registering cores if not already declared
    if (!foreach::getDoParRegistered()){
      cl <- parallel::makePSOCKcluster(cores)
      doParallel::registerDoParallel(cl)
    }

    # Determining the splits to store for each core
    parallel.id <- split(1:nrow(final.splits), factor(sort(rank(1:nrow(final.splits))%%cores)))

    # Hack initialization
    subset.ind <- NULL
    # Parallel computation for the subsets
    splits.betas <- foreach::foreach(subset.ind=1:min(cores, length(parallel.id)), .packages=c("splitSelect"),
                                     .export=c("splitSelect", "splitSelect_coef")) %dopar% {

      # Data to store the mspes for the core
      core.splits <- parallel.id[[subset.ind]]
      if(intercept)
        core.splits.betas <- matrix(nrow=p+1, ncol=length(core.splits)) else
          core.splits.betas <- matrix(nrow=p, ncol=length(core.splits))

      # Generating the adaptive SPLIT coefficients
      for(split.id in 1:length(core.splits)){

        # Store the current split
        current.split <- final.splits[core.splits[split.id],]
        # Generate the adaptive SPLIT coefficients for current
        core.splits.betas[, split.id] <- splitSelect_coef(x=x, y=y, variables.split=as.vector(current.split),
                                                          intercept=intercept,
                                                          family=family,
                                                          group.model=group.model, 
                                                          lambdas=lambdas, alphas=alphas)
      }
      # Returning the splits.betas
      return(core.splits.betas)
    }
    # Unlisting the betas
    splits.betas <- do.call("cbind", splits.betas)
  } else{
    
    # Matrix to store the coefficients
    if(intercept)
      splits.betas <- matrix(nrow=p+1, ncol=nrow(final.splits)) else
        splits.betas <- matrix(nrow=p, ncol=nrow(final.splits))
    
    # Generating the adaptive SPLIT coefficients
    for(split.id in 1:nrow(final.splits)){
      
      # Store the current split
      current.split <- final.splits[split.id,,drop=FALSE]
      # Generate the adaptive SPLIT coefficients for current 
      splits.betas[, split.id] <- splitSelect_coef(x=x, y=y, variables.split=as.vector(current.split), 
                                                   intercept=intercept, 
                                                   family=family,
                                                   group.model=group.model, 
                                                   lambdas=lambdas, alphas=alphas)
    }
  }
  
  # Construct the output
  splitSelect.out <- list(splits=final.splits, betas=splits.betas)
  fn.call <- match.call()
  splitSelect.out <- construct.splitSelect(object=splitSelect.out, 
                                           fn_call=fn.call, 
                                           x=x, y=y, intercept=intercept, 
                                           family=family)
  
  # Returning the splits and the coefficients
  return(splitSelect.out)
}








