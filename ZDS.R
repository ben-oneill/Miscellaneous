#' Zero-Defect Sampling (ZDS) analysis for defective items in a finite population
#'
#' \code{ZDS} returns zero-defect sampling analysis for defective items in a finite population.
#'
#' This function implements Bayesian analysis of a Zero-Defect Sampling (ZDS) scheme for defective
#' items in a finite population.  The analysis assumes that the sample is obtained by simple random
#' sampling from the population and defects in the sampled items are detected with fixed probability
#' given by the \code{prob} parameter.  The output of the function is a list of model elements with
#' the class \code{'zds'}, containing the parameters and the log-prior, sampling log-probabilities,
#' log-posterior, and log-posterior-expectation of the number of defective items.  These are given 
#' for all possible sample sizes up to the population size.  The output has a custom printing and
#' plotting method to give user-friendly outputs and visualisation for the analysis.
#' 
#' By default, the analysis uses a discrete uniform prior for the number of defective items in the
#' population.  If the user wishes to use a different prior they can vary the \code{prior.rate} parameter
#' to use a geometric prior with the specified rate or they can input a custom prior as the \code{prior} input.
#' In the latter case the prior should be a vector of probabilities/log-probabilities with length \code{N+1}.
#'
#' @usage \code{ZDS(N, prob = 1, prior = NULL, prior.rate = 0, prior.warning = TRUE)}
#' @param N The population size
#' @param prob The probability of detecting a defect in a sampled item (if it is defective)
#' @param prior An optional vector of prior probabilities for the number of defective items in the population
#' @param prior.rate Rate parameter for the default prior (exponential decaying prior)
#' @param log.prior Logical; if \code{TRUE} the input \code{prior} is treated as a vector of log-probabilities
#' @param prior.warning Logical; if \code{TRUE} the function will warn you if your prior does not sum to one
#' @return A list of class \code{defect} containing prior and posterior information

ZDS <- function(N, prob = 1, prior = NULL, prior.rate = 0, log.prior = FALSE, prior.warning = TRUE) {

  #Check input N
  if (!is.numeric(N))               stop('Error: Population size N should be numeric')
  if (length(N) != 1)               stop('Error: Population size N should be a single value')
  NN <- as.integer(N)
  if (N != NN)                      stop('Error: Population size N must be an integer')
  if (N < 1)                        stop('Error: Population size cannot be less than one')
  if (N == Inf)                     stop('Error: Population must be finite')

  #Check input prob
  if (!is.numeric(prob))            stop('Error: Input prob should be numeric')
  if (length(prob) != 1)            stop('Error: Input prob should be a single value')
  if (prob < 0)                     stop('Error: Input prob cannot be less than zero')
  if (prob > 1)                     stop('Error: Input prob cannot be greater than one')

  #Check input prior.rate
  if (!is.numeric(prior.rate))      stop('Error: Input prior.rate should be numeric')
  if (length(prior.rate) != 1)      stop('Error: Input prior.rate should be a single value')
  if (min(prior.rate) < 0)          stop('Error: Input prior.rate should be a non-negative value')
  
  #Check input log.prior and prior.warning
  if (!is.logical(log.prior))       stop('Error: Input log.prior should be logical')
  if (length(log.prior) != 1)       stop('Error: Input log.prior should be a single logical value')
  if (!is.logical(prior.warning))   stop('Error: Input prior.warning should be logical')
  if (length(prior.warning) != 1)   stop('Error: Input prior.warning should be a single logical value')
  
  #Check input prior
  if (!is.null(prior)) {
    if (!is.numeric(prior))         stop('Error: Input prior should be numeric')
    if (length(prior) != NN+1)      stop('Error: Input prior should have length N+1')
    if (log.prior) {
      LSUM <- matrixStats::logSumExp(prior)
      LOGPRIOR <- prior - LSUM
      if (LSUM != 0) {
        warning('Input prior had probabilities that did not sum to one --- it has been scaled to sum to one') }
    } else {
      SUM  <- sum(prior)
      LOGPRIOR <- log(prior)/SUM
      if (SUM != 1)  {
        warning('Input prior had probabilities that did not sum to one --- it has been scaled to sum to one') } } }
  if (is.null(prior)) {
    if (prior.rate == 0) { 
      LOGPRIOR <- rep(0, NN+1) 
    } else { 
      LOGPRIOR <- dexp(0:NN, rate = prior.rate, log = TRUE) }
    LOGPRIOR <- LOGPRIOR - matrixStats::logSumExp(LOGPRIOR) }
  names(LOGPRIOR) <- sprintf('r[%s]', 0:NN)
  
  #Compute the binomial-hypergeometric log-probabilities
  LOGPROBS <- matrix(-Inf, nrow = NN+1, ncol = NN+1)
  rownames(LOGPROBS) <- sprintf('n[%s]', 0:NN)
  colnames(LOGPROBS) <- sprintf('r[%s]', 0:NN)
  for (n in 0:NN) {
  for (r in 0:NN)  {
    ss <- 0:NN
    T1 <- dbinom(0, size = ss, prob = prob, log = TRUE)
    T2 <- dhyper(ss, r, NN-r, n, log = TRUE)
    LOGPROBS[n+1, r+1] <- matrixStats::logSumExp(T1 + T2) } }

  #Compute the posterior log-probabilities
  LOGPOST <- matrix(-Inf, nrow = NN+1, ncol = NN+1)
  rownames(LOGPOST) <- sprintf('n[%s]', 0:NN)
  colnames(LOGPOST) <- sprintf('r[%s]', 0:NN)
  for (n in 0:NN) {
    LOGPOST[n+1, ] <- LOGPROBS[n+1, ] + LOGPRIOR - matrixStats::logSumExp(LOGPROBS[n+1, ] + LOGPRIOR) }
  
  #Compute the posterior log-expected-values
  LOGPOSTEV <- rep(-Inf, NN+1)
  names(LOGPOSTEV) <- sprintf('D[%s]', 0:NN)
  for (n in 0:NN) {
    LOGPOSTEV[n+1] <- matrixStats::logSumExp(log(0:NN) + LOGPOST[n+1, ]) }
  
  #Generate the output list
  OUT <- list(parameters = list(N = NN, prob = prob), 
              logprior = LOGPRIOR, logprobs = LOGPROBS, logpost = LOGPOST, logpostEV = LOGPOSTEV)
  if (is.null(prior)) {
    OUT$parameters$prior.default <- TRUE
    OUT$parameters$prior.rate    <- prior.rate
  } else {
    OUT$parameters$prior.default <- FALSE }
  class(OUT) <- c('zds', 'list')
  
  #Return the output
  OUT }


print.zds <- function(object) {
  
  #Check object class
  if (!('zds' %in% class(object)))        stop('This print method is for objects of class \'zds\'')
  
  #Extract information
  N         <- object$parameters$N
  prob      <- object$parameters$prob
  DEFAULT   <- object$parameters$prior.default
  LOGPRIOR  <- object$logprior
  LOGPOSTEV <- object$logpostEV
  
  #Print header
  cat('\n    Defect Analysis using Zero-Defect-Sampling \n \n')
  if (prob == 1) { 
    cat(paste0('Analysis is for a population of ', format(N, big.mark = ','), ' items.  ')) 
  } else {
    cat(paste0('Analysis is for a population of ', format(N, big.mark = ','), 
               ' items with sampling defect-detection probability = ', prob, '.  ')) }
  
  #Print prior 
  if (DEFAULT) {  
    cat(paste0('Analysis uses the default prior distribution \n(geometric distribution with decay rate = ', 
               object$parameters$prior.rate, '). \n'))
  } else { 
    STRTEXT <- capture.output(str(unname(exp(LOGPRIOR))))
    cat(paste0('Analysis uses prior distribution supplied by user \n(prior probabilities =', STRTEXT, ' ). \n')) }
  cat('\n')
  
  #Print posterior 
  cat('Posterior expected value of defective items in population: \n \n')
  str(unname(exp(MODEL1$logpostEV)))
  cat('\n') }


plot.zds <- function(object, sample.size = NULL, max.defect.prop = NULL) {
  
  #Check object class
  if (!('zds' %in% class(object)))        stop('This print method is for objects of class \'zds\'')
  
  #Extract information
  N         <- object$parameters$N
  prob      <- object$parameters$prob
  DEFAULT   <- object$parameters$prior.default
  LOGPRIOR  <- object$logprior
  LOGPOSTEV <- object$logpostEV
  
  #Check other inputs
  if ((!is.null(sample.size))&(!is.null(max.defect.prop))) {
    stop('Error: Input sample.size or max.defect.prop (or neither) but not both') }
  if (!is.null(sample.size)) {
    if (!is.numeric(sample.size))         stop('Error: Input sample.size must be numeric')
    nn <- as.integer(sample.size)
    if (any(nn != sample.size))           stop('Error: Input sample.size must be an integer')
    if (min(sample.size) < 0)             stop('Error: Input sample.size must be a non-negative integer')
    if (max(sample.size) > N)             stop('Error: Input sample.size cannot be greater than population size') }
  if (!is.null(max.defect.prop)) {
    if (!is.numeric(max.defect.prop))     stop('Error: Input max.defect.prop must be numeric')
    ee <- max.defect.prop
    if (min(max.defect.prop) < 0)         stop('Error: Input max.defect.prop must be a proportion')
    if (max(max.defect.prop) > 1)         stop('Error: Input max.defect.prop must be a proportion') }

  #Generate theme and plotdata
  THEME    <- ggplot2::theme(
                plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold'),
                plot.subtitle = ggplot2::element_text(hjust = 0.5),
                axis.title.y  = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)))
  PLOTDATA <- data.frame(Value = 0:N, postev = c(exp(LOGPOSTEV[1:N]), NA), prior = c(exp(LOGPRIOR)))
  PRIORMAX <- exp(max(LOGPRIOR))
  
  #Generate subtitles
  if (DEFAULT) { 
    RATE <- object$parameters$prior.rate
    if (RATE == 0) { SUB1 <- '(Uniform distribution)' }
    if (RATE >  0) { SUB1 <- paste0('(Geometric distribution with decay rate = ', RATE, ')') }
  } else { 
    SUB1 <- paste0('(Distribution supplied by user)') }
  SUB2 <- '(Zero-defect sampling --- No observed defects in sample)'
  
  #Generate figures
  FIGURE1  <- ggplot2::ggplot(ggplot2::aes(x = Value, y = prior), data = PLOTDATA) +
              ggplot2::geom_line(colour = 'red') + 
              ggplot2::scale_y_continuous(limits = c(0, min(1, 1.2*PRIORMAX))) + 
              THEME + ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)) +
              ggplot2::ggtitle('Prior Distribution') +
              ggplot2::labs(subtitle = SUB1) +
              ggplot2::xlab('Defective Items') + ggplot2::ylab('Prior Probability')
  FIGURE2  <- ggplot2::ggplot(ggplot2::aes(x = Value, y = postev), data = PLOTDATA) +
              ggplot2::geom_line(colour = 'blue') + 
              THEME + ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)) +
              ggplot2::scale_y_log10(breaks = scales::trans_breaks('log10', function(x) 10^x),
                                     labels = scales::trans_format('log10', scales::math_format(10^.x))) +
              ggplot2::ggtitle('Posterior Inference') +
              ggplot2::labs(subtitle = SUB2) +
              ggplot2::xlab('Sample Size') + ggplot2::ylab('Expected Defective Items')
  
  #Add identified sample sizes
  if (!is.null(sample.size)) {
    for (i in 1:length(nn)) {
      FIGURE2 <- FIGURE2 + 
                 ggplot2::geom_segment(x    = nn[i], y    = LOGPOSTEV[nn[i]+1]/log(10), 
                                       xend = 0,     yend = LOGPOSTEV[nn[i]+1]/log(10), 
                                       colour = 'darkgrey') + 
                 ggplot2::geom_segment(x    = nn[i], y    = LOGPOSTEV[nn[i]+1]/log(10), 
                                       xend = nn[i], yend = -Inf, 
                                       colour = 'darkgrey') +
                 ggplot2::geom_point(ggplot2::aes(x = nn[i], y = exp(LOGPOSTEV[nn[i]+1])), 
                                     colour = 'blue', size = 2) } }
  if (!is.null(max.defect.prop)) {
    for (i in 1:length(ee)) {
      nn <- sum(LOGPOSTEV > log(N) + log(ee[i]))
      if (nn < N) {
      FIGURE2 <- FIGURE2 + 
                 ggplot2::geom_segment(x    = nn, y    = LOGPOSTEV[nn+1]/log(10), 
                                       xend = 0,  yend = LOGPOSTEV[nn+1]/log(10), 
                                       colour = 'darkgrey') + 
                 ggplot2::geom_segment(x    = nn, y    = LOGPOSTEV[nn+1]/log(10), 
                                       xend = nn, yend = -Inf, 
                                       colour = 'darkgrey') +
                 ggplot2::geom_point(ggplot2::aes(x = nn, y = exp(LOGPOSTEV[nn+1])), 
                                     colour = 'blue', size = 2) } } }
  
  #Combine into single plot
  PLOT     <- suppressWarnings(cowplot::plot_grid(FIGURE1, FIGURE2, ncol = 1, align = 'v'))
  
  #Print the plot
  PLOT }

