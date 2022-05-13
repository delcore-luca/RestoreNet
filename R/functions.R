# get.polygon <- function(r, N){
#   coords <- cbind(r * cos(2*pi*(1:N)/N),
#                   r * sin(2*pi*(1:N)/N))
#   return(coords)
# }

#' Rescaling a clonal tracking dataset
#'
#' Rescales a clonal tracking dataset based on the sequencing depth.
#' @param Y A 3-dimensional array whose dimensions are the time, the cell type and the clone respectively.
#' @details This function rescales a clonal tracking dataset Y according to the formula
#' \deqn{Y_{ijk} \leftarrow Y_{ijk} \cdot \frac{min_{ij}\sum_cY_{ijc}}{\sum_cY_{ijc}}}{Y_{ijk} <- Y_{ijk}*(min_{ij}sum_cY_{ijc})/(sum_cY_{ijc})}
#' @return A rescaled clonal tracking dataset.
#' @export
get.rescaled <- function(Y){
  f <- min(setdiff(c(apply(Y, c(1,2), sum)), 0))/apply(Y, c(1,2), sum)
  f[is.infinite(f)] <- 0
  Y_res <- sapply(1:dim(Y)[3], function(cl){round(Y[,,cl]*f)}, simplify = "array")
  dimnames(Y_res) <- dimnames(Y)
  return(Y_res)
}

#' Clonal pie-chart
#'
#' Draw a clonal pie-chart of a random-effects reaction network.
#' @param re.res output list returned by fit.re().
#' @param txt logical (defaults to FALSE). If TRUE, barcode names will be printed on the pies.
#' @param legend logical (defaults to FALSE). If TRUE, the legend of the pie-chart will be printed.
#' @details This function generates a clonal pie-chart given a previously fitted random-effects model.
#' In this representation each clone \eqn{k}{k} is identified with a pie whose slices are
#' lineage-specific and weighted with \eqn{w_k}{w_k}, defined as the difference between the
#' conditional expectations of the random-effects on duplication and death parameters, that is
#' \deqn{w_k = E_{u\vert \Delta Y; \hat{\psi}}[u^k_{\alpha_{lin}}] - E_{u\vert \Delta Y; \hat{\psi}}[u^k_{\delta_{lin}}]}{Eu|y;ψ[α_lin] - Eu|y;ψ[β_lin]},
#' where \eqn{\texttt{lin}}{lin} is a cell lineage.
#' The diameter of the \eqn{k}{k}-th pie is proportional to the euclidean 2-norm of \eqn{w_k}{w_k}.
#' Therefore, the larger the diameter, the more the corresponding clone is expanding
#' into the lineage associated to the largest slice.
#' @examples
#' \donttest{
#' library(RestoreNet)
#' library(ggplot2)
#' library(scatterpie)
#' rcts <- c("A->1", "B->1", "C->1", "D->1",
#'           "A->0", "B->0", "C->0", "D->0",
#'           "A->B", "A->C", "C->D") ## set of reactions
#' nC <- 3 ## number of clones
#' S <- 100 ## trajectory length
#' tau <- 1 ## for tau-leaping algorithm
#' u_1 <- c(.2, .15, .17, .09*5,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' u_2 <- c(.2, .15, .17, .09,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' u_3 <- c(.2, .15, .17*3, .09,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' theta_allcls <- cbind(u_1, u_2, u_3) ## clone-specific parameters
#' rownames(theta_allcls) <- colnames(V)
#' s20 <- 1 ## additional noise
#' Y <- array(data = NA,
#'            dim = c(S + 1, nrow(V), nC),
#'            dimnames = list(seq(from = 0, to = S*tau, by = tau),
#'                            rownames(V),
#'                            1:nC)) ## empty array to store simulations
#' Y0 <- c(100,0,0,0) ## initial state
#' names(Y0) <- rownames(V)
#' for (cl in 1:nC) { ## loop over clones
#'   Y[,,cl] <- get.sim.tl(Yt = Y0,
#'                         theta = theta_allcls[,cl],
#'                         S = S,
#'                         s2 = s20,
#'                         tau = tau,
#'                         rct.lst = rcts)
#' }
#' null.res <- fit.null(Y = Y, rct.lst = rcts) ## null model fitting
#'
#' re.res <- fit.re(theta_0 = null.res$fit$par,
#'                      Y = Y,
#'                      rct.lst = rcts,
#'                      maxemit = 100) ## random-effects model fitting
#'
#' get.scatterpie(re.res, txt = TRUE)
#' }
#' @export
get.scatterpie <- function(re.res, txt = FALSE, legend = FALSE){
  V <- re.res$design$V ## net-effect matrix
  M <- re.res$design$M ## base-model design matrix

  linColMap <- brewer.pal(nrow(V), "Set2") ## color palette
  names(linColMap) <- rownames(V)

  K <- ncol(re.res$design$M_bdiag)/ncol(re.res$design$M) # dim(Y)[3] ## number of clones

  M_bdiag <- re.res$design$M_bdiag ## random-effects design matrix
  VEuy_sol <- re.res$fit$VEuy
  euy_sol <- VEuy_sol[["euy"]] ## conditional expectaction E[u|y] of the random effects u given the data y

  euy_allIS <- matrix(data = euy_sol, nrow = K, byrow = TRUE)
  rownames(euy_allIS) <- unique(rownames(M))
  colnames(euy_allIS) <- colnames(V)
  t(euy_allIS)

  d <- as.data.frame(euy_allIS[,1:(2*nrow(V))] %*% t(V[,1:(2*nrow(V))]))
  d[d < 0] <- 0
  n <- nrow(d)
  A <- runif(n, 0, 1)
  B <- runif(n, 0, 2*pi)
  diameter <- sqrt(sum(sapply(apply(euy_allIS %*% t(V), 1, norm2), function(d){
    d[d < 0] <- 0
    return(pi*(d/2)^2)
    }))/pi)*2
  lat <- sqrt(A)*cos(B)*diameter
  long <- sqrt(A)*sin(B)*diameter
  d[, c("long","lat")] <- cbind(long, lat)

  cellTypes <- rownames(V)
  n <- nrow(d)
  d$region <- factor(1:n)
  d$radius <- apply(d[,colnames(d) %in% cellTypes], 1, norm2)
  d <- d[, c("long","lat","radius", "region", cellTypes)]

  XLIM <- (max(d$lat) + max(d$radius))
  YLIM <- (max(d$long) + max(d$radius))
  XLGEGEND <- YLEGEND <- 100

  piechart <- (ggplot() + scatterpie::geom_scatterpie(aes(x=long, y=lat, group=region, r=radius),
                                                      data=d,
                                                      cols=cellTypes,
                                                      color="white")
               + scale_fill_manual(values=alpha(linColMap, alpha = .5)) + coord_equal()
               + xlab("") + ylab("")
               + theme(axis.line=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks=element_blank(),
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank(),
                       panel.background=element_blank(),
                       panel.border=element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       plot.background=element_blank()))
  if(txt){
    piechart <- piechart + geom_text(aes(x=long, y=lat, group=region, label = region, fontface = "bold"),
                                     data = d, size = 10)
  }
  if(legend){
    piechart <- piechart + geom_scatterpie_legend(d$radius,
                                                  x=-XLIM,
                                                  y=-YLIM)
  }
  return(piechart)
}

#' Clonal boxplots
#'
#' Draw clonal boxplots of a random-effects reaction network.
#' @param re.res output list returned by fit.re().
#' @details This function generates the boxplots of the conditional expectations
#' \deqn{w_k = E_{u\vert \Delta Y; \hat{\psi}}[u^k_{\alpha_{l}}] - E_{u\vert \Delta Y; \hat{\psi}}[u^k_{\delta_{l}}]}{Eu|y;ψ[α_l] - Eu|y;ψ[δ_l]},
#' computed from the estimated parameters \eqn{\hat{\psi}}{ψ} for the clone-specific net-duplication in each cell lineage l (different colors).
#' The whiskers extend to the data extremes.
#' @examples
#' \donttest{
#' library(RestoreNet)
#' library(ggplot2)
#' library(scatterpie)
#' rcts <- c("A->1", "B->1", "C->1", "D->1",
#'           "A->0", "B->0", "C->0", "D->0",
#'           "A->B", "A->C", "C->D") ## set of reactions
#' nC <- 3 ## number of clones
#' S <- 100 ## trajectory length
#' tau <- 1 ## for tau-leaping algorithm
#' u_1 <- c(.2, .15, .17, .09*5,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' u_2 <- c(.2, .15, .17, .09,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' u_3 <- c(.2, .15, .17*3, .09,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' theta_allcls <- cbind(u_1, u_2, u_3) ## clone-specific parameters
#' rownames(theta_allcls) <- colnames(V)
#' s20 <- 1 ## additional noise
#' Y <- array(data = NA,
#'            dim = c(S + 1, nrow(V), nC),
#'            dimnames = list(seq(from = 0, to = S*tau, by = tau),
#'                            rownames(V),
#'                            1:nC)) ## empty array to store simulations
#' Y0 <- c(100,0,0,0) ## initial state
#' names(Y0) <- rownames(V)
#' for (cl in 1:nC) { ## loop over clones
#'   Y[,,cl] <- get.sim.tl(Yt = Y0,
#'                         theta = theta_allcls[,cl],
#'                         S = S,
#'                         s2 = s20,
#'                         tau = tau,
#'                         rct.lst = rcts)
#' }
#' null.res <- fit.null(Y = Y, rct.lst = rcts) ## null model fitting
#'
#' re.res <- fit.re(theta_0 = null.res$fit$par,
#'                      Y = Y,
#'                      rct.lst = rcts,
#'                      maxemit = 100) ## random-effects model fitting
#'
#' get.boxplots(re.res)
#' }
#' @export
get.boxplots <- function(re.res){
  K <- ncol(re.res$design$M_bdiag)/ncol(re.res$design$M)
  euy_allIS <- matrix(data = re.res$fit$VEuy$euy, nrow = K, byrow = TRUE)
  rownames(euy_allIS) <- unique(rownames(re.res$design$M))
  colnames(euy_allIS) <- colnames(re.res$design$V)

  linColMap <- brewer.pal(nrow(re.res$design$V), "Set2") ## color palette
  names(linColMap) <- rownames(re.res$design$V)


  bp <- boxplot(euy_allIS %*% t(re.res$design$V),
                range = 0, lwd = 3,
                col = linColMap, border = linColMap,
                cex.axis = 2, cex.lab = 2,
                horizontal = TRUE,
                ylab = ""
                # ylab = expression('E'['u|y'] (u[alpha]^k) - 'E'['u|y'] (u[delta]^k))
  )
  return(bp)
}


#' \eqn{\tau}{tau}-leaping simulation algorithm
#'
#' Simulate a trajectory of length S for a stochastic reaction network.
#' @param Yt starting point of the trajectory
#' @param theta vector parameter for the reactions.
#' @param V A \eqn{p \times K}{p by K} dimensional net-effect matrix.
#' @param S length of the simulated trajectory.
#' @param s2 noise variance (defaults to 0).
#' @param tau time interval length (defaults to 1).
#' @details This function allows to simulate a trajectory of a single clone given an
#' initial conditions \eqn{Y_0}{Y_0} for the cell counts, and obeying to a particular cell differentiation network
#' defined by a net-effect (stoichiometric) matrix \eqn{V}{V} and an hazard function \eqn{h()}{h()}.
#' The function allows to consider only three cellular events, such as
#' cell duplication (\eqn{Y_{it} \rightarrow 1}{Y_{it} -> 1}), cell death (\eqn{Y_{it} \rightarrow \emptyset}{Y_{it} -> 0})
#' and cell differentiation (\eqn{Y_{it} \rightarrow Y_{jt}}{Y_{it} -> Y_{jt}}) for a clone-specific time counting process
#' \deqn{Y_t = (Y_{1t},\dots,Y_{Nt})}{Y_t = (Y_{1t},...,Y_{Nt})}
#' observed in $N$ distinct cell lineages.
#' In particular,  the cellular events of duplication, death and differentiation are
#' respectively coded with the character labels \eqn{\texttt{"A->1"}}{"A->1"}, \eqn{\texttt{"A->0"}}{"A->0"},
#' and \eqn{\texttt{"A->B"}}{"A->B"}, where  \eqn{\texttt{A}}{A} and  \eqn{\texttt{B}}{B} are two distinct cell types.
#' The output is a \eqn{3}{3}-dimensional array \eqn{Y}{Y} whose \eqn{ijk}{ijk}-entry \eqn{Y_{ijk}}{Y_{ijk}}
#' is the number of cells of clone \eqn{k}{k} for cell type \eqn{j}{j} collected at time \eqn{i}{i}.
#' More mathematical details can be found in the vignette of this package.
#' @return A \eqn{S \times p}{S by p} dimensional matrix of the simulated trajectory.
#' @examples
#' \donttest{
#' library(RestoreNet)
#' rcts <- c("A->1", "B->1", "C->1", "D->1",
#'           "A->0", "B->0", "C->0", "D->0",
#'           "A->B", "A->C", "C->D") ## set of reactions
#' nC <- 3 ## number of clones
#' S <- 100 ## trajectory length
#' tau <- 1 ## for tau-leaping algorithm
#' u_1 <- c(.2, .15, .17, .09*5,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' u_2 <- c(.2, .15, .17, .09,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' u_3 <- c(.2, .15, .17*3, .09,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' theta_allcls <- cbind(u_1, u_2, u_3) ## clone-specific parameters
#' rownames(theta_allcls) <- colnames(V)
#' s20 <- 1 ## additional noise
#' Y <- array(data = NA,
#'            dim = c(S + 1, nrow(V), nC),
#'            dimnames = list(seq(from = 0, to = S*tau, by = tau),
#'                            rownames(V),
#'                            1:nC)) ## empty array to store simulations
#' Y0 <- c(100,0,0,0) ## initial state
#' names(Y0) <- rownames(V)
#' for (cl in 1:nC) { ## loop over clones
#'   Y[,,cl] <- get.sim.tl(Yt = Y0,
#'                         theta = theta_allcls[,cl],
#'                         S = S,
#'                         s2 = s20,
#'                         tau = tau,
#'                         rct.lst = rcts)
#' }
#' }
#' @export
get.sim.tl <- function(Yt, theta, S, s2 = 0, tau = 1, rct.lst){
  # ct.lst <- unique(setdiff(c(sapply(rct.lst, function(r){
  #   as.vector(unlist(strsplit(r, split = "->", fixed = T)))
  # }, simplify = "array")), c("0", "1")))

  # pattern <- "[A-Z]+[0-9]*->([A-Z]+[0-9]*|[0-1]+)"
  pattern <- "[A-Z0-9]{1,}->[A-Z0-9]{0,}[0-1]{0,}"

  if(sum(grepl(pattern,
               rct.lst,
               ignore.case = FALSE)) < length(rct.lst)){
    cat(paste("Unvalid reactions provided in 'rct.lst': ",
              paste(rct.lst[which(!grepl(pattern,
                                         c(rct.lst),
                                         ignore.case = FALSE))], collapse = ", "), sep = ""))
    return()
  }

  # ct.lst <- unique(setdiff(c(sapply(rct.lst, function(r){
  #   as.vector(unlist(strsplit(r, split = "->", fixed = T)))
  # }, simplify = "array")), c("0", "1")))
  # ct.lst <- setdiff(as.vector(unlist(c(sapply(rct.lst, function(r){
  #   as.vector(unlist(strsplit(r, split = "->", fixed = T)))
  # }, simplify = "array")))), c("0","1"))
  ct.lst <- setdiff(as.vector(unlist(c(sapply(rct.lst, function(r){
    as.vector(unlist(str_split(r, "->")))
  }, simplify = "array")))), c("0","1"))

  if(sum(!(ct.lst %in% colnames(Y))) > 0){
    cat(paste("Cell types not present in 'Y': ",
              paste(as.vector(unlist(lapply(ct.lst[which(!(ct.lst %in% colnames(Y)))],
                                            function(l){paste("'", l, "'", sep = "")}))), collapse = ", "),
              "\nPlease provide compatible cell types in 'Y' and 'rct.lst'.",
              sep = ""))
    return()
  }

  compile.h(rct.lst = rct.lst, envir = environment())
  V <- get.V(ct.lst = ct.lst, rct.lst = rct.lst)

  Y0 <- Yt
  Ysim <- matrix(data = NA, nrow = S, ncol = nrow(V))
  colnames(Ysim) <- names(Yt)
  for (s in 1:S) {
    if(s %% 100 == 0) cat(s, "\n") ## massage to keep track of the status

    th_tr <- tau*theta
    N_tr <- as.numeric(get.h(Yt, th_tr))

    Yt <- as.numeric(Yt + V %*% N_tr) ## update the molecule counts
    names(Yt) <- head(LETTERS,4)

    Ysim[s,] <- Yt ## store the current simulated counts
  }

  Ysim <- rbind(Y0, Ysim)
  Ysim <- Ysim + rnorm(prod(dim(Ysim)), mean = 0, sd = sqrt(s2))
  Ysim[Ysim < 0] <- 0
  rownames(Ysim) <- seq(from = 0, to = S*tau, by = tau)

  return(Ysim)
}

#' Net-effect matrix
#'
#' This function builds the net-effect matrix V.
#' @param ct.lst list of cell types
#' @param rct.lst list of biochemical reactions.
#' A differentiation move from cell type "A" to cell type "B" must be coded as "A->B"
#' Duplication of cell "A" must be coded as "A->1"
#' Death of cell "A" must be coded as "A->0"
#' @return The net-effect matrix V
#' @keywords internal
get.V <- function(ct.lst, rct.lst){
  V <- matrix(data = 0, nrow = length(ct.lst), ncol = length(rct.lst))
  rownames(V) <- ct.lst
  colnames(V) <- rct.lst

  for (r in rct.lst) {
    rgts_prod <- as.vector(unlist(strsplit(r, split = "->", fixed = T)))
    V[rgts_prod[1],r] <- -1
    if(rgts_prod[2] != "0"){
      if(rgts_prod[2] == "1"){
        V[rgts_prod[1],r] <- 1
      }else{
        V[rgts_prod[2],r] <- 2
      }
    }
  }
  return(V)
}

#' Generate function for hazard
#'
#' This function dynamically builds the function get.h() that will be used to get the hazard vector.
#' @param rct.lst list of biochemical reactions.
#' A differentiation move from cell type "A" to cell type "B" must be coded as "A->B"
#' Duplication of cell "A" must be coded as "A->1"
#' Death of cell "A" must be coded as "A->0"
#' @return The function doesn't return anything.
#' @keywords internal
compile.h <- function(rct.lst, envir){
  constr.lst <- c()
  get.h.string <- paste("get.h <- function(x, theta){
                h <- c(",
                        paste(sapply(rct.lst, function(r){
                          rgts <- as.vector(unlist(strsplit(r, split = "->", fixed = T)))[1]
                          prds <- as.vector(unlist(strsplit(r, split = "->", fixed = T)))[2]
                          if(prds == "0"){
                            hx <- paste("x['", rgts, "']^2*theta['", r, "']", sep = "")
                          }else{
                            hx <- paste("x['", rgts, "']*theta['", r, "']", sep = "")
                          }

                        }, simplify = "array"), collapse = ",\n"),
                        ")
return(h)
}", sep = "")

  for (constr in constr.lst) {
    constr <- as.vector(unlist(strsplit(constr, split = "=", fixed = TRUE)))
    pattern <- constr[1]
    replacement <- constr[2]

    get.h.string <- gsub(pattern = pattern,
                         replacement = replacement,
                         x = get.h.string)
  }

  eval(parse(text=get.h.string), envir = envir) # envir = .GlobalEnv
}

#' Random generator from a Multivariate Normal distribution
#'
#' This function generate a random vector from a multivariate normal distribution.
#' @param d sample d.
#' @param p vector dimension.
#' @param Mu mean vector.
#' @param Sigma covariance matrix.
#' @return A \eqn{p \times d}{p by d} matrix. Each column is a p-variate random vector from
#' a multivariate normal with mean vector Mu and covariance matrix Sigma.
#' @keywords internal
rmvNorm <- function(d, p, Mu, Sigma){
  Z <- matrix(rnorm(d*p), nrow = d, ncol = p)
  X <- t(Z%*%chol(Sigma)) + matrix(Mu, ncol = 1)[,rep(1,d)]

  return(X)
}

#' Time-adjacent increments of the cell counts
#'
#' This function generates time-adjacent increments from a cell counts matrix y.
#' @param y clone-specific \eqn{t \times p}{t by p} cell count matrix, where t is the number of time-points and p the number of cell types.
#' @return A \eqn{(t-1) \times p}{(t - 1) by p} dimensional vector containing the time-adjacent increments.
#' @keywords internal
get.dx <- function(y){
  nCells <- ncol(y)
  T<-nrow(y)
  compdat<-t(y)
  DX <- as.matrix(compdat[,-1]-compdat[,-T])
  colnames(DX) <- as.character(apply(DX, 2, function(dx){sum(dx != 0)>0}))
  colnames(DX) <- substr(x = colnames(DX), start = 1, stop = 1)
  dx <- as.vector(DX)
  names(dx) <- rep(colnames(DX), each = nCells)
  return(dx)
}

#' Design matrix M
#'
#' This function generates time-adjacent increments from a cell counts matrix y.
#' @param y clone-specific \eqn{t \times p}{t by p} cell count matrix, where t is the number of time-points and p the number of cell types.
#' @param V A \eqn{p \times K}{p by K} dimensional net-effect matrix.
#' @param dT A (t-1)-dimensional vector of the time-increments.
#' @return A \eqn{tp \times K}{t*p by K} dimensional (design) matrix.
#' @keywords internal
get.M <- function(y, V, dT, get.h){
  theta <- rep(1, ncol(V))
  names(theta) <- colnames(V)
  Mt <- lapply(matrix(1:(nrow(y))), FUN = function(t){
    Mt <- Matrix(V%*%diag(get.h(x = y[t,], theta)), sparse = TRUE);
    Mt_01 <- Mt;
    Mt_01[Mt_01>0] <- 1;
    if(sum(Mt_01 != 0)>0){
      rownames(Mt) <- rep("T", nrow(Mt));
    }else{
      rownames(Mt) <- rep("F", nrow(Mt));
    }
    return(Mt)
  })
  M <- Reduce(rbind, Mt)
  colnames(M) <- colnames(V)
  M <- M*dT
  return(M)
}

#' Stochastic covariance matrix W
#'
#' This function returns the stochastic covariance matrix W associated to a design matrix M, a net-effect matrix V, and a vector parameter theta.
#' @param M A \eqn{n \times K}{n by K} dimensional (design) matrix.
#' @param V A \eqn{p \times K}{p by K} dimensional net-effect matrix.
#' @param theta p-dimensional vector parameter.
#' @return A \eqn{n \times n}{n by n} dimensional stochastic covariance matrix.
#' @keywords internal
get.W <- function(M, V, VCNs, theta, nObs){
  p <- ncol(M)
  beta <- theta[1:p]
  NN <- nrow(M)
  nC <- nrow(V)
  nLin <- nrow(V)

  idx <- cbind(rep(1:nC, times = nC), rep(rep(1:nC, each = nC)))
  idx <- idx[rep(1:nrow(idx), NN/nC),]
  idx_unit <- matrix(data = rep((0:(NN/nC - 1))*nC, each = prod(dim(cbind(rep(1:nC, times = nC), rep(rep(1:nC, each = nC)))))), nrow = nrow(idx), ncol = 2, byrow = TRUE)
  idx <- idx + idx_unit

  nObs <- rep(nObs, times = table(rownames(M)))
  nObsMat <- matrix(data = rep(nObs, times = nLin), ncol = nLin, byrow = FALSE)

  W <- Matrix(data = 0, nrow = nrow(M), ncol = nrow(M), sparse = TRUE)
  W[idx] <- t(((M%*%Matrix(diag(beta), sparse = TRUE)%*%t(V))*nObsMat))

  diag(W) <- diag(W)*VCNs
  rownames(W) <- colnames(W) <- rownames(M)

  return(W) ## PLEASE NOTE: dT is already included in M, we don't need "W <- W*dT"
}

#' Negative log-likelihood of the base (null) model
#'
#' Compute the negative log-likelihood of the null model,
#' namely the linear model without rondom effects at all.
#' @param theta p-dimensional vector parameter.
#' @param M A \eqn{n \times K}{n by K} dimensional (design) matrix.
#' @param y n-dimensional vector of the time-adjacent cellular increments
#' @param V A \eqn{p \times K}{p by K} dimensional net-effect matrix.
#' @param dW p-dimensional list of the partial derivatives of W w.r.t. theta.
#' @return the value of the negative log-likelihood.
#' @keywords internal
get.nl  <- function(theta, M, y, V, VCNs, nObs, dW){
  N <- nrow(M)
  p <- ncol(M)
  b <- head(theta, p)
  s2 <- tail(theta, 1)

  W <- get.W(M, V, VCNs, theta, nObs)
  S <- W
  diag(S) <- diag(S) + s2
  Si <- chol2inv(chol(S))
  nll_temp <- as.numeric(ldet(S) + t(y - M%*%b) %*% Si %*% (y - M%*%b))

  gcRes <- gc()
  return(nll_temp)
}

#' Gradient of the negative log-likelihood of the base (null) model
#'
#' Compute the gradient of the negative log-likelihood of the null model,
#' namely the linear model without rondom effects at all.
#' @param theta p-dimensional vector parameter.
#' @param M A \eqn{n \times K}{n by K} dimensional (design) matrix.
#' @param y n-dimensional vector of the time-adjacent cellular increments
#' @param V A \eqn{p \times K}{p by K} dimensional net-effect matrix.
#' @param dW p-dimensional list of the partial derivatives of W w.r.t. theta.
#' @return p-dimensional vector of the gradient of the negative log-likelihood.
#' @keywords internal
get.gnl <- function(theta, M, y, V, VCNs, nObs, dW){
  N <- nrow(M)
  p <- ncol(M)
  b <- head(theta, p)
  s2 <- tail(theta, 1)

  W <- get.W(M, V, VCNs, theta, nObs)
  S <- W
  diag(S) <- diag(S) + s2
  Si <- chol2inv(chol(S))

  gnq <- rep(NA, p + 1)
  gnq[1:p] <-  -2*as.numeric(t(y)%*%Si%*%M)  + 2*as.numeric(t(M) %*% Si%*%M%*%b)
  for (j in 1:p) {
    djSi <- -Si%*%dW[[j]]%*%Si
    djldetSi <- sum(t(Si)*dW[[j]])
    gnq[j] <- gnq[j] + t(y)%*%djSi%*%y - 2*t(y)%*%djSi%*%M%*%b + t(b)%*%t(M)%*%djSi%*%M%*%b + djldetSi #
  }
  gnq[p+1] <-  - t(y)%*%Si%*%Si%*%y + 2*t(b)%*%t(M)%*%Si%*%Si%*%y - t(b)%*%t(M)%*%Si%*%Si%*%M%*%b + tr(Si) #

  gcRes <- gc()
  return(gnq)
}

#' Fit base model
#'
#' Fit the base (null) model to the given data using a maximum likelihood approach.
#' @param theta_start p-dimensional vector parameter used as initial guess in the inference procedure.
#' @param M A \eqn{n \times K}{n by K} dimensional (design) matrix.
#' @param y n-dimensional vector of the time-adjacent cellular increments
#' @param V A \eqn{p \times K}{p by K} dimensional net-effect matrix.
#' @param psiLB p-dimensional vector of lower bound values for theta.
#' @param psiUB p-dimensional vector of upper bound values for theta.
#' @param maxit maximum number of iterations for the optimization step.
#' This argument is passed to optim() function. Details on "maxit" can be found in "optim()" documentation page.
#' @param factr controls the convergence of the "L-BFGS-B" method.
#' Convergence occurs when the reduction in the objective is within this factor of the machine tolerance.
#' Default is 1e7, that is a tolerance of about 1e-8.
#' This argument is passed to optim() function.
#' @param pgtol helps control the convergence of the "L-BFGS-B" method.
#' It is a tolerance on the projected gradient in the current search direction.
#' This defaults to zero, when the check is suppressed.
#' This argument is passed to optim() function.
#' @param lmm is an integer giving the number of BFGS updates retained in the "L-BFGS-B" method, It defaults to 5.
#' This argument is passed to optim() function.
#' @param VCNs A n-dimensional vector including values of the vector copy number corresponding to the cell counts of y.
#' @param nObs A K-dimensional vector including the frequencies of each clone k (\eqn{k = 1,\dots,K}{k = 1,..,K}).
#' @param trace Non-negative integer. If positive, tracing information on the progress of the optimization is produced.
#' This parameter is also passed to the optim() function.
#' Higher values may produce more tracing information: for method "L-BFGS-B" there are six levels of tracing.
#' (To understand exactly what these do see the source code: higher levels give more detail.)
#' @return A 3-length list. First element is the output returned by "optim()" function (see "optim()" documentation for details).
#' Second element is a vector of statistics associated to the fitted null model:
#' \itemize{
#'  \item{"nPar"}{number of parameters of the base(null) model}
#'  \item{"cll"}{value of the conditional log-likelihood, in this case just the log-likelihood}
#'  \item{"mll"}{value of the marginal log-likelihood, in this case just the log-likelihood}
#'  \item{"cAIC"}{conditional Akaike Information Criterion (cAIC), in this case simply the AIC.}
#'  \item{"mAIC"}{marginal Akaike Information Criterion (mAIC), in this case simply the AIC.}
#'  \item{"Chi2"}{value of the \eqn{\chi^2}{Chi-squared} statistic \eqn{(y - M\theta)'S^{-1}(y - M\theta)}{(y - Mθ)'S^-1(y - Mθ)}.}
#'  \item{"p-value"}{p-value of the \eqn{\chi^2}{Chi-squared} test for the null hypothesis that Chi2
#'  follows a \eqn{\chi^2}{Chi-squared} distribution with n - nPar degrees of freedom.}
#' }
#' The third element, called "design", is a list including:
#' \itemize{
#'  \item{"M"}{A \eqn{n \times K}{n by K} dimensional (design) matrix.}
#'  \item{"V"}{A \eqn{p \times K}{p by K} dimensional net-effect matrix.}
#' }
#' @keywords internal
nullModelFitting <- function(theta_start,
                             M,
                             y,
                             V,
                             psiLB,
                             psiUB,
                             maxit,
                             factr,
                             pgtol,
                             lmm,
                             VCNs,
                             nObs,
                             trace = TRUE){

  p <- ncol(M)
  N <- nrow(M)
  dW <- apply(diag(1, nrow = p, ncol = p), 2, function(j){get.W(M, V, VCNs, j, nObs)})

  res <- list();

  res_null <- try(optim(par = theta_start,
                        fn = get.nl,
                        gr = get.gnl,
                        M = M, y = y, V = V, VCNs = VCNs, nObs = nObs, dW = dW,
                        method = "L-BFGS-B",
                        lower = psiLB,
                        upper = psiUB,
                        control = list(maxit = maxit,
                                       factr = factr,
                                       pgtol = pgtol,
                                       lmm = lmm,
                                       trace = trace)), silent = TRUE);

  res$fit <- res_null


  if(!inherits(res_null, "try-error")){

    theta_res <- res_null$par
    b_res <- head(theta_res, p)
    s2_res <- tail(theta_res, 1)

    W_res <- get.W(M, V, VCNs, theta_res, nObs)
    S_res <- W_res
    diag(S_res) <- diag(S_res) + s2_res
    Si_res <- chol2inv(chol(S_res))

    ll_res <- dmvNorm(x = y, Mu = M %*% b_res, Sigma = S_res, SigmaInv = Si_res) ## marginal distribution p(y)
    nll_res <- -ll_res

    p_null <- length(theta_res)
    AIC_null <- 2*nll_res + 2*p_null ## AIC
    AICc_null <- 2*nll_res + 2*(p_null + 1)*N/(N - p_null - 2) ## corrected AIC
    BIC_null <- log(N)*p_null + 2*nll_res ## BIC

    chi2Stat <- as.numeric(t(y - M%*%b_res)%*%S_res%*%(y - M%*%b_res))
    chi2pVal <- 1 - pchisq(q = chi2Stat, df = N - p_null)

    res$stats <- c(as.numeric(p_null), as.numeric(ll_res), as.numeric(ll_res), as.numeric(AIC_null), as.numeric(AICc_null), as.numeric(chi2Stat), as.numeric(chi2pVal))
    names(res$stats) <- c("nPar","cll","mll","cAIC","mAIC","Chi2","p-value")
  }

  design <- list()
  design$M <- M
  design$V <- V
  res$design <- design

  return(res)
}

#' Fit random-effects model
#'
#' Fit the random-effects model to the given data using an expectation-maximization algorithm.
#' @param theta_0 p-dimensional vector parameter used as initial guess in the inference procedure.
#' @param V A \eqn{p \times K}{p by K} dimensional net-effect matrix.
#' @param M A \eqn{n \times K}{n by K} dimensional (design) matrix.
#' @param M_bdiag A\eqn{n \times Jp}{n by J*p} dimensional block-diagonal design matrix.
#' Each j-th block (\eqn{j = 1,\dots,J}{j = 1,..,J}) is a \eqn{n_j \times p}{n_j by p} dimensional design matrix for the j-th clone.
#' @param y n-dimensional vector of the time-adjacent cellular increments
#' @param VCNs A n-dimensional vector including values of the vector copy number corresponding to the cell counts of y.
#' @param nObs A K-dimensional vector including the frequencies of each clone k (\eqn{k = 1,\dots,K}{k = 1,..,K}).
#' @param maxit maximum number of iterations for the optimization step.
#' This argument is passed to optim() function. Details on "maxit" can be found in "optim()" documentation page.
#' @param maxemit maximum number of iterations for the expectation-maximization algorithm.
#' @param eps relative error for the value x and the objective function f(x) that has to be optimized
#' in the expectation-maximization algorithm.
#' @param thetaLB p-dimensional vector of lower bound values for theta.
#' @param thetaUB p-dimensional vector of upper bound values for theta.

#' @param factr controls the convergence of the "L-BFGS-B" method.
#' Convergence occurs when the reduction in the objective is within this factor of the machine tolerance.
#' Default is 1e7, that is a tolerance of about 1e-8.
#' This argument is passed to optim() function.
#' @param pgtol helps control the convergence of the "L-BFGS-B" method.
#' It is a tolerance on the projected gradient in the current search direction.
#' This defaults to zero, when the check is suppressed.
#' This argument is passed to optim() function.
#' @param lmm is an integer giving the number of BFGS updates retained in the "L-BFGS-B" method, It defaults to 5.
#' This argument is passed to optim() function.

#' @param trace Non-negative integer. If positive, tracing information on the progress of the optimization is produced.
#' This parameter is also passed to the optim() function.
#' Higher values may produce more tracing information: for method "L-BFGS-B" there are six levels of tracing.
#' (To understand exactly what these do see the source code: higher levels give more detail.)
#' @return The output returned by "optim()" function (see "optim()" documentation for details) along with
#' the conditional expectation \eqn{E[u \vert y]}{E[u|y]} and variance \eqn{V[u \vert y]}{V[u|y]}
#' of the latent states u given the observed states y from the last step of the expectation-maximization algorithm.
#' @keywords internal
rndEffModelFitting <- function(theta_0,
                               V,
                               M,
                               M_bdiag,
                               y,
                               VCNs,
                               nObs,
                               maxit,
                               maxemit,
                               eps = 1e-5,
                               thetaLB,
                               thetaUB,
                               factr,
                               pgtol,
                               lmm,
                               trace = TRUE){

  ## for EM:
  nIT <- 0 ## number of iterations
  deltaQ <- 10^3 ## |Q(theta_new) - Q(theta_curent)|
  deltaTheta <- 10^3  ## |theta_new - theta_curent|

  p <- ncol(M)
  theta_curr <- theta_start <- rep(1, 2*p + 1)
  theta_curr[1:p] <- theta_start[1:p] <- head(theta_0, -1)

  VEuy_curr <- VEuy(theta_curr, M, M_bdiag, y, V, VCNs, nObs)
  euy_curr <- VEuy_curr[["euy"]] ## conditional expectaction \eqn{E[u \vert y]}{E[u|y]} of the random (latent) effects u given the observed data y
  vuy_curr <- VEuy_curr[["vuy"]] ## conditional variance \eqn{V[u \vert y]}{V[u|y]} of the random (latent) effects u given the observed data y

  dW <- apply(diag(1, nrow = p, ncol = p), 2, function(j){get.W(M, V, VCNs, j, nObs)})

  nQ_curr <- nQ(theta_curr, euy_curr = euy_curr, vuy_curr = vuy_curr, M, M_bdiag, y, V, VCNs = VCNs, nObs = nObs, dW = dW)
  ngQ_curr <- ngQ(theta_curr, euy_curr = euy_curr, vuy_curr = vuy_curr, M, M_bdiag, y, V, VCNs = VCNs, nObs = nObs, dW = dW)

  while (nIT < maxemit & (deltaQ > eps | deltaTheta > eps)) { ## EM-algorithm
    # while (nIT < maxIT & deltaQ > eps & deltaTheta > eps) { ## EM-algorithm
    if(trace){cat(paste("EM-alg. It. n. ", nIT, "\terr_Q = ", deltaQ, ";\terr_theta = ", deltaTheta, "...\n", sep = ""))}

    res <- try(optim(par = theta_curr,
                            fn = nQ,
                            gr = ngQ,
                            euy_curr = euy_curr, vuy_curr = vuy_curr, M = M, M_bdiag = M_bdiag, y = y, V = V, VCNs = VCNs, nObs = nObs, dW = dW,
                            method = "L-BFGS-B",
                            lower = thetaLB,
                            upper = thetaUB,
                            control = list(maxit = maxit,
                                           factr = factr,
                                           pgtol = pgtol,
                                           lmm = lmm,
                                           trace = trace)), silent = TRUE);

    theta_old <- theta_curr
    theta_curr <- res$par

    nQ_old <- nQ_curr

    VEuy_curr <- VEuy(theta_curr, M, M_bdiag, y, V, VCNs, nObs)
    euy_curr <- VEuy_curr[["euy"]] ## conditional expectaction \eqn{E[u \vert y]}{E[u|y]} of the random (latent) effects u given the observed data y
    vuy_curr <- VEuy_curr[["vuy"]] ## conditional variance \eqn{V[u \vert y]}{V[u|y]} of the random (latent) effects u given the observed data y

    nQ_curr <- nQ(theta_curr, euy_curr = euy_curr, vuy_curr = vuy_curr, M, M_bdiag, y, V, VCNs = VCNs, nObs = nObs, dW = dW)

    deltaQ <- relErr(x = nQ_curr, y = nQ_old)
    deltaTheta <- relErr(x = theta_old, y = theta_curr)
    nIT <- nIT + 1
  }

  res$VEuy <- VEuy_curr
  return(res)
}

#' E-step function Q
#'
#' Negative E-step function -Q of the expectation-maximization algorithm
#' @param theta p-dimensional vector parameter.
#' @param euy_curr current value of the conditional expectation \eqn{E[u \vert y]}{E[u|y]} of u given y,
#' where u and y are the latent and observed states respectively.
#' @param vuy_curr current value of the conditional variance \eqn{V[u \vert y]}{V[u|y]} of u given y,
#' where u and y are the latent and observed states respectively.
#' @param M A \eqn{n \times K}{n by K} dimensional (design) matrix.
#' @param M_bdiag A\eqn{n \times Jp}{n by J*p} dimensional block-diagonal design matrix.
#' Each j-th block (\eqn{j = 1,\dots,J}{j = 1,..,J}) is a \eqn{n_j \times p}{n_j by p} dimensional design matrix for the j-th clone.
#' @param y n-dimensional vector of the time-adjacent cellular increments
#' @param V A \eqn{p \times K}{p by K} dimensional net-effect matrix.
#' @param VCNs A n-dimensional vector including values of the vector copy number corresponding to the cell counts of y.
#' @param nObs A K-dimensional vector including the frequencies of each clone k (\eqn{k = 1,\dots,K}{k = 1,..,K}).
#' @param dW p-dimensional list of the partial derivatives of W w.r.t. theta.
#' @return The current value of the negative E-step function -Q.
#' @keywords internal
nQ <- function(theta, euy_curr, vuy_curr, M, M_bdiag, y, V, VCNs, nObs, dW){
  p <- ncol(V)
  K <- length(unique(rownames(M)))
  b <- theta[1:p]
  d2 <- theta[(p + 1): (2*p)]
  s2 <- tail(theta, 1)

  W <- get.W(M, V, VCNs, theta = c(b,s2), nObs)
  S <- W
  diag(S) <- diag(S) + s2
  Si <- chol2inv(chol(S))

  D2 <- kronecker(X = Diagonal(n = K, x = 1), Y = Diagonal(n = p, x = d2)) ## covariance matrix of the random effects u, using the free parameters
  D2i <- kronecker(X = Diagonal(n = K, x = 1), Y = Diagonal(n = p, x = d2^-1))
  bu <- kronecker(X = Matrix(data = 1, nrow = K, ncol = 1, sparse = TRUE), Y = Matrix(data = b, ncol = 1, sparse = TRUE)) ## expected values of the random effects u, using the free parameters

  Q <- ( - .5*ldet(A = S)
         - .5*( t(y)%*%Si%*%y
                -2*t(euy_curr)%*%t(M_bdiag)%*%Si%*%y
                + sum( t(t(M_bdiag)%*%Si) * (M_bdiag%*%vuy_curr))
                + sum((t(M_bdiag)%*%Si%*%M_bdiag%*%euy_curr)*euy_curr) )
         + t(euy_curr)%*%D2i%*%bu
         -.5*t(bu)%*%D2i%*%bu
         - .5*ldet(A = D2)
         - .5*(sum(D2i*vuy_curr) + sum((D2i%*%euy_curr)*(euy_curr))) );

  Q <- as.numeric(Q)

  gcRes <- gc()
  return(-Q)
}

#' Gradient of the E-step function Q
#'
#' Gradient -\eqn{\nabla Q}{\nabla Q} of the negative E-step function -Q of the expectation-maximization algorithm
#' @param theta p-dimensional vector parameter.
#' @param euy_curr current value of the conditional expectation \eqn{E[u \vert y]}{E[u|y]} of u given y,
#' where u and y are the latent and observed states respectively.
#' @param vuy_curr current value of the conditional variance \eqn{V[u \vert y]}{V[u|y]} of u given y,
#' where u and y are the latent and observed states respectively.
#' @param M A \eqn{n \times K}{n by K} dimensional (design) matrix.
#' @param M_bdiag A\eqn{n \times Jp}{n by J*p} dimensional block-diagonal design matrix.
#' Each j-th block (\eqn{j = 1,\dots,J}{j = 1,..,J}) is a \eqn{n_j \times p}{n_j by p} dimensional design matrix for the j-th clone.
#' @param y n-dimensional vector of the time-adjacent cellular increments
#' @param V A \eqn{p \times K}{p by K} dimensional net-effect matrix.
#' @param VCNs A n-dimensional vector including values of the vector copy number corresponding to the cell counts of y.
#' @param nObs A K-dimensional vector including the frequencies of each clone k (\eqn{k = 1,\dots,K}{k = 1,..,K}).
#' @param dW p-dimensional list of the partial derivatives of W w.r.t. theta.
#' @return p-dimensional vector of the gradient -\eqn{\nabla Q}{\nabla Q} of the negative E-step function -Q.
#' @keywords internal
ngQ <- function(theta, euy_curr, vuy_curr, M, M_bdiag, y, V, VCNs, nObs, dW){
  p <- ncol(V)
  K <- length(unique(rownames(M)))
  b <- theta[1:p]
  d2 <- theta[(p + 1): (2*p)]
  s2 <- tail(theta, 1)

  W <- get.W(M, V, VCNs, theta = c(b,s2), nObs)
  S <- W
  diag(S) <- diag(S) + s2
  Si <- chol2inv(chol(S))

  D2 <- kronecker(X = Diagonal(n = K, x = 1), Y = Diagonal(n = p, x = d2)) ## covariance matrix of the random effects u, using the free parameters
  D2i <- kronecker(X = Diagonal(n = K, x = 1), Y = Diagonal(n = p, x = d2^-1)) ## covariance matrix of the random effects u, using the free parameters
  bu <- kronecker(X = Matrix(data = 1, nrow = K, ncol = 1, sparse = TRUE), Y = Matrix(data = b, ncol = 1, sparse = TRUE)) ## expected values of the random effects u, using the free parameters
  dbu <- dBu(K, p)
  dd2 <- dD2(K, p)

  ngq <- rep(NA, length(theta)) ## dQ/dtheta

  ngq[1:p] <- apply(matrix(data = 1:p, ncol = 1), MARGIN = 1,
                    FUN = function(j){ as.numeric( -.5*sum(t(Si)*dW[[j]])
                                                   -.5*( -t(y)%*%Si%*%dW[[j]]%*%Si%*%y
                                                         + 2*t(euy_curr)%*%t(M_bdiag)%*%Si%*%dW[[j]]%*%Si%*%y
                                                         + sum( t(t(M_bdiag)%*%(-Si%*%dW[[j]]%*%Si))*(M_bdiag%*%vuy_curr)  )
                                                         + sum(  (t(M_bdiag)%*%(-Si%*%dW[[j]]%*%Si)%*%M_bdiag%*%euy_curr)  *  euy_curr  ))
                                                   + t(euy_curr)%*%D2i%*%dbu[[j]]
                                                   -t(bu)%*%D2i%*%dbu[[j]] )}) ## dQ/dbeta

  ngq[(p + 1): (2*p)] <- 1*apply(matrix(data = 1:p, ncol = 1), MARGIN = 1,
                                 FUN = function(j){ as.numeric( -.5*sum(solve(D2)*t(dd2[[j]]))
                                                                -.5*(-sum(  t(solve(D2)%*%dd2[[j]]%*%solve(D2))*vuy_curr  )
                                                                     - sum(  (solve(D2)%*%dd2[[j]]%*%solve(D2)%*%euy_curr)*(euy_curr)  ))
                                                                -t(euy_curr)%*%solve(D2)%*%dd2[[j]]%*%solve(D2)%*%bu
                                                                + .5*t(bu)%*%solve(D2)%*%dd2[[j]]%*%solve(D2)%*%bu   )}) ## dQ/dd2

  ngq[2*p + 1] <- as.numeric(-.5*tr(Si)
                             -.5*(-t(y)%*%Si%*%Si%*%y
                                  + 2*t(euy_curr)%*%t(M_bdiag)%*%Si%*%Si%*%y
                                  + (-sum(  t(t(M_bdiag)%*%Si%*%Si)*(M_bdiag%*%(vuy_curr))  )
                                     - sum(  (t(M_bdiag)%*%Si%*%Si%*%M_bdiag%*%euy_curr)*euy_curr  ))))

  gcRes <- gc()
  return(-ngq)
}

#' Gradient of \eqn{E[u \vert \theta]}{E[u|θ]}
#'
#' Gradient of the expected values of the random effects u w.r.t. the free parameters
#' @param K number of clones being analyzed.
#' @param p number of free parameters.
#' @return the p-dimensional gradient of u.
#' @keywords internal
dBu <- function(K, p){
  dbu <- apply(diag(1, p, p), 2, function(cc){kronecker(X = Matrix(data = 1, nrow = K, ncol = 1, sparse = TRUE),
                                                        Y = Matrix(data = cc, ncol = 1, sparse = TRUE))})
  gcRes <- gc()
  return(dbu)
}

#' Gradient of \eqn{V[u \vert \theta]}{V[u|θ]}
#'
#' Gradient of the covariance matrix of the random effects u w.r.t. the free parameters
#' @param K number of clones being analyzed.
#' @param p number of free parameters.
#' @return the gradient of the covariance matrix of u.
#' @keywords internal
dD2 <- function(K, p){

  dd2 <- apply(diag(1, p, p), 2, function(cc){kronecker(X = Diagonal(n = K, x = 1),
                                                        Y = Diagonal(n = p, x = cc))})
  gcRes <- gc()
  return(dd2)
}

#'  \eqn{E[u \vert y]}{E[u|y]} and \eqn{V[u \vert y]}{V[u|y]}
#'
#' Conditional expectation \eqn{E[u \vert y]}{E[u|y]} and variance \eqn{V[u \vert y]}{V[u|y]} of the latent states u given the observed states y
#' @param theta_curr current p-dimensional vector parameter.
#' @param M A \eqn{n \times K}{n by K} dimensional (design) matrix.
#' @param M_bdiag A\eqn{n \times Jp}{n by J*p} dimensional block-diagonal design matrix.
#' Each j-th block (\eqn{j = 1,\dots,J}{j = 1,..,J}) is a \eqn{n_j \times p}{n_j by p} dimensional design matrix for the j-th clone.
#' @param y n-dimensional vector of the time-adjacent cellular increments
#' @param V A \eqn{p \times K}{p by K} dimensional net-effect matrix.
#' @param VCNs A n-dimensional vector including values of the vector copy number corresponding to the cell counts of y.
#' @param nObs A K-dimensional vector including the frequencies of each clone k (\eqn{k = 1,\dots,K}{k = 1,..,K}).
#' @return the conditional expectation \eqn{E[u \vert y]}{E[u|y]} and variance \eqn{V[u \vert y]}{V[u|y]} of the latent states u given the observed states y.
#' @keywords internal
VEuy <- function(theta_curr, M, M_bdiag, y, V, VCNs, nObs){
  p <- ncol(V)
  K <- length(unique(rownames(M)))
  beta_curr <- theta_curr[1:p]
  d2_curr <- theta_curr[(p + 1): (2*p)]
  s2_curr <- tail(theta_curr, 1)

  W_curr <- get.W(M, V, VCNs, theta = c(beta_curr,s2_curr), nObs)
  S_curr <- W_curr
  diag(S_curr) <- diag(S_curr) + s2_curr
  Si_curr <- chol2inv(chol(S_curr))

  D2_curr <- kronecker(X = Diagonal(n = K, x = 1), Y = Diagonal(n = p, x = d2_curr))
  D2i_curr <- kronecker(X = Diagonal(n = K, x = 1), Y = Diagonal(n = p, x = d2_curr^-1))
  bu_curr <- kronecker(X = Matrix(data = 1, nrow = K, ncol = 1, sparse = TRUE), Y = Matrix(data = beta_curr, ncol = 1, sparse = TRUE))

  vuy <- chol2inv(chol(t(M_bdiag)%*%Si_curr%*%M_bdiag + D2i_curr))
  euy <- vuy %*% ( t(M_bdiag)%*%Si_curr%*%y + D2i_curr%*%bu_curr)

  euy[euy<0] <- 0

  veuy <- list()
  veuy$euy <- euy
  veuy$vuy <- vuy

  gcRes <- gc()
  return(veuy)
}

#' Multivariate normal density function
#'
#' Proability density function of a n-variate normal distribution
#' @param x n-dimensional vector.
#' @param Mu n-dimensional mean vector.
#' @param Sigma n-dimensional covariance matrix.
#' @param SigmaInv n-dimensional precision matrix.
#' @return the probability density of a n-variate normal distribution with mean vector Mu and covariance Sigma evaluated at x.
#' @keywords internal
dmvNorm <- function(x, Mu, Sigma, SigmaInv){
  n <- length(Mu)
  res <- as.numeric(-n/2*log(2*pi) - .5*ldet(Sigma) - .5*t(x - Mu)%*%SigmaInv%*%(x - Mu))
  return(res)
}

#' Matrix trace
#'
#' Trace of a matrix
#' @param M a matrix.
#' @return the trace of M.
#' @keywords internal
tr <- function(M){
  sum(diag(M))
}

#' Euclidean 2-norm
#'
#' Compute the euclidean 2-norm of a vector x.
#' @param x a vector.
#' @return the euclidean 2-norm of x.
#' @export
norm2 <- function(x){
  return(sqrt(sum((x)^2)))
}

#' Multivariate relative error
#'
#' Compute the multivariate relative error between two vectors.
#' @param x a vector.
#' @param y a vector.
#' @return the multivariate relative error between x and y.
#' @export
relErr <- function(x,y){
  norm2(x - y)/min(norm2(x), norm2(y))
}

#' Matrix determinant
#'
#' Compute the determinant of a matrix.
#' @param A a matrix
#' @param logarithm logical; if TRUE (default) return the logarithm of the modulus of the determinant.
#' @return the determinant of A.
#' @keywords internal
ldet <- function(A, logarithm = TRUE){
  if(logarithm){
    detA <- as.numeric(determinant(x = A, logarithm = TRUE)$modulus)
  }else{
    detA <- as.numeric(determinant(x = A, logarithm = FALSE)$modulus)
  }
  return(detA)
}

#' Base and random-effects model statistics
#'
#' Main statistics from the fitted base and random-effects models
#' @param theta_null the estimated p-dimensional vector parameter for the base (null) model.
#' @param theta_rndEff the estimated p-dimensional vector parameter for the random-effects model.
#' @param M A \eqn{n \times K}{n by K} dimensional (design) matrix.
#' @param M_bdiag A\eqn{n \times Jp}{n by J*p} dimensional block-diagonal design matrix.
#' Each j-th block (\eqn{j = 1,\dots,J}{j = 1,..,J}) is a \eqn{n_j \times p}{n_j by p} dimensional design matrix for the j-th clone.
#' @param y n-dimensional vector of the time-adjacent cellular increments
#' @param VCNs A n-dimensional vector including values of the vector copy number corresponding to the cell counts of y.
#' @param nObs A K-dimensional vector including the frequencies of each clone k (\eqn{k = 1,\dots,K}{k = 1,..,K}).
#' @param trace Non-negative integer. If positive, tracing information on the progress of the optimization is produced.
#' This parameter is also passed to the optim() function.
#' Higher values may produce more tracing information: for method "L-BFGS-B" there are six levels of tracing.
#' (To understand exactly what these do see the source code: higher levels give more detail.)
#' @return A vector of statistics associated to the fitted base and random-effects models:
#' \itemize{
#'  \item{"nPar"}{number of parameters of the base(null) model}
#'  \item{"cll"}{value of the conditional log-likelihood, in this case just the log-likelihood}
#'  \item{"mll"}{value of the marginal log-likelihood, in this case just the log-likelihood}
#'  \item{"cAIC"}{conditional Akaike Information Criterion (cAIC), in this case simply the AIC.}
#'  \item{"mAIC"}{marginal Akaike Information Criterion (mAIC), in this case simply the AIC.}
#'  \item{"Chi2"}{value of the \eqn{\chi^2}{Chi-squared} statistic \eqn{(y - M\theta)'S^{-1}(y - M\theta)}{(y - Mθ)'S^-1(y - Mθ)}.}
#'  \item{"p-value"}{p-value of the \eqn{\chi^2}{Chi-squared} test for the null hypothesis that Chi2
#'  follows a \eqn{\chi^2}{Chi-squared} distribution with n - nPar degrees of freedom.}
#'  \item{"KLdiv"}{Kullback-Leibler divergence of the random-effects model from the null model.}
#'  \item{"KLdiv/N"}{Rescaled Kullback-Leibler divergence of the random-effects model from the null model.}
#'  \item{"BhattDist_nullCond"}{Bhattacharyya distance between the random-effects model and the null model.}
#'  \item{"BhattDist_nullCond/N"}{Rescaled Bhattacharyya distance between the random-effects model and the null model.}
#' }
#' @keywords internal
rndEffModelStats <- function(theta_null,
                             theta_rndEff,
                             V,
                             M,
                             M_bdiag,
                             y,
                             VCNs,
                             nObs,
                             trace = TRUE){

  K <- length(unique(rownames(M)))
  p <- ncol(M)
  N <- nrow(M)

  d2 <- (theta_rndEff[(p + 1): (2*p)])
  D2 <- kronecker(X = Diagonal(n = K, x = 1), Y = Diagonal(n = p, x = d2)) ## covariance matrix of the random effects u, using the free parameter
  q <- ncol(M_bdiag)

  theta_res <- theta_rndEff
  VEuy_res <- VEuy(theta_res, M, M_bdiag, y, V, VCNs, nObs)
  euy_res <- VEuy_res[["euy"]] ## conditional expectaction \eqn{E[u \vert y]}{E[u|y]} of the random (latent) effects u given the observed data y
  vuy_res <- VEuy_res[["vuy"]] ## conditional variance \eqn{V[u \vert y]}{V[u|y]} of the random (latent) effects u given the observed data y

  thetaRndEffRes <- c(theta_res[1:p], as.numeric(euy_res), tail(theta_res,1))
  if(trace){cat("computing AIC...\n")}
  cAIC <- condAIC(X = M, Z = M_bdiag, y = y, theta = thetaRndEffRes, Delta = D2, V, VCNs = VCNs, nObs = nObs, trace = trace)

  b_re <- theta_res[1:p]
  W_re <- get.W(M, V, VCNs, theta = c(thetaRndEffRes[1:p],tail(thetaRndEffRes,1)), nObs)
  S_re <- W_re
  diag(S_re) <- diag(S_re) + tail(thetaRndEffRes,1)
  Si_re <- chol2inv(chol(S_re))


  if(trace){cat("computing Chi2...\n")}
  chi2Stat <- as.numeric(t(y - M_bdiag %*% thetaRndEffRes[(p+1):(p+q)])%*%Si_re%*%(y - M_bdiag %*% thetaRndEffRes[(p+1):(p+q)]))
  chi2pVal <- 1 - pchisq(q = chi2Stat, df = N - (cAIC$rho + 1))

  b_null <- head(theta_null,p)
  s2_null <- tail(theta_null,1)
  muNull <- M%*%b_null

  W_null <- get.W(M, V, VCNs, theta = c(b_null, s2_null), nObs)
  S_null <- W_null
  diag(S_null) <- diag(S_null) + s2_null
  Si_null <- chol2inv(chol(S_null))

  mu_condModel <- M%*%b_re + M_bdiag %*% thetaRndEffRes[(p+1):(p+q)]
  if(trace){cat("computing KLdiv...\n")}
  klDiv_nullCond <- KLDiv(mu1 = muNull, S1 = S_null, mu2 = mu_condModel, S2 = S_re)

  if(trace){cat("computing BhattDist...\n")}
  BhattDist_nullCond <- BhattDist(mu1 = muNull, S1 = S_null, mu2 = mu_condModel, S2 = S_re)

  stats <- c(cAIC$rho + 1, NA, NA, NA, cAIC$cAIC, chi2Stat, chi2pVal, klDiv_nullCond, klDiv_nullCond/N, BhattDist_nullCond, BhattDist_nullCond/N)
  names(stats) <- c("nPar","cll","mll","cAIC","mAIC","Chi2","p-value","KLdiv","KLdiv/N", "BhattDist_nullCond", "BhattDist_nullCond/N")

  return(stats)
}

#' Conditional AIC (cAIC)
#'
#' Conditional AIC (cAIC) of the conditional log-likelihood \eqn{l(y \vert u)}{l(y|u)} of y given the random effects u
#' @param X A \eqn{n \times K}{n by K} dimensional (design) matrix.
#' @param Z A\eqn{n \times Jp}{n by J*p} dimensional block-diagonal design matrix.
#' Each j-th block (\eqn{j = 1,\dots,J}{j = 1,..,J}) is a \eqn{n_j \times p}{n_j by p} dimensional design matrix for the j-th clone.
#' @param y n-dimensional vector of the time-adjacent cellular increments
#' @param theta p-dimensional vector parameter.
#' @param Delta covariance matrix of the random effects u
#' @param V A \eqn{p \times K}{p by K} dimensional net-effect matrix.
#' @param VCNs A n-dimensional vector including values of the vector copy number corresponding to the cell counts of y.
#' @param nObs A K-dimensional vector including the frequencies of each clone k (\eqn{k = 1,\dots,K}{k = 1,..,K}).
#' @param trace Non-negative integer. If positive, tracing information on the progress of the optimization is produced.
#' This parameter is also passed to the optim() function.
#' Higher values may produce more tracing information: for method "L-BFGS-B" there are six levels of tracing.
#' (To understand exactly what these do see the source code: higher levels give more detail.)
#' @return Conditional AIC (cAIC) of the conditional log-likelihood \eqn{l(y \vert u)}{l(y|u)} of y given the random effects u.
#' @keywords internal
condAIC <- function(X, Z, y, theta, Delta, V, VCNs, nObs, trace = TRUE){
  p <- ncol(X)
  q <- ncol(Z)
  b <- theta[1:p]
  s2 <- tail(theta, 1)

  W <- get.W(X, V, VCNs, theta = c(b,s2), nObs)
  S <- W
  diag(S) <- diag(S) + s2
  S <- (S + t(S))/2
  Si <- chol2inv(chol(S))

  n <- nrow(X)

  D_fit <- Delta
  D_fiti <- D_fit
  diag(D_fiti) <- diag(D_fiti)^-1

  if(trace){cat("ny_u...\n")}
  n2ll <- 2*(ny_u(theta, Z, X, y, V, VCNs, nObs = nObs)-n/2*log(2*pi))

  if(trace){cat("computing A and B...\n")}
  A11 <- B11 <- t(X)%*%Si%*%X
  A12 <- B12 <- ((t(X)%*%Si)%*%Z)
  A21 <- B21 <- ((t(Z)%*%Si)%*%X)
  B22 <- (t(Z)%*%Si)%*%Z
  A22 <- B22 + D_fiti


  if(trace){cat("computing A^-1...\n")}
  A22 <- (A22 + t(A22))/2
  A22inv <- chol2inv(chol(A22))
  A22inv_X_A21 <- A22inv%*%A21
  Ei <- A11 - A12%*%(A22inv_X_A21)
  Ei <- (Ei + t(Ei))/2
  E <- solve(Ei)
  A11i <- E
  A12i <- -E%*%(A12%*%A22inv)
  A21i <- -A22inv%*%(A21%*%E)

  if(trace){cat("computing rho...\n")}
  rho <- sum(A11i*B11) + 2*sum(t(A12i)*B21)
  rho <- rho + sum(A22inv*B22) - sum(t(A22inv_X_A21)*(A12i%*%B22))

  cAIC <- n2ll + 2*(rho + 1)

  AICres <- list()
  AICres$cAIC <- cAIC
  AICres$rho <- rho

  return(AICres)
}

#' Conditional negative log-likelihood
#'
#' Conditional negative log-likelihood \eqn{-l(y \vert u)}{-l(y|u)} of y given the random effects u
#' @param theta p-dimensional vector parameter.
#' @param Z A\eqn{n \times Jp}{n by J*p} dimensional block-diagonal design matrix.
#' @param X A \eqn{n \times K}{n by K} dimensional (design) matrix.
#' Each j-th block (\eqn{j = 1,\dots,J}{j = 1,..,J}) is a \eqn{n_j \times p}{n_j by p} dimensional design matrix for the j-th clone.
#' @param y n-dimensional vector of the time-adjacent cellular increments
#' @param V A \eqn{p \times K}{p by K} dimensional net-effect matrix.
#' @param VCNs A n-dimensional vector including values of the vector copy number corresponding to the cell counts of y.
#' @param nObs A K-dimensional vector including the frequencies of each clone k (\eqn{k = 1,\dots,K}{k = 1,..,K}).
#' @return Conditional negative log-likelihood \eqn{-l(y \vert u)}{-l(y|u)} of y given the random effects u.
#' @keywords internal
ny_u <- function(theta, Z, X, y, V, VCNs, nObs){
  p <- ncol(V)
  q <- ncol(Z)
  b <- theta[1:p]
  u <- theta[(p+1):(p+q)]

  s2 <- tail(theta, 1)

  W <- get.W(X, V, VCNs, theta = c(b,s2), nObs)
  S <- W
  diag(S) <- diag(S) + s2
  Si <- chol2inv(chol(S))

  n <- nrow(Z)

  nll <- as.numeric(-.5*ldet(S) -.5*t(y - Z%*%u)%*%Si%*%(y - Z%*%u))
  return(-nll)
}

#' Kullback-Leibler divergence
#'
#' Kullback-Leibler divergence \eqn{KL(P \Vert Q)}{KL(P || Q)} of Q from P for multivariate normal distributions
#' @param mu1 mean vector of P.
#' @param S1 covariance matrix of P.
#' @param mu2 mean vector of Q.
#' @param S2 covariance matrix of Q.
#' @return Kullback-Leibler divergence \eqn{KL(P \Vert Q)}{KL(P || Q)} of Q from P.
#' @keywords internal
KLDiv <- function(mu1, S1, mu2, S2){
  d <- length(mu1)
  klDiv <- as.numeric(.5*(ldet(S2) - ldet(S1) - d + tr(solve(S2)*S1) + t(mu2 - mu1)%*%(solve(S2)%*%(mu2 - mu1))))
  return(klDiv)
}

#' Bhattacharyya distance
#'
#' Bhattacharyya distance between P and Q for multivariate normal distributions
#' @param mu1 mean vector of P.
#' @param S1 covariance matrix of P.
#' @param mu2 mean vector of Q.
#' @param S2 covariance matrix of Q.
#' @return Bhattacharyya distance between P and Q.
#' @keywords internal
BhattDist <- function(mu1, S1, mu2, S2){
  d <- length(mu1)
  S <- (S1 + S2)/2
  res <- as.numeric(1/8*t(mu1 - mu2)%*%S%*%(mu1 - mu2) +.5*ldet(S) -1/4*ldet(S1) -1/4*ldet(S2))
  return(res)
}

#' Fit the base (null) model and get the results
#'
#' This function builds the design matrix of the null model and returns the fitted values and the corresponding statistics.
#' @param Y A 3-dimensional array whose dimensions are the time, the cell type and the clone respectively.
#' @param rct.lst list of biochemical reactions.
#' A differentiation move from cell type "A" to cell type "B" must be coded as "A->B"
#' Duplication of cell "A" must be coded as "A->1"
#' Death of cell "A" must be coded as "A->0"
#' @param maxit maximum number of iterations for the optimization step.
#' This argument is passed to optim() function. Details on "maxit" can be found in "optim()" documentation page.
#' @param factr controls the convergence of the "L-BFGS-B" method.
#' Convergence occurs when the reduction in the objective is within this factor of the machine tolerance.
#' Default is 1e7, that is a tolerance of about 1e-8.
#' This argument is passed to optim() function.
#' @param pgtol helps control the convergence of the "L-BFGS-B" method.
#' It is a tolerance on the projected gradient in the current search direction.
#' This defaults to zero, when the check is suppressed.
#' This argument is passed to optim() function.
#' @param lmm is an integer giving the number of BFGS updates retained in the "L-BFGS-B" method, It defaults to 5.
#' This argument is passed to optim() function.
#' @param trace Non-negative integer. If positive, tracing information on the progress of the optimization is produced.
#' This parameter is also passed to the optim() function.
#' Higher values may produce more tracing information: for method "L-BFGS-B" there are six levels of tracing.
#' (To understand exactly what these do see the source code: higher levels give more detail.)
#' @return A 3-length list. First element is the output returned by "optim()" function (see "optim()" documentation for details).
#' Second element is a vector of statistics associated to the fitted null model:
#' \itemize{
#'  \item{"nPar"}{number of parameters of the base(null) model}
#'  \item{"cll"}{value of the conditional log-likelihood, in this case just the log-likelihood}
#'  \item{"mll"}{value of the marginal log-likelihood, in this case just the log-likelihood}
#'  \item{"cAIC"}{conditional Akaike Information Criterion (cAIC), in this case simply the AIC.}
#'  \item{"mAIC"}{marginal Akaike Information Criterion (mAIC), in this case simply the AIC.}
#'  \item{"Chi2"}{value of the \eqn{\chi^2}{Chi-squared} statistic \eqn{(y - M\theta)'S^{-1}(y - M\theta)}{(y - Mθ)'S^-1(y - Mθ)}.}
#'  \item{"p-value"}{p-value of the \eqn{\chi^2}{Chi-squared} test for the null hypothesis that Chi2
#'  follows a \eqn{\chi^2}{Chi-squared} distribution with n - nPar degrees of freedom.}
#' }
#' The third element, called "design", is a list including:
#' \itemize{
#'  \item{"M"}{A \eqn{n \times K}{n by K} dimensional (design) matrix.}
#'  \item{"V"}{A \eqn{p \times K}{p by K} dimensional net-effect matrix.}
#' }
#' @examples
#' \donttest{
#' library(RestoreNet)
#' library(ggplot2)
#' library(scatterpie)
#' rcts <- c("A->1", "B->1", "C->1", "D->1",
#'           "A->0", "B->0", "C->0", "D->0",
#'           "A->B", "A->C", "C->D") ## set of reactions
#' nC <- 3 ## number of clones
#' S <- 100 ## trajectory length
#' tau <- 1 ## for tau-leaping algorithm
#' u_1 <- c(.2, .15, .17, .09*5,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' u_2 <- c(.2, .15, .17, .09,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' u_3 <- c(.2, .15, .17*3, .09,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' theta_allcls <- cbind(u_1, u_2, u_3) ## clone-specific parameters
#' rownames(theta_allcls) <- colnames(V)
#' s20 <- 1 ## additional noise
#' Y <- array(data = NA,
#'            dim = c(S + 1, nrow(V), nC),
#'            dimnames = list(seq(from = 0, to = S*tau, by = tau),
#'                            rownames(V),
#'                            1:nC)) ## empty array to store simulations
#' Y0 <- c(100,0,0,0) ## initial state
#' names(Y0) <- rownames(V)
#' for (cl in 1:nC) { ## loop over clones
#'   Y[,,cl] <- get.sim.tl(Yt = Y0,
#'                         theta = theta_allcls[,cl],
#'                         S = S,
#'                         s2 = s20,
#'                         tau = tau,
#'                         rct.lst = rcts)
#' }
#' null.res <- fit.null(Y = Y, rct.lst = rcts) ## null model fitting
#' }
##' @export
fit.null <- function(Y,
                         rct.lst,
                         maxit = 10000,
                         factr = 1e7,
                         pgtol = 1e-8,
                         lmm = 100,
                         trace = TRUE){

  # pattern <- "[A-Z]+[0-9]*->([A-Z]+[0-9]*|[0-1]+)"
  pattern <- "[A-Z0-9]{1,}->[A-Z0-9]{0,}[0-1]{0,}"

  if(sum(grepl(pattern,
               rct.lst,
               ignore.case = FALSE)) < length(rct.lst)){
    cat(paste("Unvalid reactions provided in 'rct.lst': ",
              paste(rct.lst[which(!grepl(pattern,
                                         c(rct.lst),
                                         ignore.case = FALSE))], collapse = ", "), sep = ""))
    return()
  }

  # ct.lst <- unique(setdiff(c(sapply(rct.lst, function(r){
  #   as.vector(unlist(strsplit(r, split = "->", fixed = T)))
  # }, simplify = "array")), c("0", "1")))
  # ct.lst <- setdiff(as.vector(unlist(c(sapply(rct.lst, function(r){
  #   as.vector(unlist(strsplit(r, split = "->", fixed = T)))
  # }, simplify = "array")))), c("0","1"))
  ct.lst <- setdiff(as.vector(unlist(c(sapply(rct.lst, function(r){
    as.vector(unlist(str_split(r, "->")))
  }, simplify = "array")))), c("0","1"))

  if(sum(!(ct.lst %in% colnames(Y))) > 0){
    cat(paste("Cell types not present in 'Y': ",
              paste(as.vector(unlist(lapply(ct.lst[which(!(ct.lst %in% colnames(Y)))],
                                            function(l){paste("'", l, "'", sep = "")}))), collapse = ", "),
              "\nPlease provide compatible cell types in 'Y' and 'rct.lst'.",
              sep = ""))
    return()
  }

  compile.h(rct.lst = rct.lst, envir = environment())
  V <- get.V(ct.lst = ct.lst, rct.lst = rct.lst)

  Y <- Y[,rownames(V),]
  Y_rbind <- Reduce(rbind, lapply(1:dim(Y)[3], function(cl){return(head(Y[,,cl],-1))}))
  dT <- rep(rep(diff(as.numeric(rownames(Y))), each = ncol(Y)), times = dim(Y)[3])

  if(trace){cat("creating design matrix...")}
  M <- get.M(y = Y_rbind, V = V, dT = dT, get.h = get.h)
  dx <- unlist(lapply(1:dim(Y)[3], function(cl){return(get.dx(Y[,,cl]))}))
  VCNs <- rep(1, length(dx))
  clones <- as.vector(unlist(dimnames(Y)[3]))
  y <- dx
  y_all <- y
  M_all <- M
  clones_all <- clones <- rep(clones, each = (nrow(Y) - 1)*ncol(Y))
  VCNs_all <- VCNs
  yfull <- y
  Vfull <- V
  yM <- cbind(y,M)

  yM <- yM[rownames(M) == "T" | names(y) == "T",]
  clones <- clones[rownames(M) == "T" | names(y) == "T"]
  VCNs <- VCNs[rownames(M) == "T" | names(y) == "T"]
  dT <- dT[rownames(M) == "T" | names(y) == "T"]

  y <- dx <- yM[,1]
  Mfull <- M <- yM[,-1]

  rownames(M) <- names(y) <- names(VCNs) <- rownames(yM) <- clones
  nrow(M) == length(dx) & length(dx) == length(VCNs) & length(VCNs) == length(dT)
  cat(" DONE\n")
  nObs <- apply(apply(Y != 0, 3, rowSums) != 0, 2, sum)/max(apply(apply(Y != 0, 3, rowSums) != 0, 2, sum))

  if(trace){cat(paste("Fitting null model...\n", sep = ""))}
  p <- ncol(V)
  epsLB <- 1e-7

  resNullFullV <- nullModelFitting(theta_start = rep(1, p + 1),
                                   M = M,
                                   y = y,
                                   V = V,
                                   psiLB = rep(epsLB, p + 1),
                                   psiUB = rep(+Inf, p + 1),
                                   maxit = maxit,
                                   factr = factr,
                                   pgtol = pgtol,
                                   lmm = lmm,
                                   VCNs = VCNs,
                                   nObs = nObs,
                                   trace = trace);
  if(trace){cat(" DONE\n")}

  return(resNullFullV)
}

#' Fit the random-effects model and get the results
#'
#' This function builds the design matrix of the random-effects model and returns the fitted values and the corresponding statistics.
#' @param theta_0 A p-dimensional vector parameter as the initial guess for the inference.
#' @param Y A 3-dimensional array whose dimensions are the time, the cell type and the clone respectively.
#' @param rct.lst list of biochemical reactions.
#' A differentiation move from cell type "A" to cell type "B" must be coded as "A->B"
#' Duplication of cell "A" must be coded as "A->1"
#' Death of cell "A" must be coded as "A->0"
#' @param maxit maximum number of iterations for the optimization step.
#' This argument is passed to optim() function. Details on "maxit" can be found in "optim()" documentation page.
#' @param factr controls the convergence of the "L-BFGS-B" method.
#' Convergence occurs when the reduction in the objective is within this factor of the machine tolerance.
#' Default is 1e7, that is a tolerance of about 1e-8.
#' This argument is passed to optim() function.
#' @param pgtol helps control the convergence of the "L-BFGS-B" method.
#' It is a tolerance on the projected gradient in the current search direction.
#' This defaults to zero, when the check is suppressed.
#' This argument is passed to optim() function.
#' @param lmm is an integer giving the number of BFGS updates retained in the "L-BFGS-B" method, It defaults to 5.
#' This argument is passed to optim() function.
#' @param maxemit maximum number of iterations for the expectation-maximization algorithm.
#' @param eps relative error for the value x and the objective function f(x) that has to be optimized
#' in the expectation-maximization algorithm.
#' @param trace Non-negative integer. If positive, tracing information on the progress of the optimization is produced.
#' This parameter is also passed to the optim() function.
#' Higher values may produce more tracing information: for method "L-BFGS-B" there are six levels of tracing.
#' (To understand exactly what these do see the source code: higher levels give more detail.)
#' @return A 3-length list. First element is the output returned by "optim()" function (see "optim()" documentation for details)
#' along with the conditional expectation \eqn{E[u \vert y]}{E[u|y]} and variance \eqn{V[u \vert y]}{V[u|y]}
#' of the latent states u given the observed states y from the last step of the expectation-maximization algorithm.
#' Second element is a vector of statistics associated to the fitted random-effects model:
#' \itemize{
#'  \item{"nPar"}{number of parameters of the base(null) model}
#'  \item{"cll"}{value of the conditional log-likelihood, in this case just the log-likelihood}
#'  \item{"mll"}{value of the marginal log-likelihood, in this case just the log-likelihood}
#'  \item{"cAIC"}{conditional Akaike Information Criterion (cAIC), in this case simply the AIC.}
#'  \item{"mAIC"}{marginal Akaike Information Criterion (mAIC), in this case simply the AIC.}
#'  \item{"Chi2"}{value of the \eqn{\chi^2}{Chi-squared} statistic \eqn{(y - M\theta)'S^{-1}(y - M\theta)}{(y - Mθ)'S^-1(y - Mθ)}.}
#'  \item{"p-value"}{p-value of the \eqn{\chi^2}{Chi-squared} test for the null hypothesis that Chi2
#'  follows a \eqn{\chi^2}{Chi-squared} distribution with n - nPar degrees of freedom.}
#' \item{"KLdiv"}{Kullback-Leibler divergence of the random-effects model from the null model.}
#' \item{"KLdiv/N"}{Rescaled Kullback-Leibler divergence of the random-effects model from the null model.}
#' \item{"BhattDist_nullCond"}{Bhattacharyya distance between the random-effects model and the null model.}
#' \item{"BhattDist_nullCond/N"}{Rescaled Bhattacharyya distance between the random-effects model and the null model.}
#' }
#' The third element, called "design", is a list including:
#' \itemize{
#'  \item{"M"}{A \eqn{n \times K}{n by K} dimensional (design) matrix.}
#'  \item{"M_bdiag"}{ A\eqn{n \times Jp}{n by J*p} dimensional block-diagonal design matrix.}
#'  \item{"V"}{A \eqn{p \times K}{p by K} dimensional net-effect matrix.}
#' }
#' @examples
#' \donttest{
#' library(RestoreNet)
#' library(ggplot2)
#' library(scatterpie)
#' rcts <- c("A->1", "B->1", "C->1", "D->1",
#'           "A->0", "B->0", "C->0", "D->0",
#'           "A->B", "A->C", "C->D") ## set of reactions
#' nC <- 3 ## number of clones
#' S <- 100 ## trajectory length
#' tau <- 1 ## for tau-leaping algorithm
#' u_1 <- c(.2, .15, .17, .09*5,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' u_2 <- c(.2, .15, .17, .09,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' u_3 <- c(.2, .15, .17*3, .09,
#'          .001, .007, .004, .002,
#'          .13, .15, .08)
#' theta_allcls <- cbind(u_1, u_2, u_3) ## clone-specific parameters
#' rownames(theta_allcls) <- colnames(V)
#' s20 <- 1 ## additional noise
#' Y <- array(data = NA,
#'            dim = c(S + 1, nrow(V), nC),
#'            dimnames = list(seq(from = 0, to = S*tau, by = tau),
#'                            rownames(V),
#'                            1:nC)) ## empty array to store simulations
#' Y0 <- c(100,0,0,0) ## initial state
#' names(Y0) <- rownames(V)
#' for (cl in 1:nC) { ## loop over clones
#'   Y[,,cl] <- get.sim.tl(Yt = Y0,
#'                         theta = theta_allcls[,cl],
#'                         S = S,
#'                         s2 = s20,
#'                         tau = tau,
#'                         rct.lst = rcts)
#' }
#' null.res <- fit.null(Y = Y, rct.lst = rcts) ## null model fitting
#'
#' re.res <- fit.re(theta_0 = null.res$fit$par,
#'                      Y = Y,
#'                      rct.lst = rcts,
#'                      maxemit = 100) ## random-effects model fitting
#' }
##' @export
fit.re <- function(theta_0,
                       Y,
                       rct.lst,
                       maxit = 10000,
                       factr = 1e7,
                       pgtol = 1e-8,
                       lmm = 100,
                       maxemit = 100,
                       eps = 1e-5,
                       trace = TRUE){

  # ct.lst <- unique(setdiff(c(sapply(rct.lst, function(r){
  #   as.vector(unlist(strsplit(r, split = "->", fixed = T)))
  # }, simplify = "array")), c("0", "1")))

  # pattern <- "[A-Z]+[0-9]*->([A-Z]+[0-9]*|[0-1]+)"
  pattern <- "[A-Z0-9]{1,}->[A-Z0-9]{0,}[0-1]{0,}"

  if(sum(grepl(pattern,
               rct.lst,
               ignore.case = FALSE)) < length(rct.lst)){
    cat(paste("Unvalid reactions provided in 'rct.lst': ",
              paste(rct.lst[which(!grepl(pattern,
                                         c(rct.lst),
                                         ignore.case = FALSE))], collapse = ", "), sep = ""))
    return()
  }

  # ct.lst <- unique(setdiff(c(sapply(rct.lst, function(r){
  #   as.vector(unlist(strsplit(r, split = "->", fixed = T)))
  # }, simplify = "array")), c("0", "1")))
  # ct.lst <- setdiff(as.vector(unlist(c(sapply(rct.lst, function(r){
  #   as.vector(unlist(strsplit(r, split = "->", fixed = T)))
  # }, simplify = "array")))), c("0","1"))
  ct.lst <- setdiff(as.vector(unlist(c(sapply(rct.lst, function(r){
    as.vector(unlist(str_split(r, "->")))
  }, simplify = "array")))), c("0","1"))

  if(sum(!(ct.lst %in% colnames(Y))) > 0){
    cat(paste("Cell types not present in 'Y': ",
              paste(as.vector(unlist(lapply(ct.lst[which(!(ct.lst %in% colnames(Y)))],
                                            function(l){paste("'", l, "'", sep = "")}))), collapse = ", "),
              "\nPlease provide compatible cell types in 'Y' and 'rct.lst'.",
              sep = ""))
    return()
  }

  compile.h(rct.lst = rct.lst, envir = environment())
  V <- get.V(ct.lst = ct.lst, rct.lst = rct.lst)

  Y <- Y[,rownames(V),]
  Y_rbind <- Reduce(rbind, lapply(1:dim(Y)[3], function(cl){return(head(Y[,,cl],-1))}))
  dT <- rep(rep(diff(as.numeric(rownames(Y))), each = ncol(Y)), times = dim(Y)[3])

  if(trace){cat("creating design matrix...")}
  M <- get.M(y = Y_rbind, V = V, dT = dT, get.h = get.h)
  dx <- unlist(lapply(1:dim(Y)[3], function(cl){return(get.dx(Y[,,cl]))}))
  VCNs <- rep(1, length(dx))
  clones <- as.vector(unlist(dimnames(Y)[3]))
  y <- dx
  y_all <- y
  M_all <- M
  clones_all <- clones <- rep(clones, each = (nrow(Y) - 1)*ncol(Y))
  VCNs_all <- VCNs
  yfull <- y
  Vfull <- V
  yM <- cbind(y,M)

  yM <- yM[rownames(M) == "T" | names(y) == "T",]
  clones <- clones[rownames(M) == "T" | names(y) == "T"]
  VCNs <- VCNs[rownames(M) == "T" | names(y) == "T"]
  dT <- dT[rownames(M) == "T" | names(y) == "T"]

  y <- dx <- yM[,1]
  Mfull <- M <- yM[,-1]

  rownames(M) <- names(y) <- names(VCNs) <- rownames(yM) <- clones
  nrow(M) == length(dx) & length(dx) == length(VCNs) & length(VCNs) == length(dT)
  cat(" DONE\n")
  nObs <- apply(apply(Y != 0, 3, rowSums) != 0, 2, sum)/max(apply(apply(Y != 0, 3, rowSums) != 0, 2, sum))

  p <- ncol(V)
  epsLB <- 1e-7

  if(trace){cat(paste("Fitting RE model with full V...\n", sep = ""))}
  p <- ncol(V)
  re.fit <- rndEffModelFitting(theta_0 = theta_0,
                               V = V,
                               M = M,
                               M_bdiag = bdiag(lapply(unique(rownames(M)), function(cl){return(M[which(rownames(M) == cl),])})),
                               y = y,
                               VCNs = VCNs,
                               nObs = nObs,
                               maxit = maxit,
                               maxemit = maxemit,
                               eps = eps,
                               thetaLB = rep(epsLB, 2*p + 1),
                               thetaUB = rep(+Inf, 2*p + 1),
                               factr = factr,
                               pgtol = pgtol,
                               lmm = lmm,
                               trace = trace)
  if(trace){cat(" DONE\n")}

  if(trace){cat(paste("Computing statistics...\n", sep = ""))}
  re.stats <- rndEffModelStats(theta_null = theta_0,
                               theta_rndEff = re.fit$par,
                               V = V,
                               M = M,
                               M_bdiag = bdiag(lapply(unique(rownames(M)), function(cl){return(M[which(rownames(M) == cl),])})),
                               y = y,
                               VCNs = VCNs,
                               nObs = nObs,
                               trace = trace)

  res <- list()
  res$fit <- re.fit
  res$stats <- re.stats

  design <- list()
  design$M <- M
  design$M_bdiag <- bdiag(lapply(unique(rownames(M)), function(cl){return(M[which(rownames(M) == cl),])}))
  design$V <- V
  res$design <- design

  return(res)
}
