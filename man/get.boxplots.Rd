% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{get.boxplots}
\alias{get.boxplots}
\title{Clonal boxplots}
\usage{
get.boxplots(re.res)
}
\arguments{
\item{re.res}{output list returned by fit.re().}
}
\description{
Draw clonal boxplots of a random-effects reaction network.
}
\details{
This function generates the boxplots of the conditional expectations
\deqn{w_k = E_{u\vert \Delta Y; \hat{\psi}}[u^k_{\alpha_{l}}] - E_{u\vert \Delta Y; \hat{\psi}}[u^k_{\delta_{l}}]}{Eu|y;ψ[α_l] - Eu|y;ψ[δ_l]},
computed from the estimated parameters \eqn{\hat{\psi}}{ψ} for the clone-specific net-duplication in each cell lineage l (different colors).
The whiskers extend to the data extremes.
}
\examples{
\donttest{
library(RestoreNet)
library(ggplot2)
library(scatterpie)
rcts <- c("A->1", "B->1", "C->1", "D->1",
          "A->0", "B->0", "C->0", "D->0",
          "A->B", "A->C", "C->D") ## set of reactions
nC <- 3 ## number of clones
S <- 100 ## trajectory length
tau <- 1 ## for tau-leaping algorithm
u_1 <- c(.2, .15, .17, .09*5,
         .001, .007, .004, .002,
         .13, .15, .08)
u_2 <- c(.2, .15, .17, .09,
         .001, .007, .004, .002,
         .13, .15, .08)
u_3 <- c(.2, .15, .17*3, .09,
         .001, .007, .004, .002,
         .13, .15, .08)
theta_allcls <- cbind(u_1, u_2, u_3) ## clone-specific parameters
rownames(theta_allcls) <- colnames(V)
s20 <- 1 ## additional noise
Y <- array(data = NA,
           dim = c(S + 1, nrow(V), nC),
           dimnames = list(seq(from = 0, to = S*tau, by = tau),
                           rownames(V),
                           1:nC)) ## empty array to store simulations
Y0 <- c(100,0,0,0) ## initial state
names(Y0) <- rownames(V)
for (cl in 1:nC) { ## loop over clones
  Y[,,cl] <- get.sim.tl(Yt = Y0,
                        theta = theta_allcls[,cl],
                        S = S,
                        s2 = s20,
                        tau = tau,
                        rct.lst = rcts)
}
null.res <- fit.null(Y = Y, rct.lst = rcts) ## null model fitting

re.res <- fit.re(theta_0 = null.res$fit$par,
                     Y = Y,
                     rct.lst = rcts,
                     maxemit = 100) ## random-effects model fitting

get.boxplots(re.res)
}
}
