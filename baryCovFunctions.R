##################################################################################################################
#Functions related to the paper 
#M. Rosario Oliveira, Diogo Pinheiro, & Lina Oliveira. Location and association measures for interval data based on Mallows' distance.
#Corresponding author: rosario.oliveira@tecnico.ulisboa.pt
#Last update: 01/07/2024
##################################################################################################################
#Library

##################################################################################################################
#' Pairs plot for Interval-valued Symbolic data.
#'
#' Adapted from "pairs.panels" (R package "psych") shows a scatter plot of matrices, with bivariate symbolic scatter plots below the main diagonal, variables' names on the main diagonal, and the symbolic covariances or correlations above the main diagonal. Useful for descriptive statistics of symbolic objects characterized by interval variables.
#' 
#' @param C The data matrix (nxp) containing the centres of the intervals. Each row represents an observation and columns represent variables.
#' @param R The data matrix (nxp) containing the ranges of the intervals. Each row represents an observation and columns represent variables.
#' @param delta The real number in [0,1/4] that weights the matrix of the second moments of the ranges when computing the symbolic covariance matrix. Default is delta = 1/12, which corresponds to assuming that the microdata follow a uniform distribution.
#' @param cor.print Should the symbolic correlations be displayed above the main diagonal, instead of symbolic covariances? Default is cor.print = FALSE.
#' @param cor.cex The character expansion factor. Default is cor.cex = 1.
#' @param palette A list of colors to represent the objects. Default is palette = grDevices::rainbow(n).
#' @param B.print Should the symbolic barycenter of each pair of variables be represented by a green rectangle? Default is B.print = FALSE.
#' @param ... Optional graphical parameters of "pairs" can be supplied.
#' 
#' @return A scatter plot matrix (SPLOM) is drawn in the graphic window. The diagonal contains the variables' names, and remaining entries are scatter plots between the variables. If cor.print = True, the entries above the main diagonal are replaced with the symbolic correlation in [Oliveira et al., 2022] and [Oliveira et al., 2024].
#' @examples
#' data(creditCardDataset)
#' pairSymb(C = creditCardDataset$creditCard.c, R = creditCardDataset$creditCard.r,
#'         delta = 1/24, cor.print = TRUE, cor.cex = 1.2, 
#'         palette = c("black", "red", "blue")[creditCardDataset$creditCard.class], B.print = TRUE)
#' 
#' @references Oliveira, M.R., Pinheiro, D., & Oliveira, L. (2024). Location and association measures for interval data based on Mallows' distance.
#' @references Oliveira, M.R., Azeitona, M., Pacheco, A., & Valadas, R. (2022). Association measures for interval variables. Advances in Data Analysis and Classification, 1-30. https://doi.org/10.1007/s11634-021-00445-8
#' @references Serrao, R.G., Oliveira, M.R., & Oliveira, L. (2023). Theoretical derivation of interval principal component analysis. Information Sciences, 621, 227-247. https://doi.org/10.1016/j.ins.2022.11.093

#' @export
pairSymb <- function(C, R, delta = 1/12, 
                     cor.print = FALSE, 
                     cor.cex = 1.0, 
                     palette = grDevices::rainbow(nrow(C)), 
                     B.print = FALSE,...) {
  
  "panel.corSymb" <- function(x, y, digits = 3, palette = palette,...) {
    usr <- graphics::par("usr")
    on.exit(graphics::par(usr))
    graphics::par(usr = c(0, 1, 0, 1))
    
    x.len <- length(x)
    y.len <- length(y)
    x.min <- x[seq(1, x.len/2)]
    x.max <- x[seq(x.len/2 + 1, x.len)]
    y.min <- y[seq(1, y.len/2)]
    y.max <- y[seq(y.len/2 + 1, y.len)]
    
    x.C <- (x.min + x.max)/2
    x.R <- (x.max - x.min)
    y.C <- (y.min + y.max)/2
    y.R <- (y.max - y.min)
    
    auxC <- cbind(x.C, y.C)
    auxR <- cbind(x.R, y.R)
    
    covC <- stats::var(auxC)
    muR <- colMeans(auxR)
    ER2 <- stats::var(auxR) + muR%*%t(muR)
    
    ### Symbolic Covariance matrices
    
    #Symbolic Covariance Matrix according to [Oliveira et al., 2022]
    covDiag <- covC + delta*diag(diag(ER2)) 
    
    #Symbolic Covariance Matrix using the barycentre approach in [Oliveira et al., 2024]
    covB <- covC + delta*stats::var(auxR)
    
    if (cor.print == TRUE) {
      ### Symbolic Correlation matrices 
      auxD <- diag(1/sqrt(diag(covDiag)))
      corDiag <- auxD%*%covDiag%*%auxD
      
      auxB <- diag(1/sqrt(diag(covB)))
      corB <- auxB%*%covB%*%auxB
      
      corDiag.val <- corDiag[1,2]
      corB.val <- corB[1,2]
      
      txtDiag <- format(round(corDiag.val, digits), nsmall = digits)
      text(0.5, 0.625, txtDiag, cex = cor.cex)
      txtB <- format(round(corB.val, digits), nsmall = digits)
      text(0.5, 0.4, txtB, cex = cor.cex)
      
    } else {
      covDiag.val <- covDiag[1,2]
      covB.val <- covB[1,2]
      txtDiag <- format(round(covDiag.val, digits), nsmall = digits)
      text(0.5, 0.625, txtDiag, cex = cor.cex)
      txtB <- format(round(covB.val, digits), nsmall = digits)
      text(0.5, 0.4, txtB, cex = cor.cex)
    }   
  }###END panel.corSymb
  
  "panel.rectSymb" <- function(x, y, n = length(x)/2, B.colour = "green",...) {
    
    x.len <- length(x)
    y.len <- length(y)
    x.min <- x[seq(1, x.len/2)]
    x.max <- x[seq(x.len/2 + 1, x.len)]
    y.min <- y[seq(1, y.len/2)]
    y.max <- y[seq(y.len/2 + 1, y.len)]
    
    xB.min <- mean(x.min)
    xB.max <- mean(x.max)  
    yB.min <- mean(y.min)
    yB.max <- mean(y.max) 
    
    R <- matrix(0, n, 4)
    
    for (i in seq(1,n)) {
      R[i,1] <- x.min[i]
      R[i,2] <- y.min[i]
      R[i,3] <- x.max[i]
      R[i,4] <- y.max[i]
      graphics::rect(R[i,1], R[i,2], R[i,3], R[i,4], border = palette[i], lwd = 1.1)
    }
    
    if (B.print == TRUE) graphics::rect(xB.min, yB.min, xB.max, yB.max, border = B.colour, col = B.colour, density = 60, lwd = 1.1)
  }###END panel.rectSymb
  
  nn <- nrow(C)
  x <- matrix(NA, 2*nn, ncol(C))
  colnames(x) <- colnames(C)
  x[seq(1,nn),] <- as.matrix(C - R/2)      ###lower limits of the intervals
  x[seq(nn+1,2*nn),] <- as.matrix(C + R/2) ###upper limits of the intervals
  
  graphics::pairs(x, upper.panel = panel.corSymb, lower.panel = panel.rectSymb, ...)
}
