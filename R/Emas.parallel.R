#' @title Epigenome-Wide Mediation Analysis Study: Parallel multiple mediation model
#'
#' @description This function can perform the parallel multiple mediation model after the Epigenome-Wide Mediation Analysis Study (EMAS).
#'
#' @details This function can perform the parallel multiple mediation model after the Epigenome-Wide Mediation Analysis Study (EMAS) to further explore the potential parallel mediating CpG sites of exposure variables affecting outcome variables.
#'
#' @param data A data.frame included id, x, y, x.cov, y.cov, m.cov.
#'
#' @param M.matrix A matrix with the CpG information screened from EMAS., maybe a M-value matrix or a beta value matrix.
#'
#' @param id Variable name of the id.
#'
#' @param x,y Variable name of exposure(x) and outcome(y).
#'
#' @param x.cov Variable names of covariates related to exposure(x).
#'
#' @param y.cov Variable names of covariates related to outcome(y).
#'
#' @param m.cov Variable names of covariates related to mediator(m).
#'
#' @param m.cor A logical value. If 'TRUE', the mediators in the parallel multiple mediation model are set to correlate with each other.
#'
#' @param boot A logical value or a numeric value. If a numeric value, the number for bootstrap.
#'
#' @param lavaan A logical value. If 'TRUE', a lavaan object will be given.
#' 
#' @param ... Other arguments passed to \code{\link{sem}} from \code{\link{lavaan}} package.
#'
#' @export
#'
#' @return \code{Emas.parallel} returns a data.frame with the average mediation effects(AME), 
#'    average direct effects(ADE), and total effects(Tot). If \code{lavaan} is 'TRUE', a lavaan object will be given.
#' 
#' @importFrom "lavaan" "sem" "summary"
#' @importFrom "stats" "lm" "na.omit" "residuals" "pnorm"
#' @importFrom "utils" "capture.output"
#' 
#' @author Xiuquan Nie, niexiuquan1995@foxmail.com
#' 
#' @references
#' Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling.
#' \emph{Journal of Statistical Software}, 48(2), 1–36. \doi{doi:10.18637/jss.v048.i02}.
#' 
#' @examples
#' \donttest{data(data.m)
#' data(Mvalue)
#' EP.result <- Emas.parallel(data.m, Mvalue, 
#'                            id = "ID", x = "x", y = "y", 
#'                            x.cov = c("age", "gender"), 
#'                            y.cov = c("age", "gender"), 
#'                            m.cov = c("age", "gender", "CD8T", "CD4T"), 
#'                            m.cor = TRUE, boot = FALSE, lavaan = FALSE)}
Emas.parallel <- function(data, M.matrix, 
                          id = "", x = "", y = "", 
                          x.cov = c(), y.cov = c(), m.cov = c(),
                          m.cor = TRUE, boot = FALSE, 
                          lavaan = FALSE, ...) {
  cn <- c(id, x, y, x.cov, y.cov, m.cov)
  cn <- cn[!duplicated(cn)]
  data.v <- na.omit(data[cn])
  colnames(data.v) <- c("ID", cn[-1])
  data.v.ID <- subset(data.v, select = c(ID))
  
  zzM <- t(M.matrix)
  Mv <- as.data.frame(zzM)
  ID <- rownames(Mv)
  MV <- cbind(ID, Mv)
  MV.ID <- subset(MV, select = c(ID))
  MV.ID$de <- 1
  
  ID.z <- merge(data.v.ID, MV.ID, by = "ID", all = F)
  ID.z <- ID.z[order(ID.z$ID),]
  ID.z$de <- NULL
  
  data.z <- merge(ID.z, data.v, by = "ID", all.x = T)
  data.z <- data.z[order(data.z$ID),]
  MV.z <- merge(ID.z, MV, by = "ID", all.x = T)
  MV.z <- MV.z[order(MV.z$ID),]
  snum <- dim(data.z)[1]
  MV.z$ID <- NULL
  
  message(paste("A total of", snum, "samples were included in the parallel multiple mediation model.", sep = " "))
  
  x.c <- x.cov[1]
  for (i in x.cov[-1]) {
    x.c <- paste(x.c, i, sep = " + ")
  }
  formula.x <- paste(x, " ~ ", x.c, sep = "")
  modx <- residuals(lm(formula.x, data = data.z))
  
  y.c <- y.cov[1]
  for (i in y.cov[-1]) {
    y.c <- paste(y.c, i, sep = " + ")
  }
  formula.y <- paste(y, " ~ ", y.c, sep = "")
  mody <- residuals(lm(formula.y, data = data.z))
  
  data.xy <- cbind(data.z, modx, mody)
  data.xy <- subset(data.xy, select = c(ID, modx, mody))
  
  m.c <- m.cov[1]
  for (i in m.cov[-1]) {
    m.c <- paste(m.c, i, sep = " + ")
  }
  formula.m <- paste("m", " ~ ", m.c, sep = "")
  
  funcan <- function(m){
    data.m <- cbind(data.z, m)
    modm <- residuals(lm(formula.m, data = data.m))
    if(length(modm) == snum){
      data.mm <- cbind(data.m, modm)
      row.names(data.mm) <- data.mm$ID
      modmm <- subset(data.mm, select = c(modm))
      colnames(modmm) <- colnames(m)
      return(modmm)
    }
  }
  cancha <- lapply(MV.z, funcan)
  
  if(class(cancha)[1] == "list"){
    cancha <- cancha[!sapply(cancha, is.null)]
  }
  
  cancha <- data.frame(cancha)
  mnum <- dim(cancha)[2]
  
  message(paste("A total of", mnum, "mediators were included in the parallel multiple mediation model.", sep = " "))
  
  cancha$ID <- row.names(cancha)
  data.c <- merge(data.xy, cancha, by = "ID", all = F)
  row.names(data.c) <- data.c$ID
  cancha$ID <- NULL
  canm <- colnames(cancha)
  
  yloop <- "mody ~ "
  xloop <- ""
  inloop <- ""
  cloop <- "Tot := ADE"
  mloop <- ""
  for (i in c(1:length(canm))) {
    yloop <- paste(yloop, "b", i, " * ", canm[i], " + ", sep = "")
  }
  for (i in c(1:length(canm))) {
    xloop <- paste(xloop, canm[i], " ~ ", "a", i, " * ", "modx\n", sep = "")
  }
  for (i in c(1:length(canm))) {
    inloop <- paste(inloop, "AME.", canm[i], " := ", "a", i, " * ", "b", i, "\n", sep = "")
  }
  for (i in c(1:length(canm))) {
    cloop <- paste(cloop, " + (a", i, " * ", "b", i, ")", sep = "")
  }
  for (i in c(1:(length(canm)-1))) {
    for (j in c((i+1):length(canm))) {
      mloop <- paste(mloop, canm[i], " ~~ ", canm[j], "\n", sep = "")
    }
  }
  if (m.cor == TRUE) {
    mpm <- paste(yloop, "ADE * modx\n", xloop, inloop, cloop, "\n", mloop, sep = "")
  } else {
    mpm <- paste(yloop, "ADE * modx\n", xloop, inloop, cloop, sep = "")
  }
  if (is.logical(boot) == T) {
    if (boot == FALSE) {
      fit <- lavaan::sem(model = mpm, data = data.c, ...)
    } else {
      stop("If boot is TRUE, please enter the bootstrap number.")
    }
  } else {
    fit <- lavaan::sem(model = mpm, data = data.c, se = "bootstrap", bootstrap = boot, ...)
  }
  capture.output(sfit <- lavaan::summary(fit))
  if (m.cor == TRUE) {
    AME <- sfit$PE[c((3*mnum+4+(mnum*(mnum-1)/2)):(4*mnum+3+(mnum*(mnum-1)/2))), c(6:9)]
    ADE.Tot <- sfit$PE[c(mnum+1, (4*mnum+4+(mnum*(mnum-1)/2))), c(6:9)]
    result <- rbind(AME, ADE.Tot)
  } else {
    AME <- sfit$PE[c((3*mnum+4):(4*mnum+3)), c(6:9)]
    ADE.Tot <- sfit$PE[c(mnum+1, (4*mnum+4)), c(6:9)]
    result <- rbind(AME, ADE.Tot)
  }
  if (!is.null(result)) {
    row.names(result) <- c(canm, c("ADE", "Tot"))
  }
  if (lavaan == TRUE) {
    result.l <- list(result = result, lavaan = fit)
    result <- result.l
  }
  return(result)
}
