#' @title Epigenome-Wide Mediation Analysis Study
#'
#' @description This function can perform the Epigenome-Wide Mediation Analysis Study (EMAS).
#'
#' @details This function can perform the Epigenome-Wide Mediation Analysis Study (EMAS) to explore the potential mediating CpG sites of exposure variables affecting outcome variables within the epigenome-wide.
#'
#' @param data A data.frame included id, x, y, x.cov, y.cov, m.cov.
#'
#' @param M.matrix A matrix with the epigenome-wide CpG information, maybe a M-value matrix or a beta value matrix.
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
#' @param ini.sims Initial number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation.
#'
#' @param boot A logical value. if 'FALSE' a quasi-Bayesian approximation is used for confidence intervals; if 'TRUE' nonparametric bootstrap will be used. Default is 'FALSE'.
#'
#' @param cl.n Number of cores used for parallel computing.
#' 
#' @param ... Other arguments passed to \code{mediated} and \code{boot}.
#'
#' @export
#'
#' @return \code{Emas} returns a data.frame with the average mediation effects(AME), 
#'    average direct effects(ADE), total effects, mediation proportion.
#'    
#'    AMEEst: Point estimates for average mediation effects under the exposure conditions.
#'    
#'    AMElow95, AMEupp95: 95 percentage confidence intervals for average mediation effects.
#'    
#'    AME.P: Two-sided p-values for average mediation effects.
#'    
#'    ADEEst: Point estimates for average direct effect under the exposure conditions.
#'    
#'    ADElow95, ADEupp95: 95 percentage confidence intervals for average direct effects.
#'    
#'    ADE.P: Two-sided p-values for average direct effects.
#'    
#'    TotEst: Point estimate for total effect.
#'    
#'    Totlow95, Totupp95: 95 percentage confidence interval for total effect.
#'    
#'    Tot.P: Two-sided p-values for total effect.
#'    
#'    PropEst: The 'proportions mediated', or the size of the average mediation effects relative to the total effect.
#'    
#'    Proplow95, Propupp95: 95 percentage confidence intervals for the proportions mediated.
#'    
#'    Prop.p: Two-sided p-values for proportions mediated.
#'
#' @importFrom "mediation" "mediate"
#' @importFrom "parallel" "makeCluster" "clusterExport" "parSapply" "stopCluster"
#' @importFrom "stats" "lm" "na.omit" "residuals"
#' 
#' @author Xiuquan Nie, niexiuquan1995@foxmail.com
#' 
#' @references
#' Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L. (2014). mediation: R package for Causal Mediation Analysis.
#' \emph{Journal of Statistical Software}, 59(5), 1–38. \doi{doi:10.18637/jss.v059.i05}.
#' 
#' @examples
#' \donttest{data(data.m)
#' data(Mvalue)
#' E.result <- Emas(data.m, Mvalue, id = "ID", x = "x", y = "y",
#'                  x.cov = c("age", "gender"),
#'                  y.cov = c("age", "gender"),
#'                  m.cov = c("age", "gender", "CD8T", "CD4T"), 
#'                  ini.sims = 100, boot = FALSE, cl.n = 1)}
Emas <- function(data, M.matrix, id = "", x = "", y = "", x.cov = c(), y.cov = c(), m.cov = c(), ini.sims = 100, boot = FALSE, cl.n = 1, ...){
  AME.P <- NULL
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
  MV.z$ID <- NULL
  
  message(paste("A total of", dim(data.z)[1], "were included in the EMAS."), sep = " ")
  
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
  
  m.c <- m.cov[1]
  for (i in m.cov[-1]) {
    m.c <- paste(m.c, i, sep = " + ")
  }
  formula.m <- paste("m", " ~ ", m.c, sep = "")
  
  funmediation <- function(m){
    data.m <- cbind(data.z, m)
    modm <- residuals(lm(formula.m, data = data.m))
    modhe <- cbind(modx, modm, mody)
    modhe <- data.frame(modhe)
    lm.m <- lm(modm ~ modx, data = modhe)
    lm.y <- lm(mody ~ modx + modm, data = modhe)
    out.1 <- mediate(lm.m, lm.y, sims = sim, boot = boot, 
                     treat = "modx", mediator = "modm", ...)
    pcom <- c(out.1$d0, out.1$d0.ci, out.1$d0.p, 
              out.1$z0, out.1$z0.ci, out.1$z0.p, 
              out.1$tau.coef, out.1$tau.ci, out.1$tau.p, 
              out.1$n0, out.1$n0.ci, out.1$n0.p)
    names(pcom)<-c("AMEEst", "AMElow95", "AMEupp95", "AME.P",
                   "ADEEst", "ADElow95", "ADEupp95", "ADE.P",
                   "TotEst", "Totlow95", "Totupp95", "Tot.P",
                   "PropEst", "Proplow95", "Propupp95", "Prop.P")
    return(pcom)
  }
  
  cl <- makeCluster(cl.n)
  clusterExport(cl, "data.z", envir = environment())
  clusterExport(cl, "formula.m", envir = environment())
  clusterExport(cl, "modx", envir = environment())
  clusterExport(cl, "mody", envir = environment())
  clusterExport(cl, "mediate", envir = environment())
  clusterExport(cl, "boot", envir = environment())
  sim <- ini.sims
  clusterExport(cl, "sim", envir = environment())
  pzhi <- parSapply(cl, MV.z, funmediation)
  stopCluster(cl)
  
  pzhi <- t(pzhi)
  pzhi <- data.frame(pzhi)
  
  pzhi.a <- pzhi[order(pzhi$AME.P),]
  pzhi.b <- pzhi.a
  i <- pzhi.a[1,4]
  shu <- 1
  
  while (i==0) {
    pzhi.b0 <- subset(pzhi.b, AME.P==0)
    pzhi.b1 <- subset(pzhi.b, AME.P>0)
    cn <- rownames(pzhi.b0)
    MV.zf <- MV.z[,cn]
    cl <- makeCluster(cl.n)
    clusterExport(cl, "data.z", envir = environment())
    clusterExport(cl, "formula.m", envir = environment())
    clusterExport(cl, "modx", envir = environment())
    clusterExport(cl, "mody", envir = environment())
    clusterExport(cl, "mediate", envir = environment())
    clusterExport(cl, "boot", envir = environment())
    sim <- ini.sims*10^shu
    clusterExport(cl, "sim", envir = environment())
    pzhi.s <- parSapply(cl, MV.zf, funmediation)
    stopCluster(cl)
    
    pzhi.s <- t(pzhi.s)
    pzhi.s <- data.frame(pzhi.s)
    pzhi.b <- rbind(pzhi.s, pzhi.b1)
    
    pzhi.b <- pzhi.b[order(pzhi.b$AME.P),]
    i <- pzhi.b[1,4]
    print(shu)
    shu <- shu+1
    
    if(i!=0){
      cl <- makeCluster(cl.n)
      clusterExport(cl, "data.z", envir = environment())
      clusterExport(cl, "formula.m", envir = environment())
      clusterExport(cl, "modx", envir = environment())
      clusterExport(cl, "mody", envir = environment())
      clusterExport(cl, "mediate", envir = environment())
      clusterExport(cl, "boot", envir = environment())
      sim <- ini.sims*10^shu
      clusterExport(cl, "sim", envir = environment())
      pzhi.s <- parSapply(cl, MV.zf, funmediation)
      stopCluster(cl)
      
      pzhi.s <- t(pzhi.s)
      pzhi.s <- data.frame(pzhi.s)
      pzhi.b <- rbind(pzhi.s, pzhi.b1)
      
      pzhi.b <- pzhi.b[order(pzhi.b$AME.P),]
    }
  }
  return(pzhi.b)
}
