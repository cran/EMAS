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
#' @param mem.sav A logical value. If 'TRUE', the memory required for the function will decrease, but the speed will also decrease.
#'
#' @param p.th Sobel indirect effects P-value threshold for subsequent nonparametric bootstrap or quasi-Bayesian approximation mediation analyses.
#'
#' @param ini.sims Initial number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation.
#'
#' @param boot A logical value. If 'FALSE' a quasi-Bayesian approximation is used for confidence intervals; if 'TRUE' nonparametric bootstrap will be used. Default is 'FALSE'.
#'
#' @param cl.n Number of cores used for parallel computing.
#' 
#' @param ... Other arguments passed to \code{\link{makeCluster}}.
#'
#' @export
#'
#' @return \code{Emas} returns a data.frame with the average mediation effects(AME), 
#'    average direct effects(ADE), total effects, mediation proportion.
#' \itemize{
#'     \item{AMEEst: }{Point estimates for average mediation effects under the exposure conditions.}
#'     \item{AMElow95, AMEupp95: }{95 percentage confidence intervals for average mediation effects.}
#'     \item{AME.P: }{Two-sided p-values for average mediation effects.}
#'     \item{ADEEst: }{Point estimates for average direct effect under the exposure conditions.}
#'     \item{ADElow95, ADEupp95: }{95 percentage confidence intervals for average direct effects.}
#'     \item{ADE.P: }{Two-sided p-values for average direct effects.}
#'     \item{TotEst: }{Point estimate for total effect.}
#'     \item{Totlow95, Totupp95: }{95 percentage confidence interval for total effect.}
#'     \item{Tot.P: }{Two-sided p-values for total effect.}
#'     \item{PropEst: }{The "proportions mediated", or the size of the average mediation effects relative to the total effect.}
#' }
#' 
#' @importFrom "mediation" "mediate"
#' @importFrom "parallel" "makeCluster" "clusterExport" "parSapply" "stopCluster"
#' @importFrom "stats" "lm" "na.omit" "residuals" "pnorm"
#' @importFrom "multilevel" "sobel"
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
#'                  p.th = 0.1, ini.sims = 100, boot = FALSE, cl.n = 1)}
Emas <- function(data, M.matrix, 
                 id = "", x = "", y = "", 
                 x.cov = c(), y.cov = c(), m.cov = c(), 
                 mem.sav = FALSE, p.th = 0.1, ini.sims = 100, 
                 boot = FALSE, cl.n = 1, ...){
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
  snum <- dim(data.z)[1]
  MV.z$ID <- NULL
  
  if(mem.sav == TRUE){
    rm("M.matrix", "zzM", "Mv", "MV")
    gc()
  }
  
  message(paste("A total of", snum, "were included in the EMAS.", sep = " "))
  
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
  
  if(mem.sav == TRUE){
    funsobel <- function(m){
      data.m <- cbind(data.z, m)
      modm <- residuals(lm(formula.m, data = data.m))
      if(length(modm) == snum){
        modhe <- cbind(modx, modm, mody)
        modhe <- data.frame(modhe)
        sobelre <- sobel(pred = modhe$modx, 
                         med = modhe$modm, 
                         out = modhe$mody)
        pcom <- c(sobelre$Indirect.Effect, 
                  sobelre$Indirect.Effect - 1.96*sobelre$SE, 
                  sobelre$Indirect.Effect + 1.96*sobelre$SE, 
                  2*pnorm(-abs(sobelre$z.value)), 
                  sobelre$`Mod2: Y~X+M`[2,1], 
                  sobelre$`Mod2: Y~X+M`[2,1] - 1.96*sobelre$`Mod2: Y~X+M`[2,2], 
                  sobelre$`Mod2: Y~X+M`[2,1] + 1.96*sobelre$`Mod2: Y~X+M`[2,2], 
                  sobelre$`Mod2: Y~X+M`[2,4], 
                  sobelre$`Mod1: Y~X`[2,1],
                  sobelre$`Mod1: Y~X`[2,1] - 1.96*sobelre$`Mod1: Y~X`[2,2],
                  sobelre$`Mod1: Y~X`[2,1] + 1.96*sobelre$`Mod1: Y~X`[2,2],
                  sobelre$`Mod1: Y~X`[2,4],
                  sobelre$Indirect.Effect/sobelre$`Mod1: Y~X`[2,1])
        rm("data.m", "modm", "sobelre", "modhe")
        gc()
        names(pcom) <- c("AMEEst", "AMElow95", "AMEupp95", "AME.P",
                         "ADEEst", "ADElow95", "ADEupp95", "ADE.P",
                         "TotEst", "Totlow95", "Totupp95", "Tot.P",
                         "PropEst")
        return(pcom)
      }
    }
  }else{
    funsobel <- function(m){
      data.m <- cbind(data.z, m)
      modm <- residuals(lm(formula.m, data = data.m))
      if(length(modm) == snum){
        modhe <- cbind(modx, modm, mody)
        modhe <- data.frame(modhe)
        sobelre <- sobel(pred = modhe$modx, 
                         med = modhe$modm, 
                         out = modhe$mody)
        pcom <- c(sobelre$Indirect.Effect, 
                  sobelre$Indirect.Effect - 1.96*sobelre$SE, 
                  sobelre$Indirect.Effect + 1.96*sobelre$SE, 
                  2*pnorm(-abs(sobelre$z.value)), 
                  sobelre$`Mod2: Y~X+M`[2,1], 
                  sobelre$`Mod2: Y~X+M`[2,1] - 1.96*sobelre$`Mod2: Y~X+M`[2,2], 
                  sobelre$`Mod2: Y~X+M`[2,1] + 1.96*sobelre$`Mod2: Y~X+M`[2,2], 
                  sobelre$`Mod2: Y~X+M`[2,4], 
                  sobelre$`Mod1: Y~X`[2,1],
                  sobelre$`Mod1: Y~X`[2,1] - 1.96*sobelre$`Mod1: Y~X`[2,2],
                  sobelre$`Mod1: Y~X`[2,1] + 1.96*sobelre$`Mod1: Y~X`[2,2],
                  sobelre$`Mod1: Y~X`[2,4],
                  sobelre$Indirect.Effect/sobelre$`Mod1: Y~X`[2,1])
        names(pcom) <- c("AMEEst", "AMElow95", "AMEupp95", "AME.P",
                         "ADEEst", "ADElow95", "ADEupp95", "ADE.P",
                         "TotEst", "Totlow95", "Totupp95", "Tot.P",
                         "PropEst")
        return(pcom)
      }
    }
  }
  cl <- makeCluster(cl.n, ...)
  clusterExport(cl, "data.z", envir = environment())
  clusterExport(cl, "snum", envir = environment())
  clusterExport(cl, "formula.m", envir = environment())
  clusterExport(cl, "modx", envir = environment())
  clusterExport(cl, "mody", envir = environment())
  clusterExport(cl, "sobel", envir = environment())
  pzhi <- parSapply(cl, MV.z, funsobel)
  stopCluster(cl)
  
  if(class(pzhi)[1] == "list"){
    pzhi <- pzhi[!sapply(pzhi,is.null)]
  }
  
  pzhi <- data.frame(pzhi)
  pzhi <- data.frame(t(pzhi))
  message(paste("A total of", dim(pzhi)[1], "were included in the EMAS.", sep = " "))
  pzhi.b <- pzhi[order(pzhi$AME.P),]
  i <- pzhi.b[1,4]
  
  if(mem.sav == TRUE){
    funmediation <- function(m){
      data.m <- cbind(data.z, m)
      modm <- residuals(lm(formula.m, data = data.m))
      modhe <- cbind(modx, modm, mody)
      modhe <- data.frame(modhe)
      lm.m <- lm(modm ~ modx, data = modhe)
      lm.y <- lm(mody ~ modx + modm, data = modhe)
      out.1 <- mediate(lm.m, lm.y, sims = sim, boot = boot, 
                       treat = "modx", mediator = "modm")
      pcom <- c(out.1$d0, out.1$d0.ci, out.1$d0.p, 
                out.1$z0, out.1$z0.ci, out.1$z0.p, 
                out.1$tau.coef, out.1$tau.ci, out.1$tau.p, 
                out.1$n0)
      rm("data.m", "modm", "modhe", "lm.m", "lm.y", "out.1")
      gc()
      names(pcom) <- c("AMEEst", "AMElow95", "AMEupp95", "AME.P",
                       "ADEEst", "ADElow95", "ADEupp95", "ADE.P",
                       "TotEst", "Totlow95", "Totupp95", "Tot.P",
                       "PropEst")
      return(pcom)
    }
  }else{
    funmediation <- function(m){
      data.m <- cbind(data.z, m)
      modm <- residuals(lm(formula.m, data = data.m))
      modhe <- cbind(modx, modm, mody)
      modhe <- data.frame(modhe)
      lm.m <- lm(modm ~ modx, data = modhe)
      lm.y <- lm(mody ~ modx + modm, data = modhe)
      out.1 <- mediate(lm.m, lm.y, sims = sim, boot = boot, 
                       treat = "modx", mediator = "modm")
      pcom <- c(out.1$d0, out.1$d0.ci, out.1$d0.p, 
                out.1$z0, out.1$z0.ci, out.1$z0.p, 
                out.1$tau.coef, out.1$tau.ci, out.1$tau.p, 
                out.1$n0)
      names(pcom) <- c("AMEEst", "AMElow95", "AMEupp95", "AME.P",
                       "ADEEst", "ADElow95", "ADEupp95", "ADE.P",
                       "TotEst", "Totlow95", "Totupp95", "Tot.P",
                       "PropEst")
      return(pcom)
    }
  }
  if(i <= p.th){
    pzhi.bin <- subset(pzhi.b, AME.P <= p.th)
    pzhi.bout <- subset(pzhi.b, AME.P > p.th)
    cn <- rownames(pzhi.bin)
    MV.zf <- data.frame(MV.z[,cn])
    colnames(MV.zf) <- cn
    cl <- makeCluster(cl.n, ...)
    clusterExport(cl, "data.z", envir = environment())
    clusterExport(cl, "formula.m", envir = environment())
    clusterExport(cl, "modx", envir = environment())
    clusterExport(cl, "mody", envir = environment())
    clusterExport(cl, "mediate", envir = environment())
    clusterExport(cl, "boot", envir = environment())
    sim <- ini.sims
    clusterExport(cl, "sim", envir = environment())
    pzhi.s <- parSapply(cl, MV.zf, funmediation)
    stopCluster(cl)
    pzhi.s <- t(data.frame(pzhi.s))
    pzhi.s <- data.frame(pzhi.s)
    pzhi.s <- pzhi.s[order(pzhi.s$AME.P),]
    i <- pzhi.s[1,4]
  }
  
  shu <- 1
  while (i == 0) {
    pzhi.ss <- subset(pzhi.s, AME.P == (9 * 10^(1 - log10(ini.sims) - shu)) 
                      | AME.P == (8 * 10^(1 - log10(ini.sims) - shu)) 
                      | AME.P == (7 * 10^(1 - log10(ini.sims) - shu)) 
                      | AME.P == (6 * 10^(1 - log10(ini.sims) - shu)) 
                      | AME.P == (5 * 10^(1 - log10(ini.sims) - shu)) 
                      | AME.P == (4 * 10^(1 - log10(ini.sims) - shu)) 
                      | AME.P == (3 * 10^(1 - log10(ini.sims) - shu)) 
                      | AME.P == (2 * 10^(1 - log10(ini.sims) - shu)) 
                      | AME.P == (1 * 10^(1 - log10(ini.sims) - shu)))
    pzhi.s0 <- subset(pzhi.s, AME.P == 0)
    pzhi.sr <- rbind(pzhi.s0, pzhi.ss)
    cn <- rownames(pzhi.sr)
    bu <- setdiff(rownames(pzhi.s), cn)
    pzhi.s1 <- pzhi.s[bu,]
    MV.zf <- data.frame(MV.z[,cn])
    colnames(MV.zf) <- cn
    cl <- makeCluster(cl.n, ...)
    clusterExport(cl, "data.z", envir = environment())
    clusterExport(cl, "formula.m", envir = environment())
    clusterExport(cl, "modx", envir = environment())
    clusterExport(cl, "mody", envir = environment())
    clusterExport(cl, "mediate", envir = environment())
    clusterExport(cl, "boot", envir = environment())
    sim <- ini.sims * 10^shu
    clusterExport(cl, "sim", envir = environment())
    pzhi.r <- parSapply(cl, MV.zf, funmediation)
    stopCluster(cl)
    
    pzhi.r <- t(pzhi.r)
    pzhi.r <- data.frame(pzhi.r)
    pzhi.s <- rbind(pzhi.r, pzhi.s1)
    pzhi.s <- pzhi.s[order(pzhi.s$AME.P),]
    i <- pzhi.s[1,4]
    shu <- shu + 1
    
    if(i != 0){
      cn <- rownames(pzhi.s0)
      bu <- setdiff(rownames(pzhi.s), cn)
      pzhi.s1 <- pzhi.s[bu,]
      MV.zf <- data.frame(MV.z[,cn])
      colnames(MV.zf) <- cn
      cl <- makeCluster(cl.n, ...)
      clusterExport(cl, "data.z", envir = environment())
      clusterExport(cl, "formula.m", envir = environment())
      clusterExport(cl, "modx", envir = environment())
      clusterExport(cl, "mody", envir = environment())
      clusterExport(cl, "mediate", envir = environment())
      clusterExport(cl, "boot", envir = environment())
      sim <- ini.sims * 10^shu
      clusterExport(cl, "sim", envir = environment())
      pzhi.r <- parSapply(cl, MV.zf, funmediation)
      stopCluster(cl)
      
      pzhi.r <- t(pzhi.r)
      pzhi.r <- data.frame(pzhi.r)
      pzhi.s <- rbind(pzhi.r, pzhi.s1)
      pzhi.s <- pzhi.s[order(pzhi.s$AME.P),]
    }
  }
  pzhi.b <- rbind(pzhi.s, pzhi.bout)
  return(pzhi.b)
}
