#' @title Plotting the manhattan plot from the EMAS results
#'
#' @description Function to plot a manhattan plot from the \code{Emas} results.
#'
#' @details This function can plot a manhattan plot from the \code{Emas} results according to the annotation from 450k or EPIC.
#'
#' @param E.result A data.frame produced by \code{Emas}.
#'
#' @param type A character string indicating the type of annotation, only "EPIC" and "450k" are available. 
#'
#' @param ... Other arguments passed to \code{\link{manhattan}}.
#'
#' @export
#'
#' @return No return value, called for side effects.
#'
#' @importFrom "minfi" "getAnnotation"
#' @importFrom "qqman" "manhattan"
#' 
#' @author Xiuquan Nie, niexiuquan1995@foxmail.com
#' 
#' @references
#' Turner, (2018). qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots.
#' \emph{Journal of Open Source Software}, 3(25), 731. \doi{doi:10.21105/joss.00731}.
#' 
#' @examples
#' \donttest{data(E.result)
#' EMAS.manhattan(E.result, type = "EPIC",
#'                genomewideline = -log10(0.05/2000),
#'                suggestiveline = -log10(1/100), ylim=c(0,5))}
EMAS.manhattan <- function(E.result, type = "EPIC", ...){
  CpG <- chr <- pos <- AME.P <- NULL
  if(type != "450k" & type != "EPIC"){
    stop("At present, this function can only draw Manhattan diagrams of 450k and EPIC.")
  }
  if(type == "EPIC"){
    annEPIC <- getAnnotation("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
    ann <- annEPIC
    ann$CpG <- row.names(ann)
    E.result$CpG <- row.names(E.result)
    manhz <- merge(E.result, ann, by="CpG", all.x=T)
    manh <- subset(manhz, select=c(CpG, chr, pos, AME.P))
    manh$CHR <- as.numeric(substr(manh$chr, 4, nchar(manh$chr))) 
    manhattan(manh, chr = "CHR",bp="pos", snp = "CpG", p="AME.P", ...)
  }
  if(type == "450k"){
    ann450k <- getAnnotation("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    ann <- ann450k
    ann$CpG <- row.names(ann)
    E.result$CpG <- row.names(E.result)
    manhz <- merge(E.result, ann, by="CpG", all.x=T)
    manh <- subset(manhz, select=c(CpG, chr, pos, AME.P))
    manh$CHR <- as.numeric(substr(manh$chr, 4, nchar(manh$chr))) 
    manhattan(manh, chr = "CHR",bp="pos", snp = "CpG", p="AME.P", ...)
  }
}