#' @title Plotting the volcano plot from the EMAS results
#'
#' @description Function to plot a volcano plot from the \code{Emas} results.
#'
#' @details This function can plot a volcano plot from the \code{Emas} results.
#'
#' @param E.result A data.frame produced by \code{Emas}.
#'
#' @param epiwideline Where to draw a "epigenome-wide sigificant" line. Default -log10(1.0e-7).
#'
#' @param suggestiveline Where to draw a "suggestive" line. Default -log10(1.0e-5). Set to FALSE to disable.
#'
#' @export
#'
#' @return No return value, called for side effects.
#'
#' @importFrom "ggplot2" "aes" "element_text" "ggplot" "xlab" "ylab" "scale_color_manual" "geom_point" "ylim" "xlim" "geom_vline" "geom_hline" "theme"
#' 
#' @author Xiuquan Nie, niexiuquan1995@foxmail.com
#' 
#' @examples
#' \donttest{data(E.result)
#' EMAS.volcano(E.result, 
#'              epiwideline = -log10(0.05/2000),
#'              suggestiveline = -log10(1/100))}
EMAS.volcano <- function(E.result, 
                         epiwideline = -log10(1.0e-07), 
                         suggestiveline = -log10(1.0e-05)){
  AMEEst <- AME.P <- threshold <- NULL
  yl <- ifelse(epiwideline > -log10(min(E.result$AME.P)), 
               ifelse(epiwideline > suggestiveline, 
                      epiwideline, 
                      suggestiveline), 
               ifelse(suggestiveline > -log10(min(E.result$AME.P)), 
                      suggestiveline, 
                      -log10(min(E.result$AME.P))))
  xl <- ifelse(abs(max(E.result$AMEEst)) > abs(min(E.result$AMEEst)), 
               abs(max(E.result$AMEEst)), 
               abs(min(E.result$AMEEst)))
  E.result$threshold <- as.factor(ifelse(E.result$AME.P < (10^(-epiwideline)), 
                                         ifelse(E.result$AMEEst >= 0 , 
                                                'Up', 'Down'), 
                                         'Not'))
  if(suggestiveline != FALSE){
    ggplot(E.result, aes(x = AMEEst, y = -log10(AME.P), color = threshold)) + 
      xlab("AMEEst") + ylab("-log10(AME.P)") + 
      geom_point(size = 2, alpha = 1) + 
      ylim(0, yl) + xlim(-xl, xl) + 
      scale_color_manual(values = c(Down = "blue", Not = "grey", Up = "red")) + 
      geom_vline(xintercept = c(0), lty = 2, color = "#000000") + 
      geom_hline(yintercept = c(epiwideline), lty = 2, colour = "red") + 
      geom_hline(yintercept = c(suggestiveline), lty = 2, color = "blue") + 
      theme(axis.text = element_text(size = 15), 
            axis.title = element_text(size = 15))
  }else{
    ggplot(E.result, aes(x = AMEEst, y = -log10(AME.P), color = threshold)) + 
      xlab("AMEEst") + ylab("-log10(AME.P)") + 
      geom_point(size = 2, alpha = 1) + 
      ylim(0, yl) + xlim(-xl, xl) + 
      scale_color_manual(values = c(Down = "blue", Not = "grey", Up = "red")) + 
      geom_vline(xintercept = c(0), lty = 2, color = "#000000") + 
      geom_hline(yintercept = c(epiwideline), lty = 2, colour = "red") + 
      theme(axis.text = element_text(size = 15), 
            axis.title = element_text(size = 15))
  }
}
