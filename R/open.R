# ==========================================================================
# EMAS package initialization
# ==========================================================================

.onAttach <- function(libname, pkgname) {
  mes <- "       _________    ___    __        _        _________
      |   ______|  |   \\  /  |      / \\      |   ______|
      |  |______   |    \\/   |     /   \\     |  |______
      |   ______|  |  |\\  /| |    /  /\\ \\    |_______  |
      |  |______   |  | \\/ | |   /  ____ \\    _______| |
      |_________|  |__|    |_|  /__/    \\_\\  |_________|
     ------------------------------
"
  mes2 <- "    If you have any question or suggestion about EMAS, please email to niexiuquan1995@foxmail.com.\n     ------------------------------\n"
  packageStartupMessage(">> Package version 0.2.2 loaded <<\n",mes,mes2)
}

