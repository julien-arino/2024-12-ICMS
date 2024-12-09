## ----set-options,echo=FALSE,warning=FALSE,message=FALSE-----------------------
# Load required libraries
required_packages = c("deSolve",
                      "dplyr", 
                      "ggplot2", 
                      "knitr", 
                      "latex2exp",
                      "lattice",
                      "magick",
                      "readr", 
                      "tidyr",
                      "viridis")

for (p in required_packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p, dependencies = TRUE)
    require(p, character.only = TRUE)
  }
}
# Knitr options
opts_chunk$set(echo = TRUE, 
               warning = FALSE, 
               message = FALSE, 
               dev = c("pdf", "png"),
               fig.width = 6, 
               fig.height = 4, 
               fig.path = "FIGS/course-02-",
               fig.keep = "high",
               fig.show = "hide")
knitr::knit_hooks$set(crop = knitr::hook_pdfcrop)
options(knitr.table.format = "latex")
# Date for front title page (if needed)
yyyy = strsplit(as.character(Sys.Date()), "-")[[1]][1]


## ----set-slide-background,echo=FALSE,results='asis'---------------------------
# Are we plotting for a dark background?
plot_blackBG = FALSE
if (plot_blackBG) {
  bg_color = "black"
  fg_color = "white"
  input_setup = "\\input{slides-setup-blackBG.tex}"
} else {
  bg_color = "white"
  fg_color = "black"
  input_setup = "\\input{slides-setup-whiteBG.tex}"
}
cat(input_setup)


## -----------------------------------------------------------------------------
pop = c(34.017, 1348.932, 1224.614, 173.593, 93.261) * 1e+06
countries = c("Canada", "China", "India", "Pakistan", "Philippines")
T = matrix(data = c(0, 1268, 900, 489, 200, 
                    1274, 0, 678, 859, 150, 
                    985, 703, 0, 148, 58, 
                    515, 893, 144, 0, 9, 
                    209, 174, 90, 2, 0), 
           nrow = 5, ncol = 5, byrow = TRUE)	


## -----------------------------------------------------------------------------
p = list()
p$M = mat.or.vec(nr = dim(T)[1], nc = dim(T)[2])
for (from in 1:5) {
  for (to in 1:5) {
    p$M[to, from] = -log(1 - T[from, to]/pop[from])
  }
  p$M[from, from] = 0
}
p$M = p$M - diag(colSums(p$M))


## -----------------------------------------------------------------------------
p$P = dim(p$M)[1]
p$epsilon = rep((1/1.5), p$P)
p$gamma = rep((1/5), p$P)
# The desired values for R_0
R_0 = rep(1.5, p$P)


## -----------------------------------------------------------------------------
p$idx_S = 1:p$P
p$idx_L = (p$P+1):(2*p$P)
p$idx_I = (2*p$P+1):(3*p$P)
p$idx_R = (3*p$P+1):(4*p$P)


## -----------------------------------------------------------------------------
# Set initial conditions. For example, we start with 2
# infectious individuals in Canada.
L0 = mat.or.vec(p$P, 1)
I0 = mat.or.vec(p$P, 1)
R0 = mat.or.vec(p$P, 1)
I0[1] = 2
S0 = pop - (L0 + I0 + R0)
# Vector of initial conditions to be passed to ODE solver.
IC = c(S = S0, L = L0, I = I0, R = R0)
# Time span of the simulation (5 years here)
tspan = seq(from = 0, to = 5 * 365.25, by = 0.1)


## -----------------------------------------------------------------------------
for (i in 1:p$P) {
  p$beta[i] = 
    R_0[i] / S0[i] * 1 / 
    ((1 - p$pi[i])/p$gammaI[i] + 
       p$pi[i] * p$eta[i]/p$gammaA[i])
}


## -----------------------------------------------------------------------------
SLIAR_metapop_rhs <- function(t, x, p) {
	with(as.list(p), {
		S = x[idx_S]
		L = x[idx_L]
		I = x[idx_I]
		A = x[idx_A]
		R = x[idx_R]
		N = S + L + I + A + R
		Phi = beta * S * (I + eta * A) / N
		dS = - Phi + M %*% S
		dL = Phi - epsilon * L + p$M %*% L
		dI = (1 - pi) * epsilon * L - gammaI * I + M %*% I
		dA = pi * epsilon * L - gammaA * A + M %*% A
		dR = gammaI * I + gammaA * A + M %*% R
		dx = list(c(dS, dL, dI, dA, dR))
		return(dx)
	})
}	


## -----------------------------------------------------------------------------
# Call the ODE solver
sol <- ode(y = IC, 
			times = tspan, 
			func = SLIAR_metapop_rhs, 
			parms = p,
			method = "ode45")


## ----convert-Rnw-to-R,warning=FALSE,message=FALSE-----------------------------
# From https://stackoverflow.com/questions/36868287/purl-within-knit-duplicate-label-error
rmd_chunks_to_r_temp <- function(file){
  callr::r(function(file, temp){
    out_file = sprintf("../CODE/%s", gsub(".Rnw", ".R", file))
    knitr::purl(file, output = out_file, documentation = 1)
  }, args = list(file))
}
rmd_chunks_to_r_temp("course-02-metapopulations-and-advanced-models.Rnw")

