# NBBQTSM
This repository contains the accompanying toolkit for "A New Bayesian Bootstrap for Quantitative Trade and Spatial Models" by Bas Sanders.

The folder "Vignette" illustrates the procedure. For an HTML page, download "Vignette.html" and open in your browser.

The folder "Code" contains the necessary functions for the Bayesian bootstrap procedure.

- The function "BB.m" is the main function. It takes as inputs a dataframe, a moment function, the polyadic order, an initial value for the structural parameter, and the number of bootstrap draws. It outputs the GMM point estimate for the structural parameter and a matrix of Bayesian bootstrap estimates.
- The functions "gmmObj.m", "justIdObj.m" and "varMat.m" are auxiliary functions.


