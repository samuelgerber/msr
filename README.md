# Morse Smale Decompositions for Regression and Visualization

Implementation of Morse-Smale complex and Morse-Smale regression in C++ with an R package wrapper.


The R package permits to compute Morse-Smale Complex decompositions and implements regression and visualization based on those decompositions.

## Install

```R
library(devtools)
devtools::install_github("samuelgerber/msr)
```

## More Information
The R package builds interface to the C++ Morse-Smale complex computation as first described in:

Gerber S, Bremer PT, Pascucci V, Whitaker R (2010).
“Visual Exploration of High Dimensional Scalar Functions.”
IEEE Transactions on Visualization and Computer Graphics, 16(6), 1271–1280.

The main functionality is implemented in NNMSComplex.h and was implemented by S. Gerber.
The code uses the ANN library by David Mount and Sunil Arya distributed under
the LGPL license. For detailed information on ANN see
http://www.cs.umd.edu/~mount/ANN.

On this basic functionality the exploratory visualization approach in

Gerber S, Bremer PT, Pascucci V, Whitaker R (2010).
“Visual Exploration of High Dimensional Scalar Functions.”
IEEE Transactions on Visualization and Computer Graphics, 16(6), 1271–1280.

and the regression approach in

Morse-Smale Regression
Gerber S, Rübel O, Bremer PT, Pascucci V, Whitaker RT
J Comput Graph Stat. 2013 Jan 1;22(1):193-214.

are implemented. A detailed description of the package can be found in:

Data analysis with the morse-smale complex: The msr package for r
S Gerber, K Potter
Journal of Statistical Software, 2011
