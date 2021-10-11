# Calibrated and Partially Calibrated Semi-generalized Homographies

## This code includes our proposed minimal solvers for computing semi-generalised homographies for calibrated and partially-calibrated cameras :
Calibrated case : sh5_2, sh5_3, sh5_4, sh4.5_2, sh4.5_3
Partially calibrated case : sh5f_2, sh5f_3


## Additionally, we also have two absolute pose solvers designed for our homography setup :
Calibrated case : P3P+N
Partially calibrated case : P5Pf+N

# Executing the Minimal Solvers
## The minimal solvers return semi-generalized homographies and plane normal vector. But after that we use standard homography decomposition methods to extract the relative pose and the scale. 
## An example of this, and a sample test script on synthetic scenes can be found in the synthetic_scenes/* folder, for both calibrated as well partially calibrated cases.

## The solvers for now are released in MATLAB programming language.
## The C++ version will soon be released.


