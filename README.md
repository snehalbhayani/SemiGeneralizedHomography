# Calibrated and Partially Calibrated Semi-generalized Homographies
This code library provides a MATLAB implementation of our proposed minimal solvers for computing semi-generalised homographies for calibrated and partially-calibrated cameras :
1. Calibrated case : sh5_2, sh5_3, sh5_4, sh4.5_2, sh4.5_3
2. Partially calibrated case : sh5f_2, sh5f_3
Additionally, we also have two absolute pose solvers designed for our homography setup :
1. Calibrated case : P3P+N
2. Partially calibrated case : P5Pf+N

# Executing the Minimal Solvers
* Input : Each of our solvers require a tuple, `(q, p, c)`, as the input where `q <-> (p, c)` denotes the 2D-2D point correspondences. 
* `q` denotes an array of size `3x5`, of 5 vectors, each denoting the viewing ray for the corresponding 2d image observation by the query pinhole camera.
* `p` denotes an array of size `3x5`, of 5 vectors, each denoting the viewing ray for the corresponding 2d image observation in the coordinate system of the global generalized camera system.
* Each `p` ray is accompanied with camera center `c`, or the position of the center of the pinhole camera within the generalized camera system.
* Output : The output is the generalized semi-generalized homographies,`Hs` and the corresponding plane vectors, `Nss`.

Please refer to our paper for more details.

* We have used the standard approaches for a homography decomposition to extract the relative pose and the scale. 
* An example of this, and a sample test script on synthetic scenes can be found in the synthetic_scenes/* folder, for both calibrated as well partially calibrated cameras.
* The solvers for now are released in the MATLAB programming language while the C++ version will soon be released.

# Citing
Please cite our paper : 
```@article{DBLP:journals/corr/abs-2103-06535,
  author    = {Snehal Bhayani and
               Torsten Sattler and
               Daniel Barath and
               Patrik Beliansky and
               Janne Heikkil{\"{a}} and
               Zuzana Kukelova},
  title     = {Calibrated and Partially Calibrated Semi-Generalized Homographies},
  journal   = {CoRR},
  volume    = {abs/2103.06535},
  year      = {2021},
  url       = {https://arxiv.org/abs/2103.06535},
  eprinttype = {arXiv},
  eprint    = {2103.06535},
  timestamp = {Tue, 16 Mar 2021 11:26:59 +0100},
  biburl    = {https://dblp.org/rec/journals/corr/abs-2103-06535.bib},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}
```


