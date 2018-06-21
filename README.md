# GP4C

This folder contains the code for the paper
    Variational Inference for Gaussian Processes with Panel Count Data

Dependency:
   1. GPML Matlab Code version 4.1, http://www.gaussianprocess.org/gpml/code/matlab/doc/
   2. Limited-memory projected quasi-Newton algorithm, http://proceedings.mlr.press/v5/schmidt09a/schmidt09a.pdf
   3. Safe computation of logarithm-determinat of large matrix by Dahua Lin, MIT.

The files in the main folder are main functions to run all the experiments.
   1. test_synthetic_fromGP.m
      Use a GP to draw an intensity function and use it as the underlying intensity function to test GP4C(b=0,1,0.3).
      
   2. test_synthetic_stepfunction.m
      Use a step function as the underlying intensity function to test GP4C(b=0,1,0.3).
      
   3. test_synthetic_pseudoinputs.m
      Use a step function as the underlying intensity function and vary the number of pseudo inputs.
   
   4. test_synthetic_fileration.m
      Use a step function as the underlying intensity function and vary the number of training files.
      
   5. test_synthetic_duplicate.m
      Change the number of duplicate points in the training files and compare the computational time.
   
   6. test_realworld_ver0.m
      Test GP4C,GP4CW and LocalEM on the three real-world data sets.
   
   7. test_realworld_demo.m
      A demo for the realworld dataset.
   
   8. test_additional_KS_GP3
      Compare kernel smoothing and GP3.

  The folder in util_plot contains the plotting functions for all the figures in the paper.
      
   
