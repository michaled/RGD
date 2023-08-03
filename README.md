# A Convex Optimization Framework for Regularized Geodesic Distances


MATLAB implementation of 
"A Convex Optimization Framework for Regularized Geodesic Distances" by 
Michal Edelstein, 
[Nestor Guillen](https://www.ndguillen.com/),
[Justin Solomon](https://people.csail.mit.edu/jsolomon/) and
[Mirela Ben-Chen](https://mirela.net.technion.ac.il/).


[Paper](https://mirelabc.github.io/publications/rgd.pdf), 
[Supplemental](https://mirelabc.github.io/publications/rgd_sup.pdf) 


![teaser](teaser.jpg)


If you have any questions or issues running the code, please contact smichale@cs.technion.ac.il


## Using Demo Scripts

We provide a few demo scripts that compute different regularized geodesic distances. You can either run the scripts locally or online using MATLAB Online. To run the scripts locally, download this repository and add the `code` folder to your MATLAB path. To run the scripts online, click the "Open in MATLAB Online" button below. This will open a new MATLAB Online project with the code. 

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=michaled/RGD&file=code/demo.m)

 - **demo** - compute regularized geodesic distances using the Dirichlet and vector field alignment regularizations. 

 - **cvx_demo** - compute regularized geodesic distances using the Dirichlet and $L_{\infty}$ regularizations.  
    To run the `cvx_demo` script, you must have [CVX](http://cvxr.com/cvx/) installed on your system. Additionally, CVX should be added to your MATLAB path.
    It is recommended to use MOSEK as the solver for CVX. For more information, see the [CVX documentation](http://cvxr.com/cvx/doc/mosek.html).

 - **hessian_demo** - compute regularized geodesic distances using the Hessian regularization.  
    To use the Hessian regularization, you must install and compile the code from [this repository](https://github.com/odedstein/ASmoothnessEnergyWithoutBoundaryDistortionForCurvedSurfaces), which implements "A Smoothness Energy without Boundary Distortion for Curved Surfaces" by Stein et al., 2020. This should result in a compiled mex file named `curved_hessian`. Make sure to add it to your MATLAB path using Matlab's `addpath` command.

 - **allpairs_demo** - compute the Dirichlet regularized geodesic distances matrix between all pairs of points.  


## Citation

If you find this code useful, please cite our paper:

```
@article{edelstein2023convex,
  title={A Convex Optimization Framework for Regularized Geodesic Distances},
  author={Edelstein, Michal and Guillen, Nestor and Solomon, Justin and Ben-Chen, Mirela},
  booktitle={ACM SIGGRAPH 2023 Conference Proceedings},
  year={2023}
}
```
