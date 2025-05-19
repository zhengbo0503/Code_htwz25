Code for the preprint

 - N. J. Higham, F. Tisseur, M. Webb, and Z. Zhou,
   Computing accurate eigenvalues using a mixed-precision Jacobi algorithm,
   arXiv:2501.03742, January 2025.
   Available at https://arxiv.org/abs/2501.03742.


## Usage 

### MATLAB 
Before executing tests from `MATLAB/test`, you need to add `MATLAB/src` into MATLAB search path by 

```matlab
addpath("../src/")
```

### Julia 
1. First you need to clone [JacobiEigen.jl](https://github.com/zhengbo0503/JacobiEigen.jl) to your working directory. 
2. In terminal, open Julia, type `]`, and add the package `JacobiEigen.jl` into your Julia environment.
   ```julia 
   (@v1.11) pkg> add "/path/to/the/package/"
   ```
3. Then you may execute the Julia code by 
   ```bash 
   julia ./test.jl
   ```
