## Main Function: `SignifiTestCmplxConfounders.m`

This MATLAB code is designed for **mediation analysis** where the mediator is a matrix, and both the exposure and outcome are scalars. 

The method requires installing the following toolboxes:

- **SparseReg Toolbox**: [SparseReg Toolbox Documentation](https://hua-zhou.github.io/SparseReg/)
- **TensorReg Toolbox**: [TensorReg Toolbox Documentation](https://hua-zhou.github.io/TensorReg/)

In addition, the files inside the `mex lifes.zip` needs to be extracted and added to the Matlab path

Key stages in the code include:

1. **Library and Path Setup**: Adds essential toolboxes and directories to the MATLAB path for tensor and regression operations.

2. **Data Loading and Simulation**: Defines whether to use simulated or real data. For simulated data, it generates complex data; for real data, it loads predefined datasets.

3. **Matrix and Parameter Initialization**: Sets up matrices and parameters for the regression, including specifying bootstrap samples and confounders.

4. **Regression Analysis**: Calls `Regressions_CovConf` to estimate model parameters `A` and `B`, storing results as complex matrices.

5. **Bootstrapping**: Optionally performs bootstrapping to generate multiple estimates of `A` and `B` for statistical analysis.

6. **Probability and Statistical Testing**: Calculates p-values for matrix elements of `A`, `B`, and their element-wise product `A * B`. Applies False Discovery Rate (FDR) thresholds to identify significant elements.

7. **Visualization**: Creates figures showing absolute values and significance levels of matrices `A`, `B`, and `A * B` based on calculated p-values.

8. **Timing**: Outputs the total processing time in minutes.

This code is designed for high-dimensional data processing, leveraging complex-valued matrices and bootstrapping to enhance the robustness of regression analysis in tensor-based data. It is applicable to fields such as neuroimaging and signal processing.

Reference: C. Lopez Naranjo et al., “EEG functional connectivity as a Riemannian mediator: An application to malnutrition and cognition,” Human Brain Mapping, vol. 45, no. 7, p. e26698, May 2024, doi: 10.1002/hbm.26698.



