# Cross-Spectral Estimator of Low-Frequency Hydrodynamic Parameters and Optimal Experimental Design

MATLAB code for numerical analysis and optimal experimental design (OED). The repository provides separate entry points for analysing a numerical model and determining an experimental design.

Optimal experimental design aims to choose experimental conditions that make the resulting measurements informative for a specified estimation or modelling task.

## Main scripts

| Script | Purpose |
| --- | --- |
| [`script_analysis_num.m`](script_analysis_num.m) | Main entry point for numerical analysis. |
| [`script_OED_num.m`](script_OED_num.m) | Main entry point for optimal experimental design. |
| [`script_LF_est_exp.m`](script_LF_est_exp.m) | Main script for data analysis and estimation based on experimental data. |

## Numerical analysis

Use `script_analysis_num.m` for the analysis workflow:

```matlab
script_analysis_num
```

Review the model definition, parameter values, and numerical settings in the source before adapting the analysis to a new problem.

## Optimal experimental design

Use `script_OED_num.m` for the optimal-design workflow:

```matlab
script_OED_num
```

## Parameter estimation from experimental data

Use `script_LF_est_exp.m` for the parameter estimation from experimental data:

```matlab
script_LF_est_exp.m
```

When adapting the design problem, review the model, experimental design variables, objective, constraints, and optimisation settings in the source. Interpret the resulting design in the context of these choices.

## Reproducibility

For each study, record the repository commit, MATLAB release, required toolbox versions, and configuration used. If a workflow uses random sampling, record the random seed as well.
