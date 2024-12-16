# Longitudinal and ordinal data analyses of cross-sectional, binary outcomes 

The R scripts in this repository contain the code necessary to reproduce simulations and data analysis of the manuscript "Longitudinal and ordinal data analyses of cross-sectional, binary outcomes" (Di Gravio, Tao, Nam, ??, Schildcrout). The repository contains two folder:

* [Simulation Studies](https://github.com/ChiaraDG/ordinal_clovers/tree/main/Simulation%20Studies): contains the code necessary to replicate the simulations in the paper

* [CLOVERS](https://github.com/ChiaraDG/ordinal_clovers/tree/main/CLOVERS): contains the code necessary to replicate the analysis of the CLOVERS study

More detailed instruction on how to reproduce the results are provided in each folder.

# Example

For each simulation scenario, there is a separate R file that specify the parameters used in the simulation. For instance:

```
N <- 1000
£ # parameters in the simulation study
```

The data can be done using the following code:

```
dat <- dataGeneratingFunction()
```

Now, the user can decide which modelling procedure to use:

* **Cross-sectional Binary Outcome**

```
code
```

* **Longitudinal Binary Outcome**

```
code
```

* **Cross-sectional Ordinal Outcome**

```
code
```

* **Longitudinal Ordinal Outcome**

```
code
```
