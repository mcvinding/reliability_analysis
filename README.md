# Reliability Analysis: a tool for assessing reliability between measurements.

This repository provides a unified tool for reliability tests of neuroimaging data of various types that can be either discrete or continuous data. The test calculates Krippendorff's Alpha or an approximation that is analogous to Krippendorff's Alpha. All code is implemented in MATLAB and requires no further dependencies. This implementation is designed to run fast and handle large datasets, which is necessary for neuroimaging data such as (f)MRI, EEG, MEG, PET, etc.

The methods and implementation is described in the companion paper:

  *TBA*

Please refer to the companion paper for information on the implemented methods, their usage, performance, and interpretation.

## Usage
Use the main function `reliability_analysis` to calculate Alpha with the appropriate type of test.

Use as:
````Matlab
ALPHA = reliability_analysis(DATA, METHOD)
````

Where `DATA` is the reliability data in the format *N* x *M*, where *N* is the number of "observers" (*N* > 2) and *M* is the number of data points. For neural time series data *M* = time. You might need to concatenate or "flatten" your data, so all observation for one observer is one vector. Each row in the data frame represents one observer. Each column represents the same unit of observation. E.g., the first column could be the first time point across all observers.

Use `METHOD` to declare the type of data: `nominal`, `ordinal`, `interval`, `ratio`, and `phase`, or use either the reduced but faster algorithm for cases where *N* = 2 and data is interval ( `N2fast` ) or the approximation algorithm ( `alphaprime` ), which is useful for large datasets with arbitrary numerical precision. For more information about methods, please refer to the companion paper (TBA) or function documentation.

It is also possible to get bootstrap confidence intervals of Alpha values based on the procedure described by Hayes & Krippendorff (2007) in the following way:
````Matlab
[ALPHA, BOOTS] = reliability_analysis(DATA, METHOD, BOOTSTRAP)
````
Where `BOOTSTRAP` indicates the size of the bootstrapping distribution (`BOOTSTRAP = 0` means no bootstrapping) and `BOOTS` is a vector of the bootstrapped alpha values.

## Examples
The folder Examples provides three examples of how to calculate Alpha for different types of data. Use these scripts as blueprints to start reliability analysis on your own data.

## Content
The main function is a wrapper that calls the different functions to calculate alpha and the bootstrapping procedure. For more documentation and options see the individual functions:
* `kripAlpha.m` : Krippendorff's Alpha for *interval*, *ordinal*, *nominal*, *ratio*, or *phase* data using exact caluclation of Alpha.
* `alphaprime.m` : Approximation of Krippendorff's Alpha for interval data with arbitary precision.
* `kripAlpha2fast.m` : Fast, exact calculation of Krippendorff's Alpha for interval data with two observers (*N* = 2).
* `bootstrap_alpha.m` : Run the bootstrapping procedure based on the output from either of the functions above. 

See the documentation and the companion paper for more details on the implementations.

## Background
For more information on Krippendorff's Alpha see: 
* Hayes, A.F. & Krippendorff, K. (2007). Answering the call for a standard reliability measure for coding data. *Communication Methods and Measures*, 1, 77-89
* Krippendorff, K. (2018). *Content analysis: An introduction to its methodology* (Fourth Edition). SAGE.

