# Reliability Analysis: A Unified Tool for Assessing Measurement Reliability

This repository provides a comprehensive MATLAB toolbox for conducting reliability analysis on neuroimaging data and other datasets, whether discrete or continuous. The toolbox implements Krippendorff's Alpha and related reliability metrics, supporting a variety of data types and analysis scenarios.

The methodology and implementation details are described in the companion paper:

> Vinding, M. C. (2025). A Unified Framework for Reliability Analysis in Neuroimaging With Krippendorff’s α. *International Journal of Imaging Systems and Technology, 35(5)*, e70192. https://doi.org/10.1002/ima.70192

Consult the companion paper for guidance on the available methods, their application, performance, and interpretation.

## Usage

The primary function is `reliability_analysis`, which computes Krippendorff's Alpha using the appropriate method for your data type.

**Basic usage:**
```matlab
ALPHA = reliability_analysis(DATA, METHOD)
```
- `DATA`: An N x M matrix, where *N* is the number of "observers" (*N* > 1) and *M* is the number of measurements or data points (e.g., time samples in time series).
- `METHOD`: A string specifying the data type. Valid options:
  - `'nominal'`
  - `'ordinal'`
  - `'interval'`
  - `'ratio'`
  - `'phase_rad'` or `'phase_deg'`

For datasets with exactly two observers (*N* = 2) and interval data, you can use the optimized fast algorithm:
```matlab
ALPHA = reliability_analysis(DATA, 'N2fast')
```

**Bootstrapping confidence intervals:**

You can obtain bootstrap confidence intervals for Alpha values as described by Hayes & Krippendorff (2007):
```matlab
[ALPHA, BOOTS] = reliability_analysis(DATA, METHOD, BOOTSTRAP)
```
- `BOOTSTRAP`: Number of bootstrap samples (e.g., 1000). Set to `0` for no bootstrapping.
- `BOOTS`: A vector of bootstrapped Alpha values.

## Examples

The `Examples` folder contains scripts demonstrating how to calculate Alpha for different data types. These include:
- `example_nominal.m`: Reliability analysis for nominal-scaled data.
- `example_interval.m`: Reliability analysis for interval-scaled data.
- `example_phase.m`: Reliability analysis for circular (phase) data.

Use these scripts as templates for analyzing your own data. All examples use the `reliability_analysis` function and correspond directly to the code provided in this repository.

## Toolbox Contents

The main function (`reliability_analysis.m`) is a wrapper that delegates computation to specialized functions:
- `kripAlpha.m`: Exact calculation of Krippendorff's Alpha for interval, ordinal, nominal, ratio, or phase data.
- `alphaprime.m`: Approximate Alpha for large datasets with arbitrary numerical precision using data binning.
- `kripAlpha2fast.m`: Fast, exact calculation for interval or ordinal data with two observers (*N* = 2).
- `bootstrap_alpha.m`: Bootstrapping procedure for estimating confidence intervals.

See the documentation in each function and the companion paper for further details.

## Background

For more information on Krippendorff's Alpha and its application, see:
- Hayes, A. F. & Krippendorff, K. (2007). Answering the call for a standard reliability measure for coding data. *Communication Methods and Measures*, 1, 77–89.
- Krippendorff, K. (2018). *Content Analysis: An Introduction to Its Methodology* (Fourth Edition). SAGE.
- Vinding, M. C. (2025). A Unified Framework for Reliability Analysis in Neuroimaging With Krippendorff’s α. *International Journal of Imaging Systems and Technology, 35(5)*, e70192. https://doi.org/10.1002/ima.70192

## Contact and Contribution

The code is actively maintained. If you have suggestions for improvements or requests for additional features, please contact me via email or open a GitHub issue. All feedback and contributions are welcome!

**Email:** mvi@psy.ku.dk
