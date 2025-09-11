# Reliability Analysis: A Tool for Assessing Reliability Between Measurements

[![MATLAB](https://img.shields.io/badge/MATLAB-R2019b+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.31234%2Fosf.io%2Fptxv6-blue)](https://doi.org/10.1002/ima.70192)

This repository provides a comprehensive MATLAB toolbox for conducting reliability analysis on neuroimaging data and other datasets, whether discrete or continuous. The toolbox implements Krippendorff's Alpha and related reliability metrics, supporting a variety of data types and analysis scenarios.

The methodology and implementation details are described in the companion paper:

> Vinding, M. C. (2025). A Unified Framework for Reliability Analysis in Neuroimaging With Krippendorff’s α. *International Journal of Imaging Systems and Technology, 35(5)*, e70192. https://doi.org/10.1002/ima.70192

Consult the companion paper for guidance on the available methods, their application, performance, and interpretation.

## Citation

If you use this tool in your research, please cite the companion paper:

```bibtex
@article{vinding2025unified,
  title={A unified framework for reliability analysis in neuroimaging with Krippendorff's α},
  author={Vinding, Mikkel C.},
  journal={PsyArXiv},
  year={2025},
  doi={10.31234/osf.io/ptxv6}
}
```

## Prerequisites

- MATLAB (tested with MATLAB R2019b and later)
- No additional toolboxes required

## Installation

1. Clone this repository or download the source code:
   ```bash
   git clone https://github.com/mcvinding/reliability_analysis.git
   ```

2. Add the reliability_analysis folder to your MATLAB path:
   ```matlab
   addpath('/path/to/reliability_analysis')
   ```

3. Verify installation by running a simple example:
   ```matlab
   data = [1 2 3; 1 2 4];  % Simple 2x3 test data
   alpha = reliability_analysis(data, 'interval');
   ```

## Quick Start

For the impatient, here's the simplest way to calculate Krippendorff's Alpha:

```matlab
% Your data: N observers × M observations
data = [1 2 3 4 5;      % Observer 1
        1 2 3 4 6];     % Observer 2

% Calculate Alpha (choose appropriate data type)
alpha = reliability_analysis(data, 'interval');   % For interval data
alpha = reliability_analysis(data, 'ordinal');    % For ordinal data  
alpha = reliability_analysis(data, 'nominal');    % For nominal data

% For faster calculation with exactly 2 observers
alpha = reliability_analysis(data, 'n2fast_interval');  % Fast interval
alpha = reliability_analysis(data, 'n2fast_nominal');   % Fast nominal
```

## Usage

Use as:
```matlab
ALPHA = reliability_analysis(DATA, METHOD)
```


For datasets with exactly two observers (*N* = 2) and interval data, you can use the optimized fast algorithm:
```matlab
ALPHA = reliability_analysis(DATA, 'N2fast')
```

**Method Selection Guide:**

| Data Type | Method | Description | Use When |
|-----------|--------|-------------|----------|
| `'nominal'` | Exact | Categorical data without order | Categories like colors, names |
| `'ordinal'` | Exact | Ranked/ordered categories | Ratings like 1-5 Likert scales |
| `'interval'` | Exact | Continuous with equal intervals | Temperature in Celsius |
| `'ratio'` | Exact | Continuous with meaningful zero | Weight, height, reaction time |
| `'angle_rad'` | Exact | Circular data in radians | Phase angles, directions |
| `'angle_deg'` | Exact | Circular data in degrees | Phase angles, directions |
| `'n2fast_interval'` | Fast | Interval data, 2 observers only | Large datasets, N=2 |
| `'n2fast_nominal'` | Fast | Nominal data, 2 observers only | Large datasets, N=2 |
| `'alphaprime'` | Approximation | Very large interval datasets | When exact calculation is too slow |

**Legacy/Shorthand Methods:**
- `'n2fast'` - Same as `'n2fast_interval'` (kept for backward compatibility)
- `'prime'` - Same as `'alphaprime'`

For more information about methods, please refer to the companion paper (Vinding, 2025) and the function documentation.

It is also possible to get bootstrap confidence intervals of Alpha values based on the procedure described by Hayes & Krippendorff (2007) in the following way:
```matlab
[ALPHA, BOOTS] = reliability_analysis(DATA, METHOD, BOOTSTRAP)
```
Where `BOOTSTRAP` indicates the size of the bootstrapping distribution (`BOOTSTRAP = 0` means no bootstrapping) and `BOOTS` is an array of the bootstrapped alpha values.

## Examples

The `examples/` folder provides three comprehensive examples demonstrating how to calculate Alpha for different types of data. Use these scripts as blueprints to start reliability analysis on your own data.

### Example 1: Basic Usage (`example1.m`)
- Simple dataset from Krippendorff (2011)
- Demonstrates interval, ordinal, nominal, and ratio data analysis
- 4 observers × 12 observations

### Example 2: Bootstrap Confidence Intervals (`example2.m`) 
- Dataset from Hayes & Krippendorff (2007)
- Shows how to calculate 95% confidence intervals using bootstrapping
- Includes visualization of bootstrap distributions
- 5 observers × 40 observations

### Example 3: Time Series Analysis (`example3.m`)
- Synthetic sine wave data with noise
- Demonstrates reliability analysis for continuous time series
- Compares exact, N=2 fast, and alpha-prime methods
- Includes plotting and statistical assessment

## Content
The main function is a wrapper that calls the different functions to calculate alpha and the bootstrapping procedure. For more documentation and options, see the individual functions:
* `kripAlpha.m` : Krippendorff's Alpha for *interval*, *ordinal*, *nominal*, *ratio*, or *phase* data using the exact calculation of Alpha.
* `alphaprime.m` : Approximation of Krippendorff's Alpha for large datasets with arbitrary numerical precision based on binning the data.
* `kripAlphaN2fast.m` : Fast, exact calculation of Krippendorff's Alpha for interval or ordinal data with two observers (*N* = 2).
* `bootstrap_alpha.m` : Run the bootstrapping procedure based on the output from either of the functions above. 

## Toolbox Contents

The main function (`reliability_analysis.m`) is a wrapper that delegates computation to specialized functions:
- `kripAlpha.m`: Exact calculation of Krippendorff's Alpha for interval, ordinal, nominal, ratio, or phase data.
- `alphaprime.m`: Approximate Alpha for large datasets with arbitrary numerical precision using data binning.
- `kripAlpha2fast.m`: Fast, exact calculation for interval or ordinal data with two observers (*N* = 2).
- `bootstrap_alpha.m`: Bootstrapping procedure for estimating confidence intervals.

See the documentation in each function and the companion paper for further details.

## Troubleshooting

### Common Issues

**"Error: input data should be a 2-dimensional NxM matrix"**
- Ensure your data is a 2D matrix with observers as rows and observations as columns
- Use `size(your_data)` to check dimensions

**Low or negative Alpha values**
- This may indicate poor reliability between observers
- Alpha values range from -1 to 1, where 1 = perfect agreement
- Values ≥ 0.8 are typically considered reliable
- Values < 0.67 suggest questionable reliability
- Negative values indicate systematic disagreement

**Memory issues with large datasets**
- Use `'alphaprime'` method for very large datasets
- Consider using `'N2fast'` for two-observer scenarios

**NaN values in data**
- Missing values are handled automatically
- Ensure NaN values represent missing observations, not measurement errors

### Performance Tips

- Use `'N2fast'` methods when you have exactly 2 observers
- Use `'alphaprime'` for datasets with millions of observations
- Bootstrap calculations can be computationally intensive - start with smaller bootstrap samples (e.g., 1000) for testing

## Background

## Contact and Contribution

The code is continuously maintained and kept up to date. If you have suggestions for improvements or additional features, you are welcome to contact me by email or open a GitHub issue. All feedback and contributions are welcome.

### Ways to Contribute
- 🐛 Report bugs via [GitHub Issues](https://github.com/mcvinding/reliability_analysis/issues)
- 💡 Suggest new features or improvements 
- 📚 Improve documentation
- 🧪 Add test cases or examples
- 🔧 Submit pull requests with bug fixes or enhancements

### Contact
**Email:** mvi@psy.ku.dk  
**GitHub:** [@mcvinding](https://github.com/mcvinding)

## Read more
For more information on Krippendorff's Alpha and its application, see:
- Hayes, A. F. & Krippendorff, K. (2007). Answering the call for a standard reliability measure for coding data. *Communication Methods and Measures*, 1, 77–89.
- Krippendorff, K. (2018). *Content Analysis: An Introduction to Its Methodology* (Fourth Edition). SAGE.
- Vinding, M. C. (2025). A Unified Framework for Reliability Analysis in Neuroimaging With Krippendorff’s α. *International Journal of Imaging Systems and Technology, 35(5)*, e70192. https://doi.org/10.1002/ima.70192
