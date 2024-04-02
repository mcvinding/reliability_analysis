# reliability_analysis
Assessing reliability between measurements in neuroscience.

The intention is to provide a test for relaiability tests of continous neural data and neural time series using Krippendorff's Alpha or an approximation that is analogous to Krippendorff's Alpha. All code is implemented in MATLAB.

Implements three different methods, see documentation and options for the individual functions.
* `kripAlpha.m` : Krippendorff's Alpha for interval, ordinal, nominal or phase data using the exact caluclation.
* `alphaprime.m` : Aproximation of Krippendorff's Alpha for interval data with arbitary precision
* `kripAlpha2fast.m` : Fast calculation of Krippendorff's Alpha for interval data with two oberservers and no missing cases.

Use the function 'reliability_analysis' to call the approbiate type of test.

For more information on Krippendorff's Alpha see: 
* Hayes, Andrew F. & Krippendorff, Klaus (2007). Answering the call for a standard reliability measure for coding data. *Communication Methods and Measures*, 1, 77-89

So far it is only work in progress and not tested if it works or is compatible.
