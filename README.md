# Functionality for quality evaluation for direct / indirect adjustment

## Aim

To provide a package for testing and quality evaluation for the indirect SA.

## Steps

1. Functions for quality criteria calculation for two time series – directly and indirectly adjusted (tsda and tsia). Could be the seasonally adjusted series, the trend-cycle estimates and the seasonal components.
    - Input: two sets of time series (two ts objects with five time series – y, sa, t, s, i).
    - Parameters:
        - Length of time series in years (could be decimal): full (NULL), defined time period (for example last 3 years).
    - Output:
        - named list of values for quality criteria.
        - plots
2. Functions for running tests on spreadsheets of time series (tssda and tssia).
    - Input: two tables with time series.
    - Output:
        - table (or list) with quality criteria.
        - plots
3. Functions for doing direct (with RJD) and indirect adjustment (with KIX) and computing quality criteria for all time series.
    - Input: time series, weights and aggregation specification (how to define in R? formulas?).
    - Output:
        - table (or list) with quality criteria.
        - plots
