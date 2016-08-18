# MSGARCH
MSGARCH GSOC project
The goal of this project is to implement a package that will give the financial community tools to estimate,
simulate, and test several MSGARCH models used in volatility (i.e., square root of conditional variance) forecasting.
By relying on a hidden/latent variable, these models are able to switch among several processes for the conditional
volatility and therefore, account for structural break in the volatility dynamics. MSGARCH have gained a huge interest 
in the financial risk management community over the recent years as they are better at forecasting volatility and provide
more accurate risk measures. The package will follow the structure of rugarch since this is one of the most used packages
for volatility modeling. The core will be implemented in C++ while simple R functions will facilitate usage of the package.
Currently, no R package is available to estimate these models.
