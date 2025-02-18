# CoMeBackV2
Co-methylation of adjacent array probes, taking into account genomic CpG background. Function cmr() constructs co-methylated regions based on pair-wise positive correlations above variable correlation threshold.  The threshold increases linearly from a low value for low-density CpG background to a high value for high density background. Modified for the EPICv2 array.

https://github.com/Treys925/CoMeBackV2_R_package

# To Install
library(devtools) <br />
install_github("Treys925/CoMeBackV2_R_package")
