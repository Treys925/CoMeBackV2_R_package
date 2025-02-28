# CoMeBackV2
Co-methylation of adjacent array probes, taking into account genomic CpG background. Function cmr() constructs co-methylated regions based on pair-wise positive correlations above variable correlation threshold.  The threshold increases linearly from a low value for low-density CpG background to a high value for high density background. Modified for the EPICv2 array and Homo sapiens (human) genome assembly GRCh38 (hg38) from Genome Reference Consortium.<br />
Now adapted for genome assembly HS1 from T2T Consortium CHM13v2.0 <br />
[https://www.science.org/doi/10.1126/science.abj6987](https://www.science.org/doi/10.1126/science.abj6987)
<br />
For EPICv2 please remove duplicate probes and use 10 character names (e.g. cg01534423)<br />
For build HS1 please load all files from /data into the global environment<br />

[https://github.com/Treys925/CoMeBackV2_R_package](https://github.com/Treys925/CoMeBackV2_R_package)

# To Install
library(devtools) <br />
install_github("Treys925/CoMeBackV2_R_package")

# CoMeBack:
Evan Gatev, Nicole Gladish, Sara Mostafavi, Michael S Kobor, CoMeBack: DNA methylation array data analysis for co-methylated regions, Bioinformatics, Volume 36, Issue 9, May 2020, Pages 2675–2683, [https://doi.org/10.1093/bioinformatics/btaa049](https://academic.oup.com/bioinformatics/article/36/9/2675/5716323)
