# Dragonfly
A collection of scripts for demographic history and inbreeding analyses of green hawker (*Aeshna viridis*)

For each major analysis, one file (or more, in which case they are numbered) named run_<analysis_name>.sh describes the overall procedure, including when to use any other smaller scripts or "helper" files.
Note that despite the file format, the run_<analysis_name>.sh files are not a bash scripts that can be run in their entirety. Rather, they contain multiple small pieces of code which should be run independently.

For example, in the folder "1. Mapping", there are two files:
- run_mapping.sh, which describes the overall procedure for mapping of reads.
- STATS.r, a small "helper script", which is used in the analysis detailed in the run_mapping.sh file to obtain mapping statistics.

Note that all run_<analysis_name>.sh-files contain a function to check job info in batch, which can be used at will, when wanting the check the status of multiple jobs. 
