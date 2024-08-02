# InSituCNV

Temporary README for project structure: 

## Figure 1 - Real Data (show that it works)

Data types:
- CosMX data
- Xenium data


## Figure 2 - Simulated Data (test technical limitations)
Data types:
* Xenium 5K (human lymph node) from 10X 

Workflow:
- 01 - Simulate CNVs in spatial data
- 02 - Generate metrics for prediction accuracy
  - Malignant cell classification
  - CNV identification
  - Subclone identification
- 03 - Downsample simulated data (apply technical limitations)
  - gene coverage
  - read depth
- 04 - Run inferred CNV analysis
  - SCEVAN
  - inferCNV
- 05 - Apply prediction metrics to resulting iCNV analysis data
- 06 - Generate illustrative figures
