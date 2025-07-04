The following code was used to generate the outputs presented in Figure 2 and some of the supplementary figures of the manuscript "*In Situ Inference of Copy Number Variations in Image-Based Spatial Transcriptomics*".

### Figure 2: Evaluation of technical constraints affecting CNV inference.
<img src="https://github.com/user-attachments/assets/88ed1195-6745-4d19-b175-23c38715ca3e" width="600">

**a**, Workflow for evaluating how technical constraints affect CNV inference performance. Three CNV subclones (A, B,
and C) were simulated alongside a normal (N) cell population lacking CNVs. The phylogenetic tree shows the gains
and losses defining each subclone, which are also represented in the CNV heatmap, where gain (red) and losses
(blue) reflect the simulated ground truth. Technical limitations of gene panel sizes and detection efficiency were
applied to generate 63 unique datasets, five of which are shown in b. **b**, CNV heatmaps (left), shows inference results
of the five selected datasets, two of which illustrate reduced detection efficiency, and two illustrate reduced gene
panels. UMAPs (right) showing distributions of cells based on their inferred CNV profiles and compares their
predicted vs true subclone identities. **c**, Performance metrics for CNV gain and loss prediction, as well as subclone
identification. Black-circled data points highlight the five datasets shown in panel b.
