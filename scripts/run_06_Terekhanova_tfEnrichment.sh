#!/bin/bash
#SBATCH --job-name=run_06_Terekhanova_tfEnrichment
#SBATCH --output=run_06_Terekhanova_tfEnrichment_%j.out
#SBATCH --error=run_06_Terekhanova_tfEnrichment_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=01:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.2.1

# Run R script
Rscript - <<'EOF'

# load libraries
library(here)
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# set genome to hg38
addArchRGenome("hg38")

# Load the project
proj <- loadArchRProject(path = "Terekhanova_hypoxia")

# Add TF Annotations from ArchRAnno database: EncodeTFBS
## note: the Anno file was downloaded seperately
print("add EncodeTFBS TF annotations from ArchR-Hg38-v1.Anno")
proj <- addArchRAnnotations(ArchRProj = proj,
						db = "/cluster/projects/wouterslab/ArchR/Terekhanova_hypoxia/anno/ArchR-Hg38-v1.Anno",
						collection = "EncodeTFBS" # specific collection to collect
						)
						
# Load Marker Peaks SE object(previously computed and saved)
## eg. markersPeaks_DAEG_UP_24 from run_04_terekhanova_peakCalling.sh
print("load markers peaks")
load(file = "Terekhanova_hypoxia/PeakCalls/markersPeaks_DAEG_UP_24.RData")

# Calculate Enrichment of EncodeTFBS with marker Peak set
print("calculate enrichment of EncodeTFBS via peakAnnoEnrichment")
enrichEncode <- peakAnnoEnrichment(
    seMarker = markersPeaks_DAEG_UP_24,
    ArchRProj = proj,
    peakAnnotation = "EncodeTFBS",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

# Save SE object results of peakAnnoEnrichment()
print("save results of peakAnnoEnrichment")
save(enrichEncode, file = "Terekhanova_hypoxia/motifs/enrichEncode.RData")
						
EOF

echo "TFBS Enrichment analysis completed"