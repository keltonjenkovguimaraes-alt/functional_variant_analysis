#Functional Variant Impact Analysis
 pipeline for interpreting the biological consequences of genomic variants. Takes annotated VCF files and performs gene ontology enrichment, pathway analysis, protein domain mapping, variant prioritization, and affector gene screening.
## Analyses

| Analysis | Input | Method | What It Reveals |
|----------|-------|--------|-----------------|
| **GO Enrichment** | Gene list | Hypergeometric test | Overrepresented biological processes |
| **Pathway Analysis** | Mutated genes | KEGG/Reactome mapping | Affected metabolic pathways |
| **Protein Domains** | Missense variants | Pfam/InterPro mapping | Domains affected by mutations |
| **Variant Prioritization** | Annotated VCF | Impact scoring | High-priority variants for follow-up |
| **Effector Gene Analysis** | Secreted proteins | SignalP detection | Pathogenicity factor mutations |
## Quick Start

```bash
git clone https://github.com/keltonjenkovguimaraes-alt/functional_variant_analysis.git
cd functional_variant_analysis
conda env create -f workflow/envs/environment.yaml
conda activate func_analysis
python workflow/scripts/run_all_analyses.py
Input
File	Description
Annotated VCF	VEP-annotated variant calls
Reference GFF	Gene annotations
Protein FASTA	Predicted proteome
Author
Kelton Guimaraes — Implementation & Analysis
