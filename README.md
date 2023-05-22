# bi-2023-tcr-tils
repository with the results on the semester project 

# Re-analysis of the public dataset: comparison and evaluation of T cells in paired tumor and normal lung single cell samples

---
## Motivation
Due to numerous mutations in the tumor genome, the neoplasm expresses mutation-associated neoantigens (MANA) that can trigger an immune response against cancer in patients. However, as the tumor evolves, malignant cells often acquire the ability to evade the immune response by expressing checkpoint proteins. Immunotherapy using antibodies against these checkpoints has been developed, but a significant number of patients do not respond to this treatment. Biomarkers, such as the expression of checkpoint proteins, tumor mutational load, and neoantigenic load, are used to determine the suitability of immunotherapy. These biomarkers are particularly relevant for highly mutated cancers like melanoma and lung cancer. However, they are not comprehensive, highlighting the need for further research on the immune component of tumors and also dysfunctional programs in MANA specific tumor-infiltrating lymphocytes (TILs).

Our project is based on open data from the article [“Transcriptional programs of neoantigen-specific TIL in anti-PD-1-treated lung cancers”](https://www.nature.com/articles/s41586-021-03752-4).


## Methods and Results

The notebook with the full analysis can be found in this repo. Here I present and summarize the main results of the work. A short recap can also be found in the presentation. 

![Pipeline of analysis](figs/pipeline.svg)

### QC and filtration

Low-quality cells were filtered out if:

> the number of detected genes was below 250 or above 2500
* the proportion of mitochondrial gene counts was higher than 10%
* the proportion of ribosomal gene counts was lower than 10%

doublets were removed with with the help of Scrublet 

Mitochondrial genes, and genes associated with poorly supported transcriptional patterns were also removed from the analysis


