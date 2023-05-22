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

* the number of detected genes was below 250 or above 2500
* the proportion of mitochondrial gene counts was higher than 10%
* the proportion of ribosomal gene counts was lower than 10%

Doublets were removed with with the help of `Scrublet`

Mitochondrial genes, and genes associated with poorly supported transcriptional patterns were also removed from the analysis

### Clustering and annotation

TCR, immunoglobulin and mitochondrial genes, as well as features that constitute Interferon I mediated pathway,  were excluded from clustering  to make sure that clustering result will not be influenced by their variability.

PCA was performed based on the 3,000 most variable genes. UMAP on PCA results have shown that cells group by samples of their origin, indicating the need for a batch effect correction. 

![UMAP of expression data before batch correction](figs/umap_before_harmony.svg)

After harmonization, cell clusters are more evenly distributed among patients, however, the biological difference of the tumor/normal immune environment is not lost

![UMAP of expression data after batch correction](figs/umap_after_harmony.svg)

Leiden clustering resulted in 14 separate clusters, that were annotated using combination of general CD4/CD8 markers with common subset specific markers:

* *FOXP3* for Tregs; 
* *MKI67* for proliferating cells;  
* *CXCL13* for Follicular helpers; 
* *GZMA*, *GZMB*, *GZMK* for effector cells; 
* *ZNF683* and *ITGAE* for memory cells; 
* *KLRC1* for NK cells; 
* *SLC4A10* for MAIT cells

![UMAP of clustering results](figs/umap_clustering.svg)
![Cell subtype marker expression](figs/clustering_markers.svg)

Differentially expressed genes were found using wilcoxon test for each cell type vs all other cells
 

![Dotplot of cell markers, differentially expressed genes and T cell checkpoint associated genes](figs/dotplot.svg)


As expected, cellular distribution in tumor and normal sample is highly different. 

* CD8+ effector2 cells are enriched in Tumor, while  CD8+ effector1 – in Normal tissue
* Both CD4+ helper subtypes are enriched in Normal, although  CD4+ Follicular helpers and CD4+ T reg are more specific for Tumor
* CD8+ mem1 is enriched in Tumor


![Cell type composition by type of tissue](figs/clustering_by_tissue.svg)


