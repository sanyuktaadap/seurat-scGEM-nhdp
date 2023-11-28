## scRNAseq-scgem-nhdp

### Libraries
- Seurat
- magrittr
- dplyr
- tibble
- msigdbr
- harmony
- biomaRt
- clusterProfiler
- org.Hs.eg.db
- lsa
##### Run [run_pbmc3k](./run_pbmc3k.R) to understand the Seurat workflow on sample PBMC data

### Training
##### Use [utils](./utils/) to [split data](.//utils/split_data.R) (for easy computing only), add [cell annotations](./utils/save_exp_batch_by_type.R), [merge the final data](./utils/merge_cell_type_data.R) together, and [split the data according to cell type](./utils/cell_annotation_utils.R). Also use essential functions from [utils](./utils/utils.R).
##### Run [training](./run_train.R) to train the [nHDP](./model_windows.R) model on  single-cell data

### Running Inference
##### Use a different single-cell dataset to infer with the trained model
##### Script for inference can be found [here](./inference/run_inference.R)
##### Analyse inferred data with [this](./inference/analyse_lc_tnk_inferred_data.R) script
##### Analyse training data with [this](./inference/analyse_crc_tnk_data.R) script
##### Create UMAP plots using Seurat with [this](./inference/tnk_plots.R) script
