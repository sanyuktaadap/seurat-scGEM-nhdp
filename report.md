# Finding Functional Similarities in Colorectal (CRC) and Lung Cancer (LC)

#### Reference
Zhang, H., Lu, X., Lu, B., Chen, L., Liu, X., Liu, L., Li, F., & Liu, X. (2023). scGEM: Unveiling the Nested Tree-Structured Gene Co-Expressing Modules in Single Cell Transcriptome Data. Cancers (Basel), 15(17), 4277. doi:10.3390/cancers15174277

#### Data
Transcriptomic scRNA-seq data obtained from [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907)

#### Model
Probabilistic Tree model – nHDP

#### Goals
 - Train the nHDP on CRC data for the model to learn genes associated with CRC in each GEM (a node of the tree model)
 - Infer using LC data to identify LC genes associated with each Gene Expression Module
 - Identify the biological processes involved in the co-expression of a group of genes
 - Identify Transcription Factors and Signaling Pathways involved

#### A High-level Overview of the Workflow

![image](https://github.com/sanyuktaadap/seurat-scGEM-nhdp/assets/126644146/16202c54-6f69-4d4c-80b9-d46df6552791)



_Gene Expression Modules (GEMs):_ Nodes of the tree model. The cells and genes that settle in the node are the ones that are involved in the function the node represents.

- **Part a** : Data that we feed the model (gene-by-cell matrix)
- **Part b** : Logic behind the nHDP
  - Takes the input and segregates it in the tree's nodes.
  - Each node represents a biological function that is carried out due to the genes expressed in the cells settled in that node.
  - The model returns 2 output matrices.
- **Part c** : Output of the model
  - gene-by-gem matrix
  - cell-by-gem matrix
- **Part d** : Identify transcription factors involved in regulating these genes and the highly enriched signaling pathways that involve the co-expression of the genes.

#### Obtaining Results

- The saved model was used to infer LC data.
  - The genes settle in the nodes of the model.
  - Each node represents a particular biological function in CRC.
  - A gene will settle in a node based on the function the node represents.
  - Each node will have a group of genes.
    - These groups of genes are involved in performing that particular function.
- Integrating, analyzing, and plotting CRC and LC data using Seurat
  - Created a Seurat object with merged CRC and LC data.
  - Performed dimensionality reduction using UMAP.
  - UMAP feature plots are created for individual gene expression patterns.
  - These UMAPs are compared with feature plots of various T-cell markers.
- 2 such GEM UMAP plots are shown below.
- Transcription Factor (TF) Analysis was performed using the ChIP-X Enrichment Analysis Tool
- Identified Signalling Pathways (SP) from the KEGG database using Entrez IDs the gene sets.
- The results could be used in designing targeted immune therapy that aims to treat CRC as well as LC, instead of having two separate immune therapies.

1. GEM3
![image](https://github.com/sanyuktaadap/seurat-scGEM-nhdp/assets/126644146/c3197982-157b-4c81-a2a3-36bac24d236d)

- TF - SP140L: DNA-binding transcription factor activity
- SP - Proteasome: Ubiquitination of proteins
- T-cell marker - CTLA-4: Regulator of T-cell homeostasis and self-tolerance

1. GEM23

![image](https://github.com/sanyuktaadap/seurat-scGEM-nhdp/assets/126644146/fb3c8f1e-b818-45fb-acf1-ab3bb825ab1d)

- TF - FOXM1: Proliferation-associated transcription factor activity
- SP - Cell Cycle: Division of Cells
- T-cell marker - MKI67: Regulation of chromosome segregation mitotic nuclear division
