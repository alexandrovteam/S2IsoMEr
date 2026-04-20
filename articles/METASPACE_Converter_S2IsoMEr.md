# METASPACE to S2IsoMEr

In this vignette, we’ll demonstrate how to use the [METASPACE
converter](https://github.com/metaspace2020/metaspace-converter) to
transform [METASPACE](https://metaspace2020.org/) datasets into a format
compatible with S2IsoMEr. Note that these datasets will be in a
single-pixel format rather than single-cell unless an additional mapping
from a cell-segmentation output links pixels to cells.

Since the METASPACE converter is written in Python, we’ll provide a
snippet demonstrating how to first export a METASPACE dataset as an
[AnnData](https://anndata.readthedocs.io/en/stable/) object and then
import it into R using the
[anndataR](https://github.com/scverse/anndataR) package. This vignette
provides only a basic overview of the METASPACE converter. For
installation instructions and more detailed usage, please refer to the
[documentation](https://metaspace2020.github.io/metaspace-converter/).

We’ll use the METASPACE datasets from the same NASH model referenced in
the [ORA
vignette](https://alexandrovteam.github.io/S2IsoMEr/articles/ORA_enrichment-single_cell_metabolomics.md).
For more details about the dataset, refer to that vignette.

## Loading packages

``` r
require(S2IsoMEr)
require(dplyr)
require(tidyr)
require(Seurat)
require(anndataR)
```

## Getting Dataset IDs

As in the [ORA
vignette](https://alexandrovteam.github.io/S2IsoMEr/articles/ORA_enrichment-single_cell_metabolomics.md),
we’ll start by downloading the datasets and mapping the condition to the
dataset ID.

``` r
NASH_scm_tmp = tempfile()
download.file("https://zenodo.org/records/13318721/files/NASH_scm_dataset.rds", destfile = NASH_scm_tmp)
NASH_scm = readRDS(NASH_scm_tmp)

condition_metadata_tmp = tempfile()
download.file("https://zenodo.org/records/13318721/files/spacem_scm_matrices.rds", destfile = condition_metadata_tmp)
condition_metadata = readRDS(condition_metadata_tmp)[["metaspace_dataset_names"]]
condition_metadata$dataset_name = sub(".ibd.*", "", condition_metadata$dataset_name)

metaspace_annotations_tmp = tempfile()
download.file("https://zenodo.org/records/13318721/files/SpaceM_metaspace_ds_annotations.rds", destfile = metaspace_annotations_tmp)
metaspace_annotations = readRDS(metaspace_annotations_tmp)
```

### Get metadata

``` r
scm = NASH_scm$scm %>% 
  as.matrix() %>% 
  t()

conds = NASH_scm$metadata %>%
  column_to_rownames("Cell")
conds = conds[colnames(scm),]

conds_unique = conds %>% 
  dplyr::distinct()

metaspace_annotations = metaspace_annotations %>% 
  dplyr::left_join(condition_metadata, by = c("ds_name" = "dataset_name")) %>%
  dplyr::rename("Replicate" = "Condition") %>%
  dplyr::left_join(conds_unique) %>% 
  dplyr::select(ds_id, Replicate, Condition) %>% 
  dplyr::distinct()
```

### Filter conditions

``` r
cond_x = "U"
cond_y = "F"

Dss_id_conds = metaspace_annotations %>% 
  dplyr::filter(Condition %in% c(cond_x, cond_y))

#Change "path/to/dir" with your local directory
write.csv(Dss_id_conds, "path/to/dir/Ds_ids_for_METASPACE_converter.csv", row.names = F)
```

## METASPACE to AnnData

Now you can use the exported metadata with dataset IDs and use METASPACE
converter to get the merged AnnData object. You’ll need to run the
following python chunk in your own python environment to export the
AnnData object to an H5AD file, which you’ll then import back into R in
the next step.

``` python
from metaspace_converter import metaspace_to_anndata
import scanpy as sc
import pandas as pd

#Change /path/to/dir with local path
dss = pd.read_csv("/path/to/dir/Ds_ids_for_METASPACE_converter.csv")

dss_list = list(dss['ds_id'])

FDR_level = 0.1
DB = ("HMDB", "v4")
Adatas = []

for ds in dss_list:
    adata = metaspace_to_anndata(
        dataset_id=ds,
        fdr=FDR_level,
        database=DB)
    Replicate = dss['Replicate'][dss['ds_id'] == ds].values[0]
    adata.obs.index = Replicate + "_" + adata.obs.index
    adata.obs['Condition'] = dss['Condition'][dss['ds_id'] == ds].values[0]
    adata.obs['Replicate'] = Replicate
    Adatas.append(adata)

merged = sc.concat(Adatas, axis=0, join='outer', merge='first')
merged.var["offSample"] = merged.var["offSample"].astype(str)

#Change /path/to/dir with local path
merged.write("/path/to/dir/NASH_FU_Anndata_FDR_0.1_HMDB_v4.h5ad")
```

## Read AnnData into R

``` r
adata = anndataR::read_h5ad("/path/to/dir/NASH_FU_Anndata_FDR_0.1_HMDB_v4.h5ad")
seurat_obj = anndataR::to_Seurat(adata)

scm = seurat_obj@assays$RNA@data %>% as.matrix()
grps =  seurat_obj@meta.data$Condition %>% as.character()
```

## Initialize S2IsoMEr object for ORA

Now that we have the input single-pixel matrix and the corresponding
condition for each pixel, we can create the `S2IsoMEr` object for ORA
enrichment as before.

``` r
cond_x = "U"
cond_y = "F"

ORA_boot_obj = initEnrichment(scmatrix = scm, conditions = grps,
                     enrichment_type = "ORA",annot_db = "HMDB",
                     consider_isomers = T, consider_isobars = T,
                     polarization_mode = "positive",
                     background_type = "sub_class",
                     molecule_type = "Metabo",
                     condition.x = cond_x,
                     condition.y = cond_y)
```
