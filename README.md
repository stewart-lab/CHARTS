# CHARTS: Characterizing Tumor Subpopulations

CHARTS is a web application, and associated backened analysis pipeline, for exploring publicly available single-cell RNA-seq data from human tumors. 

## Running the CHARTS web application locally

### Install dependencies

The web application's dependencies are outlined in `charts_web/requirements.txt`.  To install these dependencies in a Python virtual environment run the following commands:

```
python3 -m venv charts_web_env 
source charts_web_env/bin/activate
pip install -r charts_web/requirements.txt  
```

### Run the web application

To run the server locally, perform the following steps:
1. Download the associated backend database. This database is [available to download from Box](https://uwmadison.box.com/s/e2nnhzwgiuww4uid199rshsjuysf04bx).
2. Unpack the database with the following command: `tar -zxf charts_db.tar.gz`
3. Set the absolute path of the `charts_db` directory in `charts_web/config.json`
4. Run the server: `python charts_web/index.py` 

## Running the pipeline on your own single-cell data

Here we document the steps required to run the backend pipeline on your own data and add the results to the CHARTS database for exploring via a local instance of the CHARTS web application. To provide an example, we'll demonstrate how to run the CHARTS pipeline on the data in 
[example/GSE70630_MGH36_MGH53.tsv.gz](https://github.com/stewart-lab/CHARTS/blob/master/example/GSE70630_MGH36_MGH53.tsv.gz). This dataset comprises cells from two tumors in [GSE70630](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70630) by [Tirosh et al. (2016)](https://pubmed.ncbi.nlm.nih.gov/27806376/) (Specifically, the supplemental file "GSE70630_OG_processed_data_v2"). See the Jupyter notebook [example/create_toy_data_GSE70630.ipynb](https://github.com/stewart-lab/CHARTS/blob/master/example/create_toy_data_GSE70630.ipynb) to see how this dataset was generated. 

### Install dependencies

The backend pipeline uses a mix of both Python and R (v3.6.2) and thus, packages for both languages are required. 

Most of the Python dependencies are outlined in `charts_pipeline/requirements.txt`.  To install these dependencies in a Python virtual environment run the following commands:

```
python3 -m venv charts_pipeline_env
source charts_pipeline_env/bin/activate
pip install -r charts_pipeline/requirements.txt
```

In addition, the pipeline requires that [CellO](https://github.com/deweylab/CellO) be installed. To intall CellO, run: 
```
bash charts_pipeline/install_cello.sh
```

The  R dependencies are as follows:
* [GSVA (v1.34.0)](https://www.bioconductor.org/packages/release/bioc/html/GSVA.html)


### Running the pipeline on your custom data

There are two main steps required to run the pipeline on your own data.  The first, involves storing the raw data into the backend HDF5 database, located in `charts_db/charts.h5`.  The second involves running the computational pipeline on this raw data.

#### 1. Set pipeline parameters

The file [charts_pipeline/config](https://github.com/stewart-lab/CHARTS/blob/master/charts_pipeline/config.json) points the pipeline to the database and other required resources. Update these fields to point to the `charts_db` directory, a location to store temporary artifacts from the pipeline, and a directory to store the pipeline's log files.

#### 2. Load expression data into the HDF5 database

Each tumor's expression matrix is stored separately in `charts_db/charts.h5`.  These expression matrices must be in log transcriptipts per million (TPM). Specifically, log(TPM+1) where log is the natural logarithm.  For each tumor, there are three datasets that must be instantiated:
* `/per_tumor/<TUMORNAME>/log1_tpm`: An NxM matrix of log(TPM+1) expression values where N are the number of cells and M are the number of genes. 
* `/per_tumor/<TUMORNAME>/cell`: An N length vector of strings encoding each cell's ID
* `/per_tumor/<TUMORNAME>/gene_name`: An M length vector of strings encoding each gene's symbol

These datasets may be populated using your favorite programming language. We use Python's [h5py](http://www.h5py.org) package for dealing with HDF5 files. As an example, see [example/add_expression_data_to_charts_db.ipynb](https://github.com/stewart-lab/CHARTS/blob/master/example/add_expression_data_to_charts_db.ipynb) where we add the [example/GSE70630_MGH36_MGH53.tsv.gz](https://github.com/stewart-lab/CHARTS/blob/master/example/GSE70630_MGH36_MGH53.tsv.gz) dataset to the CHARTS HDF5 file. 

#### 3. Enter metadata for new tumors

For each new tumor added to the dataset, update the tumor metadata file located at `charts_db/tumor_metadata.json` to include information describing this tumor. This JSON file maps each tumor name to its characteristics:

```
{
        "pub_url": "",
        "pub_name": "",
        "age": "",
        "sex": "",
        "cancer_type": "",
        "cancer_type_abbrev": "",
        "grade": "",
        "stage": "",
        "genomic_alteration": "",
        "lesion_type": ""
    }
```

Note, the pipeline will run fine if a field is blank.  For the example data in [example/GSE70630_MGH36_MGH53.tsv.gz](https://github.com/stewart-lab/CHARTS/blob/master/example/GSE70630_MGH36_MGH53.tsv.gz), see [example/GSE70630_MGH36_MGH53_metadata.json](https://github.com/stewart-lab/CHARTS/blob/master/example/GSE70630_MGH36_MGH53_metadata.json) for example metadata to add to `charts_db/tumor_metadata.json`.

#### 4. Define tumor-specific pipeline parameters

The file [charts_pipeline/tumor_parameters.json](https://github.com/stewart-lab/CHARTS/blob/master/charts_pipeline/tumor_parameters.json) stores tumor-specific parameters. There is one particular piece of data that must be added to this file for the CHARTS pipeline to run: the group of tumors to consider when calculating the malignancy scores.  Specifically, the malignancy score calculation looks at multiple tumors from a single study and attempts to find cells whose genomic alterations are unique to its tumor. The idea here is that each individual tumor is associated with a relatively unique copy-number profile and thus, malignant cells' copy number profiles should cluster with only cells from the same tumor.  Thus, we must supply the CHARTS pipeline with the group of tumors that will be considered in unison. Each group is stored in the `study_to_tumors` field. In this field, each key is the name of a study and the values are all of the datasets included in that study that will be considered in unison for detecting malignant cells.

For the example dataset in [example/GSE70630_MGH36_MGH53.tsv.gz](https://github.com/stewart-lab/CHARTS/blob/master/example/GSE70630_MGH36_MGH53.tsv.gz), one would add the following to the `study_to_tumors` field:

```
"my_study": [
    "My_Tumor_1",
    "My_Tumor_2"
]
``` 

#### 5. Running the computational pipeline

Once the data is loaded into the database, and all metadata and parameters are specified, the CHARTS pipeline can then be executed. The CHARTS pipeline is implemented using Snakemake and stored in `charts_pipeline/Snakefile`. To run the Snakemake pipeline, execute the following command:

`snakemake`

The pipeline will automatically detect which tumors have yet to be processed in the `charts.h5` file and will process those tumors. Processing each tumor entails the following steps:

1. Perform dimensionality reduction on each tumor using both [UMAP](https://arxiv.org/abs/1802.03426) and [PHATE](https://doi.org/10.1038/s41587-019-0336-3) in two and three dimensions.
2. Cluster each tumor. Each tumor is clustered at three levels of granularity. 
3. Classify the cell type of each cluster. Classification is performed for clusters at all three cluster granularities.
4. Compute the gene set enrichment of each cluster using [GSVA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7). GSVA is performed on all clusters at all three cluster granularities.
5. Compute the malignant score of each cell. 
6. Compute differentially expressed (DE) genes between each cluster and the remaining clusters in the tumor. DE genes are computed for all clusters at all three cluster granularities.

All results are written back into the `charts.h5` file. Note, that we can tell the pipeline to re-process *all* of the tumors in the charts.h5 if we set the `OVERWRITE` variable to True at the top of the Snakefile. Reprocessing all of the tumors may take several hours.
  
### Deleting a tumor from the database

To delete a tumor from the CHARTS database, one can use the `delete_tumor.py` script. Specifically, to delete a tumor, run:

`python delete_tumor.py <charts.h5 location> <tumor name>`




 
