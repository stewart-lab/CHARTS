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

Here we document the steps required to run the backend pipeline on your own data and add the results to the CHARTS database for exploring via a local instance of the CHARTS web application.

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

#### Loading scRNA-seq data into the database

Each tumor's expression matrix is stored separately in `charts_db/charts.h5`.  These expression matrices must be in log transcriptipts per million (TPM). Specifically, log(TPM+1) where log is the natural logarithm.  For each tumor, there are three datasets that must be instantiated:
* `/per_tumor/<TUMORNAME>/log1_tpm`: An NxM matrix of log(TPM+1) expression values where N are the number of cells and M are the number of genes. 
* `/per_tumor/<TUMORNAME>/cell`: An N length vector of strings encoding each cell's ID
* `/per_tumor/<TUMORNAME>/gene_name`: An M length vector of strings encoding each gene's symbol

These datasets may be populated using your favorite programming language. We use Python's [h5py](http://www.h5py.org) package for dealing with HDF5 files. 

#### Add the metadata for new tumors

For each new tumor added to the dataset, update the tumor metadata file located at `charts_db/tumor_metadata.json` to include information describing this tumor. This JSON file maps each tumor name to its characteristics.  Currently, the following entries are required for each tumor:
* `cancer_type`: An english description of the cancer type (e.g. melanoma).
* `cancer_type_abbrev`: The TCGA abbreviation for the cancer type (from [https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations)).

#### Running the computational pipeline

Once the data is loaded in the database, the CHARTS pipeline can then be executed. The CHARTS pipeline is implemented via the Snakemake pipeline specified in `charts_pipeline/Snakefile`.


 
