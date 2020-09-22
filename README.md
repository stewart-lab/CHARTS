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


