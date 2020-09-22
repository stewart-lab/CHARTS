git clone https://github.com/deweylab/CellO.git
cwd=$(pwd)
cd CellO
bash download_resources.sh 
export PYTHONPATH=$(pwd):$PYTHONPATH`
cd $(cwd)
