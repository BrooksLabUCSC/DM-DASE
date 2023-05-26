
### JuncBASE
wget http://compbio.berkeley.edu/proj/juncbase/juncBASE_v0.6.tgz
tar zxvf juncBASE_v0.6.tgz
cd /juncBASE_v0.6

### MESA
git clone https://github.com/BrooksLabUCSC/mesa.git
cd mesa/
pip install .

### DRIMSeq
packages.install("devtools")
devtools::install_github("markrobinsonuzh/DRIMSeq")

### GSEApy
conda install -c conda-forge -c bioconda gseapy 
# OR
pip install gseapy