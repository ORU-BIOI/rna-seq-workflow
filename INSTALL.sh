wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-Linux-x86_64.sh
chmod u+x Miniconda-3.7.0-Linux-x86_64.sh
pref=/proj/b2014206/miniconda 
./Miniconda-3.7.0-Linux-x86_64.sh -b -p $pref
export PATH="$pref/bin:$PATH"
echo $PATH
conda create -n sci2 python=2
conda create -n sci3 python=3

source activate sci2
conda install pip argcomplete numpy scipy scikit-learn pandas ipython-notebook matplotlib binstar
pip install ipdb
pip install cutadapt
pip install htseq
source deactivate

source activate sci3
conda install pip argcomplete numpy scipy scikit-learn pandas ipython-notebook matplotlib binstar
pip install ipdb
pip install snakemake




