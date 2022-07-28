HiggsAnalysis-CombinedLimit
===========================

### Official documentation

[Manual to run combine](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/)

### Standalone compilation in `lxplus`
```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
source env_standalone.sh 
make -j 8; make # second make fixes compilation error of first
```

### Standalone compilation with `conda`
This recipe will work both for linux and MacOS
```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git checkout 112x

conda install --name base mamba # faster conda
mamba env create -f conda_env.yml

conda activate combine
source set_conda_env_vars.sh
# Need to reactivate
conda deactivate
conda activate combine

make CONDA=1 -j 8
```

Using combine from then on should only require sourcing the conda environment 
```
conda activate combine
```

**Note:** on OS X, `combine` can only accept workspaces, so run `text2workspace.py` first.
This is due to some ridiculous issue with child processes and `LD_LIBRARY_PATH` (see note in Makefile)
