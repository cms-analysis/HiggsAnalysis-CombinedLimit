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
```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git checkout 112x

conda install --name base mamba # faster conda
mamba env create -f conda_env.yml

conda activate combine
bash set_conda_env_vars.sh
# Need to reactivate
conda deactivate
conda activate combine

make -f Makefile_conda -j 8
```

Using combine from then on should only require sourcing the conda environment 
```
conda activate combine
```
