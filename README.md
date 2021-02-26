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
git fetch origin pull/648/head:r622conda # until merged
git checkout r622conda

conda install --name base mamba # faster conda
mamba env create -f conda_env.yml 
conda activate combine2
source env_conda.sh 
make -f Makefile_conda -j 8
```

