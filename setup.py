import setuptools

# This setup script lets you install the Python part of combine like a standard
# Python package. The big advantage of it is you can make an in-place
# installation into your virtual environment, which means you can change the
# Python files in this repo, and the changes will be effective immediately
# without any recompilation or reinstallation.

# Before installing the package, run in this repo:
#   mkdir -p .python/HiggsAnalysis
#   ln -s $PWD/python $PWD/.python/HiggsAnalysis/CombinedLimit
#   touch $PWD/.python/HiggsAnalysis/__init__.py
#   touch $PWD/.python/HiggsAnalysis/CombinedLimit/__init__.py

# To do an in-place installation, it's best to create and activate a virtual
# environment, if you don't have done so already:
#   python3 -m venv <directory>
#   source <directory>/bin/activate

# Finally, the in-place installation (for regular installation, drop the "-e"):
#   python3 -m pip install -e .

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

scripts = [
    "combineCards.py",
    "plotGof.py",
    "pruneUncerts.py",
    "combineTool.py",
    "plotImpacts.py",
    "text2workspace.py",
    "commentUncerts.py",
    "plot1DScan.py",
    "plotLimitGrid.py",
    "plotBSMxsBRLimit.py",
    "plotLimits.py",
]

setuptools.setup(
    name="HiggsAnalysisCombinedLimit",
    version="1.0",
    description="CMS Combine",
    author="The CMS Collaboration",
    long_description=long_description,
    long_description_content_type="text/markdown",
    package_dir={"": ".python"},
    packages=setuptools.find_packages(where=".python"),
    python_requires=">=3.8",
    scripts=["scripts/" + filename for filename in scripts],
)
