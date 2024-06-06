# Contributing

Contributions to <span style="font-variant:small-caps;">Combine</span> of all sizes, from minor documentation updates to big code improvements, are welcome and encouraged.

To ensure good development of the tool, we try to coordinate contributions.
However, we are happy to help overcome any steps that may pose issues for contributors.

Guidelines for contributions are below.

In general you should always perform the following steps when contributing.

1. Make your changes;
2. Test your changes to ensure they work and don't break anything;
3. If this changes any user-facing aspect, ensure that the appropriate documentation is updated;
4. Check code style (see below), and fix any issues;
5. Create your pull request;
6. Iterate with the maintainers until the request is merged.

## Code Style


We ask that you put some effort into the readability of your code.
We are, however, always happy to help if there is an issue.

We use linting as part of our ci/cd for the python code. 
That means your code will be checked automatically, and you must make sure it conforms to certain rules.

Currently no linting or automatic checks of C++ code are implemented.
Although we do not have a well-defined style guide for C++, we always appreciate readable and well-formatted code.

### Technical details on linting

we use `flake8` and `black` for linting. 
To run the linting locally before making your pull request, or before making a commit, you can do the following.

ensure `flake8` and `black` are installed: 

```
 python -m pip install -q flake8 black
```

and then from the main directory of this repository run


`flake8`:
```
flake8 .
```

and `black`:

```
black -l 160 --check --diff .
```

If you'd like to see the details of the configuration `flake8` is using, check the `.flake8` file in the main directory.
The `black` linting uses the default [black style](https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html) (for v23.3.0), with only the command line options shown above.


## Updating Documentation

It is crucial to our user base and developers that the documentation is well-maintained.
For that reason, whenever you make a change you should consider whether this requires a corresponding documentation update.

If the change is user-facing it almost certainly does require a documentation update.

Documentation is **very important** to us. 
Therefore, we will be meticulous and make sure it is done well!
However, we don't want to put extra burden on you, so we are happy to help and will make our own edits and updates to improve the documentation of your change.

We appreciate you putting in some effort and thought to ensure:

- your documentation is understandable to the audience being targeted;
- your documentation is clear and concise; and
- your documentation fits in properly with the overall documentation structure.

### Technical details of the documentation

We use [mkdocs](www.mkdocs.org) to produce the static website that documents <span style="font-variant:small-caps;">Combine</span>.

The documentation files are all under the `docs/` folder.
Which pages get included in the site, and other configuration details are set in the `mkdocs.yml` file.

In order to check the documentation rendering (features such as latex math rendering, etc) locally, you can generate the site on your local computer and check it in your browser.
To do so, after [installing mkdocs](https://www.mkdocs.org/getting-started/), [pymdown extensions](https://facelessuser.github.io/pymdown-extensions/installation/), the [material theme](https://squidfunk.github.io/mkdocs-material/), and [mike](https://github.com/jimporter/mike) (for versioning) via:

```
python -m pip install mkdocs pymdown-extensions mkdocs-cinder mike
```

 you can do:

```
mkdocs build
mkdocs serve
```

from the main repository directory. mkdocs will then print a link you can open to check the page generated in your browser.

This will serve the current version of the documentation. For the versioned documentation, where the documentation corresponds to a certain tag of the code, mike is used. Mike performs pushes to the `gh-pages` branch of the repository, and sets up certain commits to be accessible as the tagged version of the docs. For adding a new tagged version to the documentation you need to do:

```
mike deploy <version> <alias>
```

Note that this will only affect your *local* copy of the gh-pages documentation. In order to view the pages locally you can use `mike serve` instead of `mkdocs serve`.

**NOTE:** mkdocs builds that use internal links (or images, etc.) with absolute paths will work for local deployment, but will break when deployed to the public documentations pages. 
Please ensure you use relative paths. Currently, this is the only known feature where the behvaiour differs between local mkdocs and public page deployment. 
If you'd like to test the deployment directly, the suggested method is to set up a docs page using your personal github account; this should mimic the exact settings of the official page.

## Github Actions 

There are several automated processes defined in the [`.github/workflows`](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/tree/main/.github/workflows) folder which are triggered for every pull request or when a new commit is pushed into the `main` branch. 

 - `CI with CVMFS` workflow defined in [.github/workflows/cvmfs-ci.yml](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/.github/workflows/cvmfs-ci.yml) compiles <span style="font-variant:small-caps;">Combine</span> with the CMSSW which is setup from `cvmfs`. This will also run several test <span style="font-variant:small-caps;">Combine</span> commands using tutorial datacards contained in the repository. 
 - `CI/CD` workflow defined in [.github/workflows/ci.yml](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/.github/workflows/ci.yml) runs:
   
    1. [Linting checks](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/.github/workflows/ci.yml#L11) with`flake8` and `black`.
    2. <span style="font-variant:small-caps;">Combine</span> compilation using conda environment with several combinations of (ROOT, python) versions:  (6.26.4, 3.10), (6.24,3.9), (6.22, 3.9) and (6.22, 2.7.18) 
    
 - `Docs` workflow defined in [.github/workflows/docs.yml](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/.github/workflows/docs.yml) manages the documentation deployment with `mike` as described above.

In addition there are two conditional workflows that can be triggered with [github labels](https://docs.github.com/en/issues/using-labels-and-milestones-to-track-work/managing-labels): 

- `trigger gitlab job` workflow defined in [.github/workflows/gitlabci.yml](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/.github/workflows/gitlabci.yml) can be triggered with the "safe to test" label. Please note that this gitlab workflow uses `CMSSW_10_2_13`, therefore it might fail with later versions of <span style="font-variant:small-caps;">Combine</span> and it is likely to be deprecated soon. 

- `Port to <branch>` workflows defined in [.github/workflows/port_to_branch.yml](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/.github/workflows/port_to_branch.yml), e.g. the `Port to 112x-comb2022` can be triggered with the "port to 112x-comb2022" label to port the changes from a pull request to the `112x-comb2022` branch and create a new pull request. 


## Big Contributions

We welcome large contributions to <span style="font-variant:small-caps;">Combine</span>. 
Note, however, that we also follow long term planning, and there is a dedicated group stewarding the overall direction and development of the code.

This means that the code development should fit in with our long term vision;
if you have an idea for a big improvement or change it will be most efficient if you [contact us](mailto:cms-cat-stats-conveners@cern.ch) first, in order to ensure that we can integrate it as seamlessly as possible into our plans.
This will simplify any potential conflicts when you make your pull request.

## Requested Contributions

As part of the long term planning, we have a number of changes we are targeting, but have not yet had a chance to implement.
As well as the [issues](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/issues) listed on the github issues tracker, you can see the [projects](https://github.com/orgs/cms-analysis/projects) listed under the general cms-analysis organization, where we have defined several projects and areas we are targeting.
If you're interested in getting involved in any of these projects please contact us at [cms-cat-stats-conveners@cern.ch](mailto:cms-cat-stats-conveners@cern.ch).


