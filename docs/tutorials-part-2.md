

### Datacard for Shape analyses

The [datacard](#DataCard) has to be supplemented with two extensions: 1 a new block of lines defining how channels and processes are mapped into shapes 1 the block for systematics that can contain also rows with shape uncertainties.

The expected shape can be parametric or not parametric. In the first case the parametric pdfs have to be given as input to the tool. In the latter case, for each channel, histograms have to be provided for the expected shape of each process. For what concerns data, they have to be provided as input to the tool as a histogram to perform a *binned* shape analysis and as a tree to perform an *unbinned* shape analysis.

##### Not parametric shapes and uncertainties

For each channel, histograms have to be provided for the observed shape and for the expected shape of each process.

-   Within each channel, all histograms must have the same binning.
-   The normalization of the data histogram must correspond to the number of observed events
-   The normalization of the expected histograms must match the expected yields

The combine tool can take as input histograms saved as TH1 or as RooAbsHist in a RooFit workspace (an example of how to create a RooFit workspace and save histograms is available in [github](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/data/benchmarks/shapes/make_simple_shapes.cxx)).

Shape uncertainties can be taken into account by vertical interpolation of the histograms, like in the HIG-10-002 analysis. The shapes are interpolated quadratically for shifts below 1σ and linearly beyond. The normalizations are interpolated linearly in log scale just like we do for log-normal uncertainties.

For each shape uncertainty and process/channel affected by it, two additional input shapes have to be provided, obtained shifting that parameter up and down by one standard deviation. When building the likelihood, each shape uncertainty is associated to a nuisance parameter taken from a unit gaussian distribution, which is used to interpolate or extrapolate using the specified histograms.

For each given source of shape uncertainty, in the part of the datacard containing shape uncertainties (last block), there must be a row

-   \_*name* %RED%shape<span class="twiki-macro ENDCOLOR"></span> *effect for each process and channel\_*

The effect can be "-" or 0 for no effect, 1 for normal effect, and possibly something different from 1 to test larger or smaller effects (in that case, the unit gaussian is scaled by that factor before using it as parameter for the interpolation)
fre
The block of lines defining the mapping (first block in the datacard) contains one or more rows in the form

-   **<span class="twiki-macro RED"></span> shapes <span class="twiki-macro ENDCOLOR"></span> *process* *channel* *file* *histogram* \[ \_histogram\_with*systematics* \]**

In this line

-   ***process*** is any one the process names, or **\*** for all processes, or **data\_obs** for the observed data
-   ***channel*** is any one the process names, or **\*** for all channels
-   \_*file*, *histogram* and \_histogram\_with*systematics\_* identify the names of the files and of the histograms within the file, after doing some replacements (if any are found):
    -   **$PROCESS** is replaced with the process name (or "**data\_obs**" for the observed data)
    -   **$CHANNEL** is replaced with the channel name
    -   **$SYSTEMATIC** is replaced with the name of the systematic + (**Up, Down**)
    -   **$MASS** is replaced with the higgs mass value which is passed as option in the command line used to run the limit tool

The datacard in [simple-shapes-TH1.txt](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/data/benchmarks/shapes/simple-shapes-TH1.txt) is a clear example of how to include shapes in the datacard. In the first block the following line specifies the shape mapping:

    shapes * * simple-shapes-TH1.root $PROCESS $PROCESS_$SYSTEMATIC

The last block concerns the treatment of the systematics affecting shapes. In this part the two uncertainties effecting on the shape are listed.

    alpha  shape    -           1   uncertainty on background shape and normalization
    sigma  shape    0.5         -   uncertainty on signal resolution. Assume the histogram is a 2 sigma shift, 
    #                                so divide the unit gaussian by 2 before doing the interpolation

There are two options for the interpolation algorithm in the "shape" uncertainty. Putting **`shape`** will result in a quadratic interpolation (within +/-1 sigma) and a linear extrapolation (beyond +/-1 sigma) of the **fraction of events in each bin** - i.e the histograms are first normalised before interpolation. Putting **`shapeN`** while instead base the interpolation on the logs of the fraction in each bin. The total normalisation is interpolated using an asymmetric log-normal so that the effect of the systematic on both the shape and normalisation are accounted for. The following image shows a comparison of those two algorithms for this datacard.

<span class="twiki-macro TWISTY" mode="div" showlink="Show comparison " hidelink="Hide " firststart="hide" showimgright="%ICONURLPATH{toggleopen-small}%" hideimgright="%ICONURLPATH{toggleclose-small}%"></span>

&lt;img alt="compare\_shape\_algo.png" src="%ATTACHURLPATH%/compare\_shape\_algo.png" /&gt;

<span class="twiki-macro ENDTWISTY"></span>

In this case there are two processes, *signal* and *background*, and two uncertainties affecting background (*alpha*) and signal shape (*sigma*). Within the root file 2 histograms per systematic have to be provided, they are the shape obtained, for the specific process, shifting up and down the parameter associated to the uncertainty: \_background*alphaUp* and \_background\_alphaDown, \_signal*sigmaUp* and \_signal*sigmaDown*. This is the content of the root file [simple-shapes-TH1.root ](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/data/benchmarks/shapes/simple-shapes-TH1.root) associated to the datacard [simple-shapes-TH1.txt](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/data/benchmarks/shapes/simple-shapes-TH1.txt):

    root [0] 
    Attaching file simple-shapes-TH1.root as _file0...
    root [1] _file0->ls()
    TFile**     simple-shapes-TH1.root  
     TFile*     simple-shapes-TH1.root  
      KEY: TH1F signal;1    Histogram of signal__x
      KEY: TH1F signal_sigmaUp;1    Histogram of signal__x
      KEY: TH1F signal_sigmaDown;1  Histogram of signal__x
      KEY: TH1F background;1    Histogram of background__x
      KEY: TH1F background_alphaUp;1    Histogram of background__x
      KEY: TH1F background_alphaDown;1  Histogram of background__x
      KEY: TH1F data_obs;1  Histogram of data_obs__x
      KEY: TH1F data_sig;1  Histogram of data_sig__x

<span class="twiki-macro TWISTY" mode="div" showlink="Show an example " hidelink="Hide " firststart="hide" showimgright="%ICONURLPATH{toggleopen-small}%" hideimgright="%ICONURLPATH{toggleclose-small}%"></span>

For example, without shape uncertainties you could have just one row with
`shapes * * shapes.root $CHANNEL/$PROCESS`
Then for a simple example for two channels "e", "mu" with three processes "higgs", "zz", "top" you should create a rootfile that contains the following

| histogram     | meaning                                      |
|:--------------|:---------------------------------------------|
| `e/data_obs`  | observed data in electron channel            |
| `e/higgs`     | expected shape for higgs in electron channel |
| `e/zz`        | expected shape for ZZ in electron channel    |
| `e/top`       | expected shape for top in electron channel   |
| `mu/data_obs` | observed data in muon channel                |
| `mu/higgs`    | expected shape for higgs in muon channel     |
| `mu/zz`       | expected shape for ZZ in muon channel        |
| `mu/top`      | expected shape for top in muon channel       |

If you also have one uncertainty that affects the shape, e.g. jet energy scale, you should create shape histograms for the jet energy scale shifted up by one sigma, you could for example do one folder for each process and write a like like
`shapes * * shapes.root $CHANNEL/$PROCESS/nominal  $CHANNEL/$PROCESS/$SYSTEMATIC`
or just attach a postifx to the name of the histogram
`shapes * * shapes.root $CHANNEL/$PROCESS  $CHANNEL/$PROCESS_$SYSTEMATIC`

<span class="twiki-macro ENDTWISTY"></span>

**Parametric shapes and uncertainties**

If you want to combine a channel with shapes with one that is a simple counting experiment, you have to declare some fake shapes in the datacard of the counting experiment. This can be done simply by adding a line to the datacard

-   **<span class="twiki-macro RED"></span> shapes <span class="twiki-macro ENDCOLOR"></span> *process* *channel* %RED%FAKE%ENDCOLOR%**

<span class="twiki-macro TWISTY" mode="div" showlink="Implementation details and status " hidelink="Hide " firststart="hide" showimgright="%ICONURLPATH{toggleopen-small}%" hideimgright="%ICONURLPATH{toggleclose-small}%"></span>

As of tag **`T01-06-00`** <span class="twiki-macro COMPLETE5"></span> (note the initial **`T`** instead of **`V`**, since it's still only for testing)

-   **Shapes:** <span class="twiki-macro ICON">choice-yes</span>
    -   <span class="twiki-macro X"></span> If doing toy mc generation, must run with **`--generateBinnedWorkaround`** (or **`-U`** / **`--unbinned`**), due to a bug in RooFit
-   **Shape uncertainties:** <span class="twiki-macro ICON">choice-yes</span>
    -   shapes with vertical shape morphing (quadratic up to 1 sigma, then linear)
    -   exponential morphing on the normalization (just like for an asymmetric log-normal uncertainty)
    -   one can select linear morphing (`shapeL`) or multiplicative morphing (`shapeN`) as well, but with some caveats:
        -   only one morphing algorithm can be used for a given shape
        -   multiplicative morphing is applied to the *normalized* shapes, and then the normalization is interpolated separately. also, note that this is **much slower** than quadratic or linear morphing (factor 10-100)
    -   if you get results that look unreasonable, it could be that the truncation effects due to some part of the shape becoming negative are sizable. you can partially recover from this problem forcing roofit to re-normalize explicitly the function (just replace **`shape`** with ==shape\*==; note, will be much slower), but it's likely that you need some other approach to handle the shape systematics

This includes also a few other features not in this twiki:

-   Support for RooFit shapes (specifying &lt;tt&gt; *workspaceName* **:***objectName*&lt;/tt&gt; in place of *histogramName*):
    -   **Binned datasets** (RooDataHist) for data and for mc processes: <span class="twiki-macro ICON">choice-yes</span>
    -   **Unbinned datasets** (RooDataSet) for data: <span class="twiki-macro ICON">choice-yes</span>
    -   **Arbitrary shapes** (RooAbsPdf) for mc processes: <span class="twiki-macro ICON">choice-yes</span>
        -   <span class="twiki-macro ICON">choice-yes</span> Parametric shape uncertainties can be specified among the other systematics in a line
            &lt;tt&gt;*name* **param** *mean* *uncertainty* *\[range\]* &lt;/tt&gt;.
            The uncertainty can be either the sigma of a gaussian or `-xx/+yy` (with no spaces) for a bifurcated gaussian.
            The range is optional. if present, it should be `[hi,lo]` (with no spaces).
        -   <span class="twiki-macro X"></span> Normalization will always be taken from text datacard
    -   Unbinned templates (RooDataSet) for mc processes : <span class="twiki-macro ICON">choice-yes</span> (but not tested)
    -   <span class="twiki-macro X"></span> arbitrary shapes must have unique parameter names and they should match with the names of the systematics
-   Support for plain **`TTrees`** as inputs instead of RooDataSet: <span class="twiki-macro ICON">choice-yes</span> (but not tested)
-   Mixing of histograms with RooFit shapes is not supported.

<span class="twiki-macro ENDTWISTY"></span>

### Binned shape analysis

[See the 2014 Data Analysis School tutorial.](%SCRIPTURL{"view"}%auth/CMS/SWGuideCMSDataAnalysisSchool2014HiggsCombPropertiesExercise#A_shape_analysis_using_templates)

### Unbinned shape analysis

<span class="twiki-macro TWISTY" mode="div" showlink="Example from Higgs to gamma gamma" hidelink="Hide " firststart="hide" showimgright="%ICONURLPATH{toggleopen-small}%" hideimgright="%ICONURLPATH{toggleclose-small}%"></span>

*This example is taken from [2014 Data Analysis School](%SCRIPTURL{"view"}%auth/CMS/SWGuideCMSDataAnalysisSchool2014HiggsCombPropertiesExercise#A_parametric_shape_analysis_H) and use H-&gt;gg datacards as an example. This is an example of an unbinned, parametric analysis.*

In some cases, it can be convenient to describe the expected signal and background shapes in terms of analytical functions rather than templates; a typical example are the searches where the signal is apparent as a narrow peak over a smooth continuum background. In this context, uncertainties affecting the shapes of the signal and backgrounds can be implemented naturally as uncertainties on the parameters of those analytical functions. It is also possible to adapt an agnostic approach in which the parameters of the background model are left freely floating in the fit to the data, i.e. only requiring the background to be well described by a smooth function.

Technically, this is implemented by means of the RooFit package, that allows writing generic probability density functions, and saving them into ROOT files. The pdfs can be either taken from RooFit's standard library of functions (e.g. Gaussians, polynomials, ...) or hand-coded in C++, and combined together to form even more complex shapes.

A prototypical case for this kind of analysis is H→γγ analysis. For this excercise, we will use a datacard that contains only the 8 TeV data for four event categories: **cat0** and **cat1** are categories of untagged events containing good quality diphotons, the former of the two with a higer purity but lower event yield obtained preferentially selecting high p&lt;sub&gt;T&lt;/sub&gt; diphotons; **cat4** and **cat5** are categories of di-jet events with different levels of tightness (and so also different level of signal contamination from gluon fusion).

The datacard is the following:

    imax 4 number of bins
    jmax 5 number of processes minus 1
    kmax * number of nuisance parameters
    ----------------------------------------------------------------------------------------------------------------------------------
    shapes WH        cat0      hgg.inputsig_8TeV_MVA.root wsig_8TeV:hggpdfrel_wh_cat0
    shapes ZH        cat0      hgg.inputsig_8TeV_MVA.root wsig_8TeV:hggpdfrel_zh_cat0
    shapes bkg_mass  cat0      hgg.inputbkgdata_8TeV_MVA.root cms_hgg_workspace:pdf_data_pol_model_8TeV_cat0
    shapes data_obs  cat0      hgg.inputbkgdata_8TeV_MVA.root cms_hgg_workspace:roohist_data_mass_cat0
    shapes ggH       cat0      hgg.inputsig_8TeV_MVA.root wsig_8TeV:hggpdfrel_ggh_cat0
    shapes qqH       cat0      hgg.inputsig_8TeV_MVA.root wsig_8TeV:hggpdfrel_vbf_cat0
    shapes ttH       cat0      hgg.inputsig_8TeV_MVA.root wsig_8TeV:hggpdfrel_tth_cat0
    [... same as above for cat1, cat4, cat5 ...]
    ----------------------------------------------------------------------------------------------------------------------------------
    bin          cat0         cat1         cat4         cat5       
    observation  -1.0         -1.0         -1.0         -1.0       
    ----------------------------------------------------------------------------------------------------------------------------------
    bin                                      cat0         cat0         cat0         cat0         cat0         cat0         cat1         cat1         cat1         cat1         cat1         cat1         cat4         cat4         cat4         cat4         cat4         cat4         cat5         cat5         cat5         cat5         cat5         cat5       
    process                                  ZH           qqH          WH           ttH          ggH          bkg_mass     ZH           qqH          WH           ttH          ggH          bkg_mass     ZH           qqH          WH           ttH          ggH          bkg_mass     ZH           qqH          WH           ttH          ggH          bkg_mass   
    process                                  -4           -3           -2           -1           0            1            -4           -3           -2           -1           0            1            -4           -3           -2           -1           0            1            -4           -3           -2           -1           0            1          
    rate                                     6867.0000    19620.0000   12753.0000   19620.0000   19620.0000   1.0000       7259.4000    19620.0000   12360.6000   19620.0000   19620.0000   1.0000       7063.2000    19620.0000   12556.8000   19620.0000   19620.0000   1.0000       3924.0000    19620.0000   15696.0000   19620.0000   19620.0000   1.0000     
    ----------------------------------------------------------------------------------------------------------------------------------
    CMS_eff_j               lnN              0.999125     0.964688     0.999125     0.998262     0.996483     -            0.999616     0.980982     0.999616     0.99934      0.999012     -            1.02         1.02         1.02         1.02         1.02         -            1.02         1.02         1.02         1.02         1.02         -          
    CMS_hgg_JECmigration    lnN              -            -            -            -            -            -            -            -            -            -            -            -            0.853986     0.995971     0.853986     0.846677     0.927283     -            1.025        1.005        1.025        1.025        1.025        -          
    CMS_hgg_UEPSmigration   lnN              -            -            -            -            -            -            -            -            -            -            -            -            0.737174     0.991941     0.737174     0.724019     0.86911      -            1.045        1.01         1.045        1.045        1.045        -          
    CMS_hgg_eff_MET         lnN              -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -          
    CMS_hgg_eff_e           lnN              -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -          
    CMS_hgg_eff_m           lnN              -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -            -          
    CMS_hgg_eff_trig        lnN              1.01         1.01         1.01         1.01         1.01         -            1.01         1.01         1.01         1.01         1.01         -            1.01         1.01         1.01         1.01         1.01         -            1.01         1.01         1.01         1.01         1.01         -          
    CMS_hgg_n_id            lnN              1.034/0.958  1.039/0.949  1.034/0.958  1.053/0.915  1.035/0.958  -            1.042/0.936  1.038/0.948  1.042/0.936  1.053/0.909  1.038/0.954  -            1.016/0.972  1.022/0.963  1.016/0.972  1.017/0.967  1.019/0.964  -            1.024/0.959  1.030/0.948  1.024/0.959  1.014/0.975  1.023/0.962  -          
    CMS_hgg_n_pdf_1         lnN              -            0.998/0.996  -            -            1.002/0.998  -            -            0.992/0.999  -            -            1.001/1.000  -            -            0.999/0.998  -            -            1.003/0.999  -            -            0.996/1.000  -            -            1.002/0.999  -          
    CMS_hgg_n_pdf_10        lnN              -            1.028/1.004  -            -            0.998/0.998  -            -            1.051/0.997  -            -            1.004/0.999  -            -            1.020/1.000  -            -            1.000/0.998  -            -            1.010/0.999  -            -            0.999/1.000  -          
    [... and more rows like this ...]
    CMS_hgg_n_sc_gf         lnN              -            -            -            -            0.842/1.123  -            -            -            -            -            0.976/1.031  -            -            -            -            -            0.858/1.095  -            -            -            -            -            0.880/1.083  -          
    CMS_hgg_n_sc_vbf        lnN              -            0.993/1.010  -            -            -            -            -            0.994/0.994  -            -            -            -            -            0.997/1.001  -            -            -            -            -            0.999/0.996  -            -            -            -          
    CMS_hgg_n_sigmae        lnN              0.950/1.099  0.943/1.114  0.950/1.099  0.956/1.089  0.944/1.112  -            0.954/1.093  0.948/1.104  0.954/1.093  0.971/1.057  0.919/1.165  -            0.996/1.006  0.992/1.016  0.996/1.006  0.995/1.010  0.994/1.012  -            0.991/1.019  0.985/1.031  0.991/1.019  0.995/1.010  0.988/1.026  -          
    CMS_id_eff_eb           lnN              1.01999      1.020047     1.01999      1.02002      1.020022     -            1.018535     1.019332     1.018535     1.018642     1.019654     -            1.017762     1.017952     1.017762     1.018824     1.01812      -            1.016723     1.016872     1.016723     1.01884      1.017172     -          
    CMS_id_eff_ee           lnN              1.000284     1.000137     1.000284     1.000206     1.0002       -            1.004035     1.001979     1.004035     1.003759     1.001148     -            1.006058     1.005556     1.006058     1.003308     1.005127     -            1.008752     1.008351     1.008752     1.003254     1.007586     -          
    JEC                     lnN              0.995187     0.938204     0.995187     0.990442     0.980656     -            0.997891     0.966719     0.997891     0.996368     0.994567     -            1.11         1.035        1.11         1.11         1.11         -            1.11         1.035        1.11         1.11         1.11         -          
    QCDscale_VH             lnN              0.982/1.021  -            0.982/1.021  -            -            -            0.982/1.021  -            0.982/1.021  -            -            -            0.982/1.021  -            0.982/1.021  -            -            -            0.982/1.021  -            0.982/1.021  -            -            -          
    QCDscale_ggH            lnN              -            -            -            -            0.918/1.076  -            -            -            -            -            0.918/1.076  -            -            -            -            -            0.918/1.076  -            -            -            -            -            0.918/1.076  -          
    QCDscale_qqH            lnN              -            0.992/1.003  -            -            -            -            -            0.992/1.003  -            -            -            -            -            0.992/1.003  -            -            -            -            -            0.992/1.003  -            -            -            -          
    QCDscale_ttH            lnN              -            -            -            0.906/1.041  -            -            -            -            -            0.906/1.041  -            -            -            -            -            0.906/1.041  -            -            -            -            -            0.906/1.041  -            -          
    UEPS                    lnN              0.988624     0.858751     0.988624     0.977408     0.954278     -            0.995014     0.923929     0.995014     0.991416     0.987158     -            1.26         1.08         1.26         1.26         1.26         -            1.26         1.08         1.26         1.26         1.26         -          
    lumi_8TeV               lnN              1.044        1.044        1.044        1.044        1.044        -            1.044        1.044        1.044        1.044        1.044        -            1.044        1.044        1.044        1.044        1.044        -            1.044        1.044        1.044        1.044        1.044        -          
    pdf_gg                  lnN              -            -            -            0.920/1.080  0.930/1.076  -            -            -            -            0.920/1.080  0.930/1.076  -            -            -            -            0.920/1.080  0.930/1.076  -            -            -            -            0.920/1.080  0.930/1.076  -          
    pdf_qqbar               lnN              0.958/1.042  0.972/1.026  0.958/1.042  -            -            -            0.958/1.042  0.972/1.026  0.958/1.042  -            -            -            0.958/1.042  0.972/1.026  0.958/1.042  -            -            -            0.958/1.042  0.972/1.026  0.958/1.042  -            -            -          
    vtxEff                  lnN              0.991/1.025  0.989/1.030  0.991/1.025  0.993/1.020  0.989/1.030  -            0.990/1.026  0.990/1.027  0.990/1.026  0.994/1.016  0.984/1.042  -            1.000/0.999  0.997/1.007  1.000/0.999  1.000/1.003  0.998/1.006  -            0.999/1.004  0.998/1.006  0.999/1.004  1.000/1.000  0.997/1.007  -          
    CMS_hgg_nuissancedeltamcat4  param  0.0 0.001458
    CMS_hgg_nuissancedeltafracright_8TeV  param  1.0 0.002000
    CMS_hgg_nuissancedeltamcat1  param  0.0 0.001470
    CMS_hgg_nuissancedeltamcat0  param  0.0 0.001530
    CMS_hgg_nuissancedeltasmearcat4  param  0.0 0.001122
    CMS_hgg_nuissancedeltasmearcat1  param  0.0 0.001167
    CMS_hgg_nuissancedeltasmearcat0  param  0.0 0.001230
    CMS_hgg_globalscale  param  0.0 0.004717

The first difference compared to the template datacard is in the **`shapes`** line; let's take for example these two lines

    shapes ggH       cat0      hgg.inputsig_8TeV_MVA.root wsig_8TeV:hggpdfrel_ggh_cat0
    shapes data_obs  cat0      hgg.inputbkgdata_8TeV_MVA.root cms_hgg_workspace:roohist_data_mass_cat0

In the datacard using templates, the column after the file name would have been the name of the histogram. Here, instead, we found two names, separated by a colon (**`:`**): the first part identifies the name of the [RooWorkspace](http://root.cern.ch/root/htmldoc/RooWorkspace.html) containing the pdf, and the second part the name of the [RooAbsPdf](http://root.cern.ch/root/htmldoc/RooAbsPdf.html) inside it (or, for the observed data, the [RooAbsData](http://root.cern.ch/root/htmldoc/RooAbsData.html)).

Let's inspect this workspace, starting with the data

    %CODE{"cpp"}%
    TFile *fDat = TFile::Open("hgg.inputbkgdata_8TeV_MVA.root");
    RooAbsData *data = cms_hgg_workspace-&gt;data("roohist_data_mass_cat0");
    data-&gt;Print("")
    // --&gt; RooDataHist::roohist_data_mass_cat0[CMS_hgg_mass] = 160 bins (1449 weights)
    // so, we have a binned dataset, whose variable is called CMS_hgg_mass:
    RooRealVar *mass = cms_hgg_workspace-&gt;var("CMS_hgg_mass");
    mass-&gt;Print("");
    // RooRealVar::CMS_hgg_mass = 140  L(100 - 180)

    // we can make a plot of the dataset with the following
    RooPlot *plot = mass-&gt;frame();
    data-&gt;plotOn(plot);
    plot-&gt;Draw();
    %ENDCODE%

&lt;img alt="datacards\_Hgg\_data\_cat0.png" src="%ATTACHURLPATH%/datacards\_Hgg\_data\_cat0.png" /&gt;

Now let's look also at the signal

    %CODE{"cpp"}%
    TFile *sig = TFile::Open("hgg.inputsig_8TeV_MVA.root");
    RooWorkspace *wsig_8TeV = (RooWorkspace*)sig-&gt;Get("wsig_8TeV")  // necessary to Get the workspace in ROOT6 to avoid reloading
    RooAbsPdf *ggH = wsig_8TeV-&gt;pdf("hggpdfrel_ggh_cat0");
    ggH-&gt;Print("");
    // --&gt; RooAddPdf::hggpdfrel_ggh_cat0[ hist_func_frac_g0_ggh_cat0 * hgg_gaus_g0_ggh_cat0 + hist_func_frac_g1_ggh_cat0 * hgg_gaus_g1_ggh_cat0 + const_func_frac_g2_ggh_cat0 * hgg_gaus_g2_ggh_cat0 ] = 0.604097
    // this appears to be a linear combination of multiple gaussian pdfs (hgg_gaus_g0_ggh_cat0, hgg_gaus_g1_ggh_cat0, ...) with coefficents hist_func_frac_g0_ggh_cat0, hist_func_frac_g1_ggh_cat0

    // let's get the list of the parameters that describe it
    RooArgSet *params = ggH-&gt;getParameters(*data);
    params-&gt;Print("");
    // --&gt; RooArgSet::parameters = (CMS_hgg_globalscale,CMS_hgg_nuissancedeltamcat0,CMS_hgg_nuissancedeltasmearcat0,MH)

    // MH is a special parameter, which combine and text2workspace set to the Higgs mass hypothesis.
    // now since we're just looking into the input workspace, we can set it by hand
    wsig_8TeV-&gt;var("MH")-&gt;setVal(125.7);

    // Now we can make a plot of the pdf, and show also the contributions from the different gaussians from which it is composed
    RooPlot *plot =  wsig_8TeV-&gt;var("CMS_hgg_mass")-&gt;frame();
    ggH-&gt;plotOn(plot);
    ggH-&gt;plotOn(plot, RooFit::Components("hgg_gaus_g0_ggh_cat0"), RooFit::LineColor(kRed));
    ggH-&gt;plotOn(plot, RooFit::Components("hgg_gaus_g1_ggh_cat0"), RooFit::LineColor(209));
    ggH-&gt;plotOn(plot, RooFit::Components("hgg_gaus_g2_ggh_cat0"), RooFit::LineColor(222));
    plot-&gt;Draw();
    %ENDCODE%

&lt;img alt="datacards\_Hgg\_ggH\_cat0.png" src="%ATTACHURLPATH%/datacards\_Hgg\_ggH\_cat0.png" /&gt;

#### Parametric signal normalization

There is also another feature of the H→γγ datacard that should catch your eye quicky: the event yields are remarkably strange: the the signal yield is **`19620.0000`** events for ggH, qqH and ttH, the background yield is **`1.0`**, and observed event yield is **`-1.0`**.

Let's start with the simple case: an event yield of **`-1`** just instructs text2workspace and combine to take the yield from the corresponding dataset in the input rootfile, avoiding the need of writing it in the text datacard, but also making the datacard less human-readable. Incidentally, this feature can be used also for datacards that use root histograms like the H→ττ example above.

Now, the signal and background yields: 19k signal events from each production mode with just one background event would be really nice to have, but of course it can't be true. The way the H→γγ works is by relying on an additional feature of text2workspace: in addition to providing generic shapes for the signals and backgrounds, it is also possible to provide generic functions that describe the expected signal and background yields. This additional functions are multiplied by the number in the **`rate`** column of the datacard to obtain the final expected event yield. This feature allows e.g. to use the very same datacard to describe all possible different Higgs boson mass hypotheses by just parameterizing properly the expected yield as function of the **`MH`** variable.

At present, this feature is not indicated by any special line the datacard: simply, whenever a RooAbsPdf is loaded, text2workspace will search also for a RooAbsReal object with the same name but a **`_norm`** postfix, and if present it will use it to scale the event yield. <span class="twiki-macro RED"></span> WARNING <span class="twiki-macro ENDCOLOR"></span> As with **all** parameters in the workspace, if this RooAbsReal object is not set constant, (i.e is a RooRealVar which you have not set constant or a function of non-constant vars) it will be floating in the fit. This is especially problematic for the signal as the default parameter added to models for limits/p-values `r` (the floating signal strength) will be degenerate with this **`_norm`** parameter. In this H→γ&gamma example, the object has been setup to be constant to avoid this.

Note: The newest version of combine will not accept RooExtendedPdfs as an input anymore. This is to alleviate a bug that lead to improper treatment of normalization when using multiple RooExtendedPdfs to describe a single process. Simply follow the above instructions and combine will create the appropriate extended pdf as long as the name of the RooAbsReal is pdfname\_norm. Make sure the variable is not set to constant if you intend to have the normalization float.

Armed with this new piece of knowledge, we can now determine what is the expected signal yield for ggH in category 0:

    %CODE{"cpp"}%
    TFile *sig = TFile::Open("hgg.inputsig_8TeV_MVA.root");
    RooWorkspace *wsig_8TeV = (RooWorkspace*)sig-&gt;Get("wsig_8TeV")  // necessary to Get the workspace in ROOT6 to avoid reloading
    wsig_8TeV-&gt;var("MH")-&gt;setVal(125.7);
    RooAbsPdf *ggH = wsig_8TeV-&gt;pdf("hggpdfrel_ggh_cat0");
    RooAbsReal *ggH_norm = wsig_8TeV-&gt;function("hggpdfrel_ggh_cat0_norm");
    cout &lt;&lt; ggH_norm-&gt;getVal()*19620.0000 &lt;&lt; endl;
    // --&gt; 12.3902
    %ENDCODE%

This approach can also be used to make a background freely floating, by associating to them as normalization term a floating RooRealVar. However, this can be achieved in a more transparent way by putting instead a normalization uncertainty on that background using a flat pdf: e.g. to leave the background floating between 50% and 200% of its input prediction, a **`lnU`** systematic can be used with κ`2 ( ==lnU=` has a syntax like **`lnN`**, but produces a uniform pdf between 1/κ and κ rather than a log-normal; see [SWGuideHiggsAnalysisCombinedLimit](%SCRIPTURL{"view"}%auth/CMS/SWGuideHiggsAnalysisCombinedLimit#How_to_prepare_the_datacard))

#### Shape uncertainties using parameters

The part of the H→γγ datacard related to the systematics starts with many lines of log-normals that should already be familiar, except possibly for the notation with two numbers separated by a slash (e.g. 0.950/1.099). This notation is used for asymmetrical uncertainties: a log-normal with 0.950/1.099 means that at -1σ the yield is scaled down by a factor 0.95, while at +1σ the yield is scaled up by a factor 1.099.

The last part of the datacard contains some lines that use a different syntax, e.g.

    CMS_hgg_globalscale  param  0.0 0.004717

These lines encode uncertainties on the parameters of the signal and background pdfs. The example line quoted here informs text2workspace that the parameter **`CMS_hgg_globalscale`** is to be assigned a Gaussian uncertainty of ±0.004717 around its mean value of zero (0.0). One can change the mean value from 0 to 1 (or really any value, if one so chooses) if the parameter in question is multiplicative instead of additive.

The effect can be visualized from RooFit, e.g.

    %CODE{"cpp"}%
    TFile *sig = TFile::Open("hgg.inputsig_8TeV_MVA.root");
    RooWorkspace *wsig_8TeV = (RooWorkspace*)sig-&gt;Get("wsig_8TeV")  // necessary to Get the workspace in ROOT6 to avoid reloading

    wsig_8TeV-&gt;var("MH")-&gt;setVal(125.7);
    RooAbsPdf *ggH = wsig_8TeV-&gt;pdf("hggpdfrel_ggh_cat0");

    // prepare the canvas
    RooPlot *plot =  wsig_8TeV-&gt;var("CMS_hgg_mass")-&gt;frame();

    // plot nominal pdf
    ggH-&gt;plotOn(plot, RooFit::LineColor(kBlack));

    // plot minus 3 sigma pdf
    wsig_8TeV-&gt;var("CMS_hgg_globalscale")-&gt;setVal(-3*0.004717);
    ggH-&gt;plotOn(plot, RooFit::LineColor(kBlue));

    // plot plus 3 sigma pdf
    wsig_8TeV-&gt;var("CMS_hgg_globalscale")-&gt;setVal(+3*0.004717);
    ggH-&gt;plotOn(plot, RooFit::LineColor(kRed));
    plot-&gt;Draw();
    %ENDCODE%

&lt;img alt="datacards\_Hgg\_ggH\_cat0\_syst.png" src="%ATTACHURLPATH%/datacards\_Hgg\_ggH\_cat0\_syst.png" /&gt;

Note that if one wants to specify a parameter that is freely floating across its given range, and not gaussian constrained, the following syntax is used:

     CMS_my_bg_param1  flatParam 

The Hgg and HZg analyses use this syntax when adding in the parameters that correspond to their background shapes.

\#ParametricModelAndBinnedData

#### Caveat on using parametric pdfs with binned datasets

Users should be aware of a feature that affects the use of parametric pdfs together with binned datasets.

RooFit uses the integral of the pdf, computed analytically (or numerically, but disregarding the binning), to normalize it, but then computes the expected event yield in each bin evaluating only the pdf at the bin center. This means that if the variation of the pdf is sizeable within the bin then there is a mismatch between the sum of the event yields per bin and the pdf normalization, and that can cause a bias in the fits (more properly, the bias is there if the contribution of the second derivative integrated on the bin size is not negligible, since for linear functions evaluating them at the bin center is correct).

So, it is recommended to use bins that are significantly finer than the characteristic scale of the pdfs - which would anyway be the recommended thing even in the absence of this feature.

Obviously, this caveat does not apply to analyses using templates (they're constant across each bin, so there's no bias), or using unbinned datasets.

<span class="twiki-macro ENDTWISTY"></span>

### Analysis with more generic models



#### Multidimensional fits



#### Feldman-Cousins regions

The F-C procedure for a generic model is:

-   use as test statistics the profile likelihood *q(x) = - 2 ln L(data|x)/L(data|x-hat)* where *x* is a point in the parameter space, and *x-hat* are the point corresponding to the best fit (nuisance parameters are profiled both at numerator and at denominator)
-   for each point *x*:
    -   compute the observed test statistics *q&lt;sub&gt;obs&lt;/sub&gt;(x)*
    -   compute the expected distribution of *q(x)* under the hypothesis of x.
    -   accept the point in the region if *P(q(x) &lt; q&lt;sub&gt;obs&lt;/sub&gt;(x) | x) &lt; CL*

In combine, you can perform this test on each individual point (param1,param2,...) = (value1,value2,...) by doing

    combine workspace.root -M HybridNew --freq --testStat=PL --rule=CLsplusb --singlePoint  param1=value1,param2=value2,param3=value3,...   [other options of HybridNew]

The point belongs to your confidence region if CL&lt;sub&gt;s+b&lt;/sub&gt; is larger than 1-CL (e.g. 0.3173 for a 1-sigma region, CL=0.6827).

Imposing physical boundaries (such as requiring mu&gt;0) can be achieved by setting the ranges of the physics model parameters using

     --setPhysicsModelParameterRanges param1=param1_min,param1_max:param2=param2_min,param2_max ....

. If there is no upper/lower boundary, just set that value to something far from the region of interest.

As in general for HybridNew, you can split the task into multiple tasks and then merge the outputs, as described in the [HybridNew chapter](SWGuideHiggsAnalysisCombinedLimit#HybridNew_algorithm).

For uni-dimensional models only, and if the parameter behaves like a cross-section, the code is somewhat able to do interpolation and determine the values of your parameter on the contour (just like it does for the limits). In that case, the syntax is the same as per the CLs limits with [HybridNew chapter](HybridNew chapter) except that you want **`--testStat=PL --rule=CLsplusb`** .

**Extracting Contours**

There is a tool for extracting confidence intervals and 2D contours from the output of HybridNew located in **`test/makeFCcontour.py`** providing the option **`--saveToys`** was included when running HybridNew. I can be run taking as input, the toys files (or several of them) as,

    ./makeFCcontour.py  toysfile1.root toysfile2.root .... [options] -out outputfile.root

The tool has two modes (1D and 2D). For the 1D, add the option **`--d1`** and the name of the parameter of interest **`--xvar poi_name`**. For each confidence interval desired, add any confidence level of interest using **`--cl 0.68,0.95...`** The intervals corresponding to each confidence level will be printed to the terminal. The output file will contain a graph of the parameter of interest (x) vs 1-CL&lt;sub&gt;s+b&lt;/sub&gt; used to compute the intervals.

To extract 2D contours, the names of each parameter must be given **`--xvar poi_x --yvar poi_y`**. The output will be a root file containing a 2D histogram of the confidence level (1-CL&lt;sub&gt;s+b&lt;/sub&gt;) for each point which can be used to draw 2D contours. There will also be a histogram containing the number of toys found for each point.

There are several options for reducing the running time (such as setting limits on the region of interest or the minimum number of toys required for a point to be included) Finally, adding the option **`--storeToys`** will add histograms in for each point to the output file of the test-statistic distribution. This will increase the momory usage however as all of the toys will be stored.

### Signal Hypothesis separation

In some cases, instead of separating a signal from a background, you might want to separate a signal of one type from a signal of another type (e.g. scalar vs pseudo-scalar Higgs boson).

This is documented at [SWGuideHiggsCombinationSignalSeparation](HiggsWG.SWGuideHiggsCombinationSignalSeparation)



#### Throwing post-fit toys

From [here](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/pull/86):

    #build workspace for mu-mh fit
    text2workspace.py ../cards/hgg_datacard_mva_comb_bernsteins.txt -m 125 -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass --PO higgsMassRange=105,155 -o testmasshggcomb.root

    #perform s+b fit and save workspace+snapshot
    combine testmasshggcomb.root -m 125 -M MultiDimFit --saveWorkspace --verbose 9 -n mumhfit

    #throw post-fit toy with b from s+b(floating mu,mh) fit, s with r=1.0, m=best fit MH, using nuisance values and constraints re-centered on s+b(floating mu-mh) fit values (aka frequentist post-fit expected)
    #and compute post-fit expected mu uncertainty profiling MH
    combine higgsCombinemumhfit.MultiDimFit.mH125.root --snapshotName MultiDimFit -M MultiDimFit --verbose 9 -n randomtest --toysFrequentist --bypassFrequentistFit  -t -1 --expectSignal=1 -P r --floatOtherPOIs=1 --algo singles

    #throw post-fit toy with b from s+b(floating mu,mh) fit, s with r=1.0, m=128.0, using nuisance values and constraints re-centered on s+b(floating mu-mh) fit values (aka frequentist post-fit expected)
    #and compute post-fit expected significance (with MH fixed at 128 implicitly)
    combine higgsCombinemumhfit.MultiDimFit.mH125.root -m 128 --snapshotName MultiDimFit -M ProfileLikelihood --significance --verbose 9 -n randomtest --toysFrequentist --bypassFrequentistFit --overrideSnapshotMass -t -1 --expectSignal=1 --redefineSignalPOIs r --freezeNuisances MH

    #throw post-fit toy with b from s+b(floating mu,mh) fit, s with r=0.0, using nuisance values and constraints re-centered on s+b(floating mu-mh) fit values (aka frequentist post-fit expected)
    #and compute post-fit expected and observed asymptotic limit (with MH fixed at 128 implicitly)
    combine higgsCombinemumhfit.MultiDimFit.mH125.root -m 128 --snapshotName MultiDimFit -M Asymptotic --verbose 9 -n randomtest --bypassFrequentistFit --overrideSnapshotMass--redefineSignalPOIs r --freezeNuisances MH

### Modifying parameters of interest on the command line

Normally, the parameters of interest of a model are defined by the PhysicsModel used by text2workspace.

However, combine provides command line options to redefine on the fly what are the parameters of interest (useful also if the workspaces are built from external tools and not text2workspace), for setting the values of the parameters and their ranges. The MultiDimFit method already provides conventient handles to specify which parameter to analyze and how to deal with the others (floating or fixed), but these commands are more general as they apply to all statistical methods and they can affect any kind of parameters in the model, not only those that were pre-defined as parameters of interest (the wording "physics model parameter" or "nuisance" in the names of the options is there only because that's typically what they are applied to, not because their functionality is restricted to parameters of that kind)

-   **`--setPhysicsModelParameters name=value[,name2=value2,...]`** sets the starting values of the parameters, useful e.g. when generating toy MC or when also setting the parameters as fixed.
-   **`--setPhysicsModelParameterRanges name=min,max[:name2=min2,max2:...]`** sets the ranges of the parameters (useful e.g. for scanning in MultiDimFit, or for Bayesian integration)
-   **`--redefineSignalPOIs name[,name2,...]`** redefines the set of parameters of interest.
    -   if the parameters where constant in the input workspace, they are re-defined to be floating.
    -   nuisances promoted to parameters of interest are removed from the list of nuisances, and thus they are not randomized in methods that randomize nuisances (e.g. HybridNew in non-frequentist mode, or BayesianToyMC, or in toy generation with `-t` but without **`--toysFreq`**).
        This doesn't have any impact on algorithms that don't randomize nuisances (e.g. fits, Asymptotic, or HybridNew in fequentist mode) or on algorithms that treat all parameters in the same way (e.g. MarkovChainMC).
    -   Note that constraint terms for the nuisances are **dropped** after promotion to a POI using `--redefineSignalPOI`. To produce a likelihood scan for a nuisance parameter, using MultiDimFit with **`--algo grid`**, you should instead use the **`--poi`** option which will not cause the loss of the constraint term when scanning.
    -   parameters of interest of the input workspace that are not selected by this command become unconstrained nuisance parameters, but they are not added to the list of nuisances so they will not be randomized (see above)
-   **`--freezeNuisances <name>`** sets the given parameters to constant

A combination of the MultiSignalModel (defined above) and **`redefineSignalPOIs`** can be used to alter the way a datacard is interpreted to transform one signal or ordinary background process into a background with freely floating normalization, or with a fixed but different normalization (useful e.g. for cross-checks):

-   First, use the MultiSignalModel to create a workspace in which that process has its own associated signal strenght parameter (e.g. `r_B`) different from that of the nomal signal (e.g. `r`)
-   Then, fits or upper limits on `r` can be obtained running combine with **`--redefineSignalPOIs r`** so that `r_B` becomes freely floating.
-   **`setPhysicsModelParameterRanges`** can be used to redefine the range to make the process freely float only in some range (this is equivalent to adding a flat uncertainty to it, **`lnU`** or **`unif`**) while **`freezeNuisances`** and **`setPhysicsModelParameters`** can instead be used to pin it to some value.
-   fits and likelihood scans of `r_B` can be used to check whether the information on this process from the data in the signal region is consistent with the a-priori prediction for it within the respective uncertainties (and **`freezeNuisances`** can be used to shut off the uncertainties on the a-priori prediction if one wants to see only the uncertainty from the data)

### Advanced Tutorials

Follow [this link](CMS/HiggsWG/SWGuideNonStandardCombineUses) for advanced tutorials and non-standard uses of combine
