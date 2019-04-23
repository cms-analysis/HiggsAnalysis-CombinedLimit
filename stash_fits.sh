#!/bin/bash

DIR=$1
mkdir $DIR

year=$2

mv  ws_${year}_RPV_550.root $DIR
mv  D1_CMS_th1x_prefit.png $DIR
mv  D2_CMS_th1x_prefit.png $DIR
mv  D3_CMS_th1x_prefit.png $DIR
mv  D4_CMS_th1x_prefit.png $DIR
mv  D1_CMS_th1x_fit_b.png $DIR
mv  D2_CMS_th1x_fit_b.png $DIR
mv  D3_CMS_th1x_fit_b.png $DIR
mv  D4_CMS_th1x_fit_b.png $DIR
mv  covariance_fit_b.png $DIR
mv  D1_CMS_th1x_fit_s.png $DIR
mv  D2_CMS_th1x_fit_s.png $DIR
mv  D3_CMS_th1x_fit_s.png $DIR
mv  D4_CMS_th1x_fit_s.png $DIR
mv  covariance_fit_s.png $DIR
mv  combine_logger.out $DIR
mv  fitDiagnostics${year}RPV550.root $DIR
mv  higgsCombine${year}RPV550.FitDiagnostics.mH550.MODELRPV.root $DIR
mv  log_${year}RPV550_FitDiag.txt $DIR