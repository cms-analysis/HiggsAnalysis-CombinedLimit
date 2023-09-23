#coded up by a.david@cern.ch in November 2015

nuisgroups = bgexp exp sigtheo bgtheo theo bg
targets = $(patsubst %,frozen_%,$(nuisgroups) all none) frozen_all_mlf frozen_none_mlf
cardname = myanalysis

all: allfits
	@echo
	@echo "Welcome to a weird way of presenting results that is halfway between slides and jupyter."
	@echo
	@echo "0) Note the text2workspace printout regarding groups, that helps figuring out what is going on:"
	@grep "group_" logs/myanalysis.T2W.log
	@echo
	@echo "1) The best-fit snaphot has the following best-fit values:"
	@tail -n28 logs/myanalysis_bestfit.log
	@echo
	@echo "=> Note how nuisances are pulled from zero."
	@echo
	@echo "2) Now, compare what happens if you freeze nuisances to the default values directly in MLF (all_mlf) or to the best-fit values using the snapshot (all):"
	@tail -n3 logs/frozen_all{_mlf,}.log | grep "==>" -A1
	@echo
	@echo "=> It makes a big difference!"
	@echo
	@echo "3) This difference does not exist for the unconstrained fit in MLF (none_mlf) or using the snapshot (none):"
	@tail -n3 logs/frozen_none{_mlf,}.log | grep "==>" -A1
	@echo
	@echo "=> This is just a (successful) closure test."
	@echo
	@echo "4) What is the systematic uncertainty?"
	@echo "   Assuming (in quadrature notation) that tot = stat + syst, we can estimate syst = tot - stat:"
	@echo
	@echo "4.1) The (post-fit) full uncertainty is given by (frozen_none):"
	@grep "Best fit r" logs/frozen_none.log
	@echo
	@echo "4.2) The (post-fit) stat-only uncertainty is given by (frozen_all):"
	@grep "Best fit r" logs/frozen_all.log
	@echo
	@echo "4.3) The difference (the syst) is:"
	@grep -H "Best fit r" logs/frozen_{none,all}.log
	@./MLFUncDiff.py higgsCombine_frozen_{none,all}.MaxLikelihoodFit.mH120.root | grep -v "Loading FW Lite setup."
	@echo
	@echo "5) Let's get the exp. syst. uncertainty."
	@echo "   Assuming (in quadrature notation) that tot = stat + exp + theo, and given the way we defined the groups of nuisances, there are two ways:"
	@echo "   - (exp) = (tot) - (stat+theo), or"
	@echo "   - (exp) = (tot-theo) - (stat)."
	@echo "  From the fits, we have all the terms above in brackets:"
	@echo
	@echo "5.1) (exp) = (tot) - (stat+theo), which translates to (exp) = [frozen_none] - [frozen_exp]"
	@grep -H "Best fit r" logs/frozen_{none,exp}.log
	@./MLFUncDiff.py higgsCombine_frozen_{none,exp}.MaxLikelihoodFit.mH120.root | grep -v "Loading FW Lite setup."
	@echo
	@echo "5.2) (exp) = (tot-theo) - (stat), which translates to (exp) = [frozen_theo] - [frozen_all]"
	@grep -H "Best fit r" logs/frozen_{theo,all}.log
	@./MLFUncDiff.py higgsCombine_frozen_{theo,all}.MaxLikelihoodFit.mH120.root | grep -v "Loading FW Lite setup."
	@echo
	@echo "=> These are very similar as expected."
	@echo
	@echo "6) Let's get the theo. syst. uncertainty."
	@echo "   Assuming (in quadrature notation) that tot = stat + exp + theo, and given the way we defined the groups of nuisances, there are two ways:"
	@echo "   - (theo) = (tot) - (stat+exp), or"
	@echo "   - (theo) = (tot-exp) - (stat)."
	@echo "  From the fits, we have all the terms above in brackets:"
	@echo
	@echo "6.1) (theo) = (tot) - (stat+exp), which translates to (theo) = [frozen_none] - [frozen_theo]"
	@grep -H "Best fit r" logs/frozen_{none,theo}.log
	@./MLFUncDiff.py higgsCombine_frozen_{none,theo}.MaxLikelihoodFit.mH120.root | grep -v "Loading FW Lite setup."
	@echo
	@echo "6.2) (theo) = (tot-exp) - (stat), which translates to (theo) = [frozen_exp] - [frozen_all]"
	@grep -H "Best fit r" logs/frozen_{exp,all}.log
	@./MLFUncDiff.py higgsCombine_frozen_{exp,all}.MaxLikelihoodFit.mH120.root | grep -v "Loading FW Lite setup."
	@echo
	@echo "=> These are also very similar as expected."
	@echo
	@echo "7) Let's look at partial effects of groups of nuisances."
	@echo "   For this case we will proceed with (group) = (tot) - (stat+group)"
	@echo
	@echo "7.1)"
	@grep "group_bg:" logs/myanalysis.T2W.log
	@grep -H "Best fit r" logs/frozen_{none,bg}.log
	@./MLFUncDiff.py higgsCombine_frozen_{none,bg}.MaxLikelihoodFit.mH120.root | grep -v "Loading FW Lite setup."
	@echo
	@echo "7.2)"
	@grep "group_bgtheo:" logs/myanalysis.T2W.log
	@grep -H "Best fit r" logs/frozen_{none,bgtheo}.log
	@./MLFUncDiff.py higgsCombine_frozen_{none,bgtheo}.MaxLikelihoodFit.mH120.root | grep -v "Loading FW Lite setup."
	@echo
	@echo "7.2)"
	@grep "group_bgexp:" logs/myanalysis.T2W.log
	@grep -H "Best fit r" logs/frozen_{none,bgexp}.log
	@./MLFUncDiff.py higgsCombine_frozen_{none,bgexp}.MaxLikelihoodFit.mH120.root | grep -v "Loading FW Lite setup."
	@echo
	@echo "7.4)"
	@grep "group_sigtheo:" logs/myanalysis.T2W.log
	@grep -H "Best fit r" logs/frozen_{none,sigtheo}.log
	@./MLFUncDiff.py higgsCombine_frozen_{none,sigtheo}.MaxLikelihoodFit.mH120.root | grep -v "Loading FW Lite setup."
	@echo
	@echo "8) What can we learn?"
	@echo "   With these, we see that this analysis would be syst-limited (4.3 vs 4.2),"
	@echo "   with the most important systematics being from theory (6 vs 5),"
	@echo "   namely the theory part on the bg (7.2 vs 5 and 7.3)."
	@echo "   To learn more you need to look at the individual impacts <insert here next lecture>."

allfits: $(patsubst %,higgsCombine_%.MaxLikelihoodFit.mH120.root,$(targets))

%.ws.root: %.dc.txt
	@mkdir -p logs
	@echo "Creating binary datacard:"
	text2workspace.py -o $@ -v9 $< &> logs/$*.T2W.log

%.likelihood.png: %.ws.root.dot
	dot -Tpng $< -o $@

fitter_opts = --robustFit 1 --minimizerTolerance 0.00000001 --minimizerStrategy 2

# Fit and save best-fit snapshot
higgsCombine_%_bestfit.MultiDimFit.mH120.root: %.ws.root
	@echo "Running and saving the best-fit snapshot:"
	combine -M MultiDimFit $(fitter_opts) --saveWorkspace -n _$*_bestfit -v9 $< &> logs/$*_bestfit.log

# Fit freezing each group (stat + other groups)
$(patsubst %,higgsCombine_frozen_%.MaxLikelihoodFit.mH120.root,$(nuisgroups)): higgsCombine_frozen_%.MaxLikelihoodFit.mH120.root : higgsCombine_$(cardname)_bestfit.MultiDimFit.mH120.root
	@echo "Profiling with [$*] group of nuisances frozen to best-fit values:"
	combine -M MaxLikelihoodFit --freezeNuisanceGroups $*        $(fitter_opts) -d $< -w w --snapshotName "MultiDimFit" -n _frozen_$* -v9 &> logs/frozen_$*.log

#Fit freezing everything (stat-only)
higgsCombine_frozen_all.MaxLikelihoodFit.mH120.root: higgsCombine_frozen_%.MaxLikelihoodFit.mH120.root : higgsCombine_$(cardname)_bestfit.MultiDimFit.mH120.root
	@echo "Profiling with all nuisances frozen to best-fit values:"
	#TODO --profilingMode none is not working?
	combine -M MaxLikelihoodFit --freezeNuisanceGroups exp,theo  $(fitter_opts) -d $< -w w --snapshotName "MultiDimFit" -n _frozen_$* -v9 &> logs/frozen_$*.log

# Fit freezing nothing (full uncertainty)
higgsCombine_frozen_none.MaxLikelihoodFit.mH120.root: higgsCombine_frozen_%.MaxLikelihoodFit.mH120.root : higgsCombine_$(cardname)_bestfit.MultiDimFit.mH120.root
	@echo "Profiling with no nuisances frozen:"
	combine -M MaxLikelihoodFit                                  $(fitter_opts) -d $< -w w --snapshotName "MultiDimFit" -n _frozen_$* -v9 &> logs/frozen_$*.log

# Fit with everything floating directly (should be same as previous line)
higgsCombine_frozen_none_mlf.MaxLikelihoodFit.mH120.root: higgsCombine_frozen_%.MaxLikelihoodFit.mH120.root : $(cardname).ws.root
	@echo "Profiling with no nuisances frozen (and not using the snapshot):"
	combine -M MaxLikelihoodFit                                  $(fitter_opts) $<                                      -n _frozen_$* -v9 &> logs/frozen_$*.log

# Fit fixing nuisances to pre-fit
higgsCombine_frozen_all_mlf.MaxLikelihoodFit.mH120.root: higgsCombine_frozen_%.MaxLikelihoodFit.mH120.root : $(cardname).ws.root
	@echo "Profiling with all nuisances frozen (and not using the snapshot, so frozen to default values):"
	#TODO --profilingMode none is not working?
	combine -M MaxLikelihoodFit --freezeNuisanceGroups exp,theo  $(fitter_opts) $<                                      -n _frozen_$* -v9 &> logs/frozen_$*.log


.PHONY: clean
clean:
	rm -f $(patsubst %,higgsCombine_%.MaxLikelihoodFit.mH120.root,$(targets))
	rm -f $(patsubst %,mlfit_%.root,$(targets))
	rm -f higgsCombine_$(cardname)_bestfit.*.root
	rm -f $(cardname).ws.root{,.dot}
	rm -rf logs
distclean: clean
	rm -f *~

.SECONDARY: $(wildcard *.ws.root)
