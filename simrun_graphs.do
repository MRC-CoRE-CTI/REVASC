* simrun_graphs.do
* edited 20apr2026

* User-specific setup
cd "C:/ian/git/REVASC"
* End of user-specific setup

import excel "simrun_results_pretty.xlsx", sheet("run 2") firstrow case(lower) clear
replace dgm=dgm[_n-1] if mi(dgm)
replace n_tot=n_tot[_n-1] if mi(n_tot)
drop if mi(dgm)
drop if mi(analysis)
destring fwer avepower wtavepower, replace force
destring n_reps, replace force
label var n_tot "Total sample size"
label var wtavepower "Wt average"
label var analysis "Method of analysis"

* Graph powers under global alternative
line wtavepower cl* n_tot if dgm==dgm[3], by(analysis, imargin(medlarge)) legend(row(2)) lcol(black) lpattern(dash) ytitle(Power) name(simrun_power_galt, replace) saving(simrun_power_galt, replace) xsize(6) ysize(4.5) yli(.9) xtitle(,size(large)) ytitle(,size(large)) xlabel(,labsize(large)) ylabel(,labsize(large)) legend(size(large)) subtitle(,size(large))

* Zoom in on sample size needed for 90% wt average power under borrowing
* .892, .908 are MCCI for .9 with 5000 reps
solve wtavepower n_tot if dgm==dgm[3] & analysis=="borrow", value(.9)
solve wtavepower n_tot if dgm==dgm[3] & analysis=="borrow", value(.892)
solve wtavepower n_tot if dgm==dgm[3] & analysis=="borrow", value(.908)
gen mcse = sqrt(wtavepower*(1-wtavepower)/n_reps)
cigraph8 wtavepower mcse n_tot if dgm==dgm[3] & analysis=="borrow", yline(.9) s(O) col(black) pat(dash) yti(Average power) name(simrun_power_zoom,replace) saving(simrun_power_zoom,replace)
