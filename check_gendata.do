/*
check_gendata.do
Take data sets from R and estimate 5-year survival from them
*/

forvalues c=1/7 {
	use c:\temp\mydata`c', clear
	rename arr_`c' arr
	rename arm_`c' arm
	gen subgroup = `c'
	if `c'>1 append using check_gendata
	save check_gendata, replace
}
label var arr
label var arm
stset survt, fail(status)
sts if subgroup==1, by(arm) name(subgroup1, replace)
sts if arm==0, by(subgroup) name(controls, replace)

sts list if arm==0, risktable(60) noshow by(subgroup) failure

* should match Weiqi's numbers: .623,.490,.398,.609,.392,.720,.689
