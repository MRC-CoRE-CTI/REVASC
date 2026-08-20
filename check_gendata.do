/*
check_gendata.do
Take data sets from R and estimate 5-year survival from them
IW updated 20aug2026
*/

* User-specific setup
cd "C:/ian/git/REVASC"
* End of user-specific setup

stset survt, failure(status) scale(12)

sts if subgroup==1, by(rand) name(subgroup1, replace)
sts if rand==0, by(subgroup) name(controls, replace)

sts list if rand==0, risktable(5) noshow by(subgroup) failure

* should match Weiqi's numbers: .623,.490,.398,.609,.392,.720,.689
