/* 
Choi & Valladares-Esteban (2018)

This codes reads Monte-Carlo simulation data from the model and generates results
*/

cls
clear all
set more off, perm
set trace off

// Read csv files
forval f = 1/3 {
	preserve 
	import delimited simulation_`f'.csv, varnames(1) clear
	gen HHtype = `f'
	save temp, replace
	restore
	append using temp
	}
	
erase temp.dta

// Label HHtype variable
label var HHtype "Household type"
label define lab_hhtype  1 "Single Men" 2 "Single Women" 3 "Married Household"
lab val HHtype lab_hhtype 

// Married dummy
gen married = 1 if HHtype == 3
replace married = 0 if married == .
label var married "Married Dummy"
label define lab_married 0 "Single" 1 "Married"
lab val married lab_married

// Asset distribution single vs. married
tabstat wealth, by(married)
