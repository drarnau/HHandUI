/*
Choi & Valladares-Esteban (2018)

This codes reads Monte-Carlo simulation data from the model and generates results
*/

cls
clear all
set more off, perm
set trace off

// Output dir
local dir_output = "/home/arnau/Dropbox/Choi_Valladares_2015/QEresubmission/"

// Working directory
cd "/home/arnau/Dropbox/Choi_Valladares_2015/QEresubmission/code/HHandUI/"

// Auxiliary name
local aux_name = "No_UI_"

// Read csv files
forval f = 1/3 {
	preserve
	import delimited `aux_name'simulation_`f'.csv, varnames(1) clear
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

// Histagrams
	// All
	hist wealth, width(50)
	graph export "`dir_output'`aux_name'wealth_total.png", replace
	hist valuevf, width(5) //ysc(r(0 0.02)) ylabel(0(0.005)0.02)
	graph export "`dir_output'`aux_name'valuevf_total.png", replace

	// Married and single
	twoway 	(hist wealth if married==1, width(50) color(green)) ///
		(hist wealth if married==0, width(50) ///
		fcolor(none) lcolor(black)), legend(order(1 "Married" 2 "Single" ))
	graph export "`dir_output'`aux_name'wealth_bymaritalstatus.png", replace


	twoway 	(hist valuevf if married==1, width(5) color(green)) ///
		(hist valuevf if married==0, width(5) ///
		fcolor(none) lcolor(black)), legend(order(1 "Married" 2 "Single" ))
	graph export "`dir_output'`aux_name'valuevf_bymaritalstatus.png", replace

// Added worker effect
gen awe = 1 if statusfemaletoday == 3 & statusfemaletomorrow == 2
replace awe = 1 if statusfemaletoday == 3 & statusfemaletomorrow == 1
// replace awe = 0 if married == 1 & statusfemaletoday == 3 & awe == .
replace awe = 0 if married == 1 & awe == .

gen EU = 1 if statusmaletoday == 1 & statusmaletomorrow == 2
// replace EU = 0 if married == 1 & statusmaletoday == 1 & EU == .
replace EU = 0 if married == 1 & EU == .

reg awe EU


/*
// Summary with details
	// All
	su wealth, d
	su valuevf, d

	// Married and single
	forval m = 0/1 {
		su wealth if married == `m', d
		su valuevf if married == `m', d
		}
*/

// Wealth distribution single vs. married
forval m = 0/1 {
	su wealth if married == `m'
	local mean_`m' = `r(mean)'
	local cv_`m' = `r(sd)'/`r(mean)'
	fastgini wealth if married == `m'
	local gini_`m' = r(gini)
	}
	local rmean = `mean_1'/`mean_0'

	local rcv = `cv_1'/`cv_0'
	local rgini = `gini_1'/`gini_0'

	foreach z in "rmean" "rcv" "rgini" {
		if ``z'' > 1 {
			local `z' = substr("``z''",1,6)
			}
		else {
			local `z' = substr("``z''",1,5)
			local `z' = subinstr("``z''",".","0.",1)
			}
		}
forvalues z=1/1 { // necessary not to have the commands in the tex file
	qui log using "`dir_output'`aux_name'ratios_wealth.tex", text replace
	set linesize 255
	display "\begin{centering}"
	display "\begin{tabular}{lrr}"
	display " & Data & Model\tabularnewline"
	display "\hline"
	display "\hline"
	display "Average 			& 2.8481 	& `rmean' 	\tabularnewline"
	display "Coefficient of Variation 	& 0.7575 	& `rcv' 	\tabularnewline"
	display "Gini Index 			& 0.9512 	& `rgini' 	\tabularnewline"
	display "\hline"
	display "\end{tabular}"
	display "\end{centering}"
	qui log close
	} // End loop z

// Value VF distribution single vs. married
	forval m = 0/1 {
		su valuevf if married == `m'
		local mean_`m' = `r(mean)'
		local cv_`m' = `r(sd)'/`r(mean)'
		fastgini wealth if married == `m'
		local gini_`m' = r(gini)

		foreach z in "mean_`m'" "cv_`m'" "gini_`m'" {
			if ``z'' > 1 {
				local `z' = substr("``z''",1,6)
				}
			else {
				local `z' = substr("``z''",1,5)
				local `z' = subinstr("``z''",".","0.",1)
				}
			}
	} // m

	forvalues z=1/1 { // necessary not to have the commands in the tex file
		qui log using "`dir_output'`aux_name'descriptives_valuevf.tex", text replace
		set linesize 255
		display "\begin{centering}"
		display "\begin{tabular}{lrr}"
		display " & Single & Married\tabularnewline"
		display "\hline"
		display "\hline"
		display "Average 			& `mean_0' 	& `mean_1' 	\tabularnewline"
		display "Coefficient of Variation 	& `cv_0' 	& `cv_1' 	\tabularnewline"
		display "Gini Index 			& `gini_0' 	& `gini_1' 	\tabularnewline"
		display "\hline"
		display "\end{tabular}"
		display "\end{centering}"
		qui log close
		} // End loop z

// Transitions by quintiles
xtile quintile = wealth, nquantiles(5)

// Quaintile dummies
forval q = 1/5 {
	gen quint_`q' = 1 if quintile == `q'
	}
gen quint_0 = 1 // Include the whole sample

forval q = 0/5 { // Iterate over all possible quintile groups
	forval td = 1/3 { // Today
	forval tw = 1/3 { // Tomorrow
		qui gen aux = .
		foreach g in "" "male" "female" {
			qui replace aux = 1 if status`g'today == `td' & status`g'tomorrow == `tw' & quint_`q' == 1
			qui replace aux = 0 if status`g'today == `td' & aux == . & quint_`q' == 1

			}
		qui su aux
		local trans_`td'`tw'_q`q' = `r(mean)'
		drop aux
		}
		}
	}

forval q = 1/5 {
forval td = 1/3 { // Today
forval tw = 1/3 { // Tomorrow
	local t`td'`tw'_q`q' = `trans_`td'`tw'_q`q'' / `trans_`td'`tw'_q0'
	if `t`td'`tw'_q`q'' > 1 {
		local t`td'`tw'_q`q' = substr("`t`td'`tw'_q`q''",1,4)
		}
	else if `t`td'`tw'_q`q'' == 1 {
		local t`td'`tw'_q`q' = "1.00"
		}
	else {
		local t`td'`tw'_q`q' = substr("`t`td'`tw'_q`q''",1,3)
		local t`td'`tw'_q`q' = subinstr("`t`td'`tw'_q`q''",".","0.",1)
		}

	display `t`td'`tw'_q`q''
	}
	}
	}

forvalues z=1/1 { // necessary not to have the commands in the tex file
	qui log using "`dir_output'`aux_name'trans_quintiles.tex", text replace
	set linesize 255
	display "\begin{centering}"
	display "\begin{tabular}{cccccccccccc}"
	display "& \multicolumn{5}{c}{Data} &  & \multicolumn{5}{c}{Model}\tabularnewline"
	display "\cline{2-6} \cline{8-12}"
	display " & Q1 & Q2 & Q3 & Q4 & Q5 &  & Q1 & Q2 & Q3 & Q4 & Q5\tabularnewline"
	display "\hline"
	display "\hline"
	display "EE & 0.95 & 0.99 & 1.01 & 1.02 & 1.02 & EE & `t11_q1' & `t11_q2' & `t11_q3' & `t11_q4' & `t11_q5' \tabularnewline"
	display "EU & 1.82 & 1.15 & 0.87 & 0.70 & 0.54 & EU & `t12_q1' & `t12_q2' & `t12_q3' & `t12_q4' & `t12_q5' \tabularnewline"
	display "EN & 1.15 & 0.92 & 0.91 & 0.97 & 1.09 & EN & `t13_q1' & `t13_q2' & `t13_q3' & `t13_q4' & `t13_q5' \tabularnewline"
	display "UE & 0.88 & 1.08 & 1.04 & 1.10 & 1.06 & UE & `t21_q1' & `t21_q2' & `t21_q3' & `t21_q4' & `t21_q5' \tabularnewline"
	display "UU & 1.07 & 0.97 & 0.97 & 0.95 & 0.91 & UU & `t22_q1' & `t22_q2' & `t22_q3' & `t22_q4' & `t22_q5' \tabularnewline"
	display "UN & 1.06 & 0.94 & 0.99 & 0.94 & 1.04 & UN & `t23_q1' & `t23_q2' & `t23_q3' & `t23_q4' & `t23_q5' \tabularnewline"
	display "NE & 1.05 & 1.34 & 0.99 & 0.89 & 0.84 & NE & `t31_q1' & `t31_q2' & `t31_q3' & `t31_q4' & `t31_q5' \tabularnewline"
	display "NU & 1.79 & 1.37 & 0.86 & 0.63 & 0.47 & NU & `t32_q1' & `t32_q2' & `t32_q3' & `t32_q4' & `t32_q5' \tabularnewline"
	display "NN & 0.97 & 0.96 & 1.01 & 1.02 & 1.03 & NN & `t33_q1' & `t33_q2' & `t33_q3' & `t33_q4' & `t33_q5' \tabularnewline"
	display "\hline"
	display "\end{tabular}"
	display "\end{centering}"
	qui log close
	} // End loop z
