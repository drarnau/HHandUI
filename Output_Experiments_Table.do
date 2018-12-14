/*
Choi & Valladares-Esteban (2018)

This codes reads Monte-Carlo simulation data from the model and generates results
*/

cls
clear all
set more off, perm
set trace off

// Output dir
global dir_output = "/home/arnau/Dropbox/Choi_Valladares_2015/QEresubmission/"

// Working directory
global dir_work = "/home/arnau/Dropbox/Choi_Valladares_2015/QEresubmission/code/HHandUI/"

// Specify variable and experiment details
global myvar = "valuevf"
global myexp = "benchmark"
global nexp = 16 // Number of experiments
global nbm = 14 // Number of benchmark experiment

// Iterate over experiments
forval e = 1/$nexp {
	local dir_aux = "$dir_work" + "$myexp" + "`e'" + "/"

	cd `dir_aux'

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

	forval t = 1/3 {
		su $myvar if HHtype == `t', d
		foreach v in "mean" "p10" "p25" "p50" "p75" "p90" {
			global `v'_e`e'_HH`t' = `r(`v')'
			}
		fastgini $myvar if HHtype == `t'
		global gini_e`e'_HH`t' = `r(gini)'
		}
	clear

	// Load experiment info
	import delim using UI.txt

	global b_0_e`e' = v1[1]
	global b_bar_e`e' = v1[2]
	global mu_e`e' = v1[3]
	clear

	// Load aggregate variables
	import delim using aggvars.txt
	global r_e`e' = v1[1] * 100
	global tau_e`e' = v1[2] * 100
	global KLratio_e`e' = v1[3]
	global average_z_e`e' = v1[4]
	clear
}

// Print results in a table
	// Normalise to benchmark and give right format
	forval e = 1/$nexp {
		foreach v in "mean" "p10" "p25" "p50" "p75" "p90" "gini" {
		forval t = 1/3 {
			local nbm = $nbm
			local `v'_e`e'_HH`t' = ((${`v'_e`e'_HH`t'} - ${`v'_e`nbm'_HH`t'})/${`v'_e`nbm'_HH`t'})*100
				if abs(``v'_e`e'_HH`t'') > 1 {
				local `v'_e`e'_HH`t' = substr("``v'_e`e'_HH`t''",1,6)
					}
				else {
					local `v'_e`e'_HH`t' = substr("``v'_e`e'_HH`t''",1,5)
					local `v'_e`e'_HH`t' = subinstr("``v'_e`e'_HH`t''",".","0.",1)
					}
				}
				}
		foreach z in "b_0_e`e'" "b_bar_e`e'" "mu_e`e'" {
			if ${`z'} > 1 {
				local `z' = substr("$`z'",1,4)
				}
			else {
				local `z' = substr("$`z'",1,3)
				local `z' = subinstr("``z''",".","0.",1)
				}
			}
		foreach z in "r_e`e'" "tau_e`e'" "KLratio_e`e'" "average_z_e`e'" {
			if ${`z'} > 1 {
				local `z' = substr("$`z'",1,6)
				}
			else {
				local `z' = substr("$`z'",1,5)
				local `z' = subinstr("``z''",".","0.",1)
				}
			}
		}
foreach v in "mean" "p10" "p25" "p50" "p75" "p90" "gini" {
forvalues z=1/1 { // necessary not to have the commands in the tex file
	local file_aux = "$dir_output" + "$myexp" + "_" + "$myvar" + "_" + "`v'" + ".tex"
	qui log using `file_aux', text replace
	set linesize 255
	di "\begin{centering}"
	di "\begin{tabular}{cccccrrr}"
	di "$""b_0$ & $\bar{b}$ & $\mu$ & Tax ($\tau$) & Interest rate & Single Males & Single Females & Married \tabularnewline"
	di "\hline"
	di "\hline"
	forval e = 1/$nexp {
		di "`b_0_e`e'' & `b_bar_e`e'' & `mu_e`e'' & `tau_e`e'' & `r_e`e'' & ``v'_e`e'_HH1' & ``v'_e`e'_HH2' & ``v'_e`e'_HH3' \tabularnewline"
		}
	di "\hline"
	di "\end{tabular}"
	di "\par \end{centering}"
	qui log close
	} // End loop z
	}
