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

// Specify variable and experiment
global myvar = "valuevf"
global myexp = "benchmark"

// Iterate over experiments
forval e = 1/8 {
	local dir_aux = "$dir_work" + "experiment" + "`e'" + "/"

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
	global r_e`e' = v1[1]
	global tau_e`e' = v1[2]
	global KLratio_e`e' = v1[3]
	global average_z_e`e' = v1[4]
	clear
}

// Print results in a table
	// Normalise to benchmark and give right format

	forval e = 1/9 {
		forval t = 1/3 {
			local vf_e`e'_HH`t' = ((${vf_e`e'_HH`t'} - ${vf_e3_HH`t'})/${vf_e3_HH`t'})*100
				if abs(`vf_e`e'_HH`t'') > 1 {
					local vf_e`e'_HH`t' = substr("`vf_e`e'_HH`t''",1,6)
					}
				else {
					local vf_e`e'_HH`t' = substr("`vf_e`e'_HH`t''",1,5)
					local vf_e`e'_HH`t' = subinstr("`vf_e`e'_HH`t''",".","0.",1)
					}
		}
		foreach z in "b_0_e`e'" "b_bar_e`e'" "mu_e`e'" {
			// di "`z'"
			if ${`z'} > 1 {
				local `z' = substr("$`z'",1,4)
				}
			else {
				local `z' = substr("$`z'",1,3)
				local `z' = subinstr("``z''",".","0.",1)
				}
		}
	}

forvalues z=1/1 { // necessary not to have the commands in the tex file
	cd $dir_output
	qui log using "table_UI_experiments.tex", text replace
	set linesize 255
	display "\begin{centering}"
	display "\begin{tabular}{cccrrr}"
	display "$""b_0$ & $\bar{b}$ & $\mu$ & Single Males & Single Females & Married \tabularnewline"
	display "\hline"
	display "\hline"
	forval e = 1/9 {
		di "`b_0_e`e'' & `b_bar_e`e'' & `mu_e`e'' & `vf_e`e'_HH1' & `vf_e`e'_HH2' & `vf_e`e'_HH3' \tabularnewline"
		}
	display "\hline"
	display "\end{tabular}"
	display "\par \end{centering}"
	qui log close
	} // End loop z
