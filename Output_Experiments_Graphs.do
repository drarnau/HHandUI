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

// Specify variables and experiment
global myvars = "valuevf consumption"
global qs = 4 // Number of quantiles for conditional wealth and income
global myexp = "benchmark"
global nexp = 32 // Number of experiments
global nbm = 8 // Number of benchmark experiment

// Create empty data set where to store results for all experiments
// Each observation is one experiment
cd $dir_work
foreach v in "b_0" "b_bar" "mu" {
	gen `v' = .
	}

forval t = 1/3 {
foreach v of global myvars {
	foreach d in "mean" "p50"  {
		gen `v'_`d'_HH`t' = .
	}
	foreach cvar in "wealth" "income" {
		forval myq = 1/$qs {
			gen `v'_`cvar'q`myq'_HH`t' = .
		}
	}
}
}

save output_$myexp, replace
clear

// Iterate over experiments
forval e = 1/$nexp {
	local dir_aux = "$dir_work" + "$myexp" + "`e'" + "/"

	cd `dir_aux'
	// Load experiment info
	import delim using UI.txt

	global b_0 = v1[1]
	global b_bar = v1[2]
	global mu = v1[3]
	clear

	// Load aggregate variables
	import delim using aggvars.txt
	global r = ((v1[1]+1)^12)-1
	global tau = v1[2]
	global KLratio = v1[3]
	global average_z = v1[4]
	clear

	// Load model otuput
	import delim using model.txt
	local linia = 1
	forval m = 0/1 {
	forval s = 1/2 {
			forval td = 1/3 {
			forval tw = 1/3 {
				global TR`td'`tw'_m`m'_s`s' = v1[`linia']
				local linia = `linia' + 1
			}
			}
			foreach r in "ER" "UR" {
				global `r'_m`m'_s`s' = v1[`linia']
				local linia = `linia' + 1
			}
	}
	}
	clear

	// Read csv files
	forval f = 1/3 {
		preserve
		import delimited simulation_`f'.csv, varnames(1) clear
		gen HHtype = `f'
		gen income = labourincome  + benefitsreceived
		save temp, replace
		restore
		append using temp
		}
	erase temp.dta

	// Label HHtype variable
	label var HHtype "Household type"
	label define lab_hhtype  1 "Single Men" 2 "Single Women" 3 "Married Household"
	lab val HHtype lab_hhtype

	// Compute myvars moments
	forval t = 1/3 {
	foreach v of global myvars {
		su `v' if HHtype == `t', d
		foreach d in "mean" "p50" {
			gen `v'_`d'_HH`t' = `r(`d')'
		}
		foreach cvar in "wealth" "income" {
			xtile aux = `cvar', n(4)
			forval myq = 1/$qs {
				su `v' if HHtype == `t' & aux == `myq'
				gen `v'_`cvar'q`myq'_HH`t' = `r(mean)'
			}
			drop aux
		}
	}
	}

	// Compute gini
	fastgini wealth
	gen gwealth_HHall = `r(gini)'
	forval t = 1/3 {
		fastgini wealth if HHtype == `t'
		gen gwealth_HH`t' = `r(gini)'
	}

	fastgini wealth if HHtype <= 2
	gen gwealth_HHsingles = `r(gini)'

	// Compute ratio assets
	su wealth if HHtype <= 2
	local singles = `r(sum)'

	su wealth if HHtype == 3
	gen rwealth_HHall = `r(sum)' / `singles'

	// Compute UR and ER conditional on spouse state
	forval s = 1/2 {
		if (`s' == 1) {
			local main = "male"
			local sp = "female"
		}
		else if (`s' == 2) {
			local main = "female"
		  local sp = "male"
		}

		// Auxiliary spouse employed variable
		gen auxsp = 1 if status`sp'today == 1
		replace auxsp = 0 if status`sp'today != 1


		// Employment rate
		gen aux = 1 if status`main'today == 1
		replace aux = 0 if aux == .

		forval spemp = 0/1 {
			su aux if auxsp == `spemp'
			gen ER_s`s'_sp`spemp' = `r(mean)'
		}
		drop aux

		// Unemployment rate
		gen aux = 1 if status`main'today == 2
		replace aux = 0 if aux == . & status`main'today == 1

		forval spemp = 0/1 {
			su aux if auxsp == `spemp'
			gen UR_s`s'_sp`spemp' = `r(mean)'
		}
		drop aux

		// Transitions
		forval td = 1/3 {
		forval tw = 1/3 {
			gen aux = 1 if status`main'today == `td' & status`main'tomorrow == `tw'
			replace aux = 0 if aux == . & status`main'today == `td'

			forval spemp = 0/1 {
				su aux if auxsp == `spemp'
				gen TR`td'`tw'_s`s'_sp`spemp' = `r(mean)'
			}
			drop aux
		}
		}
		drop auxsp
	}

	// Keep just one observation and relevant variables
	keep if _n == 1
	keep *_HH* ER* UR* TR*

	// Add variables read in txt files
	foreach v in "b_0" "b_bar" "mu" "r" "tau" "KLratio" "average_z" {
		gen `v' = ${`v'}
		}


	forval m = 0/1 {
	forval s = 1/2 {
		forval td = 1/3 {
		forval tw = 1/3 {
			gen TR`td'`tw'_m`m'_s`s' = ${TR`td'`tw'_m`m'_s`s'}
		}
		}
		foreach r in "ER" "UR" {
			gen `r'_m`m'_s`s' = ${`r'_m`m'_s`s'}
		}
	}
	}

	// Append to output.dat
	local f = "$dir_work" + "output_" + "$myexp" + ".dta"
	ap using `f'
	cd $dir_work
	save output_$myexp, replace
	clear
	}


// Open output dataset
cd $dir_work
use output_$myexp.dta

// Normalise all variables as percentage change with respect to benchmark
sort b_0
foreach v of var *_HH* tau* KLratio* TR* ER* UR* {
	local bm = `v'[$nbm]
	gen wrtb_`v' = ((`v' - `bm')/`bm')*100
	}

// Rates and probabilities as percentage
foreach v of var TR* ER* UR* {
	replace `v' = `v' * 100
}
// Label all variables
foreach v of var *_HH1 {
	label var `v' "Single Men"
}
foreach v of var *_HH2 {
	label var `v' "Single Women"
}
foreach v of var *_HH3 {
	label var `v' "Married"
}

label var b_0 "Replacement ratio"
label var b_bar "Benefits cap"
label var mu "Duration"
label var tau "Income tax"
label var wrtb_tau "Income tax"
label var KLratio "K to L ratio"
label var wrtb_KLratio "K to L ratio"
label var average_z "Average z employed"
label var r "Interest rate"

foreach v of var *m0_s1 {
	label var `v' "Single Men"
}
foreach v of var *m0_s2 {
	label var `v' "Single Women"
}
foreach v of var *m1_s1 {
	label var `v' "Married Men"
}
foreach v of var *m1_s2 {
	label var `v' "Married Women"
}

foreach v of var *all* {
	label var `v' "All"
}

foreach v of var *s1_sp0 {
	label var `v' "Male, spouse not emp."
}
foreach v of var *s1_sp1 {
	label var `v' "Male, spouse emp."
}
foreach v of var *s2_sp0 {
	label var `v' "Female, spouse not emp."
}
foreach v of var *s2_sp1 {
	label var `v' "Female, spouse emp."
}


// Plot all variables
set scheme plotplainblind

foreach tp in "" "wrtb_" {
foreach v of global myvars {
	foreach s in "mean" "p50" {
		local file_aux = "$dir_output" + "$myexp" + "_" + "`tp'" + "`v'" + "_" + "`s'" + ".eps"
		forval t = 1/3 {
			qui su `tp'`v'_`s'_HH`t'
			local max`t' = `r(max)'
			local min`t' = `r(min)'
			}
			local ub = ceil(max(`max1',`max2',`max3'))
			local lb = floor(min(`min1',`min2',`min3'))

		if abs(`ub') + abs(`lb') > 2  {
			twoway (scatter `tp'`v'_`s'_HH1 `tp'`v'_`s'_HH2 `tp'`v'_`s'_HH3 b_0) ///
			(mspline `tp'`v'_`s'_HH1 b_0, lcolor(black)) ///
			(mspline `tp'`v'_`s'_HH2 b_0, lcolor(gray)) ///
			(mspline `tp'`v'_`s'_HH3 b_0, lcolor(ltblue)), ymtick(`lb'(1)`ub') legend(order(1 2 3))
			}
		else {
			twoway (scatter `tp'`v'_`s'_HH1 `tp'`v'_`s'_HH2 `tp'`v'_`s'_HH3 b_0) ///
			(mspline `tp'`v'_`s'_HH1 b_0, lcolor(black)) ///
			(mspline `tp'`v'_`s'_HH2 b_0, lcolor(gray)) ///
			(mspline `tp'`v'_`s'_HH3 b_0, lcolor(ltblue)),  legend(order(1 2 3))
			}
		gr export `file_aux', replace
	}


	// conditional graphs
	foreach cvar in "income" "wealth" {
	forval myq = 1/$qs {
		local file_aux = "$dir_output" + "$myexp" + "_" + "`tp'" + ///
											"`v'" + "_" + "`cvar'" + "q`myq'" + ".eps"
		forval t = 1/3 {
			qui su `tp'`v'_`cvar'q`myq'_HH`t'
			local max`t' = `r(max)'
			local min`t' = `r(min)'
			}
			local ub = ceil(max(`max1',`max2',`max3'))
			local lb = floor(min(`min1',`min2',`min3'))

		if abs(`ub') + abs(`lb') > 2  {
			twoway (scatter `tp'`v'_`cvar'q`myq'_HH1 `tp'`v'_`cvar'q`myq'_HH2 `tp'`v'_`cvar'q`myq'_HH3 b_0) ///
			(mspline `tp'`v'_`cvar'q`myq'_HH1 b_0, lcolor(black)) ///
			(mspline `tp'`v'_`cvar'q`myq'_HH2 b_0, lcolor(gray)) ///
			(mspline `tp'`v'_`cvar'q`myq'_HH3 b_0, lcolor(ltblue)), ymtick(`lb'(1)`ub') legend(order(1 2 3))
			}
		else {
			twoway (scatter `tp'`v'_`cvar'q`myq'_HH1 `tp'`v'_`cvar'q`myq'_HH2 `tp'`v'_`cvar'q`myq'_HH3 b_0) ///
			(mspline `tp'`v'_`cvar'q`myq'_HH1 b_0, lcolor(black)) ///
			(mspline `tp'`v'_`cvar'q`myq'_HH2 b_0, lcolor(gray)) ///
			(mspline `tp'`v'_`cvar'q`myq'_HH3 b_0, lcolor(ltblue)),  legend(order(1 2 3))
			}
		gr export `file_aux', replace
	}
	}
}
}

// TR, ER, and UR
foreach tp in "" "wrtb_" {
foreach v in 	"TR11" "TR12" "TR13" ///
							"TR21" "TR22" "TR23" ///
							"TR31" "TR32" "TR33" ///
							"ER" "UR" {
	local file_aux = "$dir_output" + "$myexp" + "_" + "`tp'" + "`v'" + ".eps"
	twoway (scatter `tp'`v'_m0_s1 `tp'`v'_m0_s2 `tp'`v'_m1_s1 `tp'`v'_m1_s2 b_0) ///
	(mspline `tp'`v'_m0_s1 b_0, lcolor(black)) ///
	(mspline `tp'`v'_m0_s2 b_0, lcolor(gray)) ///
	(mspline `tp'`v'_m1_s1 b_0, lcolor(ltblue)) ///
	(mspline `tp'`v'_m1_s2 b_0, lcolor(green)),  legend(order(1 2 3 4))
	gr export `file_aux', replace
}
}

// TR, ER, and UR conditional spouse
foreach tp in "" "wrtb_" {
foreach v in 	"TR11" "TR12" "TR13" ///
							"TR21" "TR22" "TR23" ///
							"TR31" "TR32" "TR33" ///
							"ER" "UR" {
	local file_aux = "$dir_output" + "$myexp" + "_" + "`tp'" + "`v'" + "_cspouse" + ".eps"
	twoway (scatter `tp'`v'_s1_sp0 `tp'`v'_s1_sp1 `tp'`v'_s2_sp0 `tp'`v'_s2_sp1 b_0) ///
	(mspline `tp'`v'_s1_sp0 b_0, lcolor(black)) ///
	(mspline `tp'`v'_s1_sp1 b_0, lcolor(gray)) ///
	(mspline `tp'`v'_s2_sp0 b_0, lcolor(ltblue)) ///
	(mspline `tp'`v'_s2_sp1 b_0, lcolor(green)),  legend(order(1 2 3 4))
	gr export `file_aux', replace
}
}

// Variables graphed alone
foreach tp in "" "wrtb_" {
foreach v in "rwealth_HHall" "KLratio" "tau" {
	local file_aux = "$dir_output" + "$myexp" + "_" + "`tp'" + "`v'" + ".eps"
	twoway (scatter `tp'`v' b_0) ///
	(mspline `tp'`v' b_0, lcolor(black)),  legend(off)
	gr export `file_aux', replace
}
}

// Gini
foreach tp in "" "wrtb_" {
	local file_aux = "$dir_output" + "$myexp" + "_" + "`tp'" + "gini" + ".eps"
	twoway (scatter `tp'gwealth_HH1 `tp'gwealth_HH2 `tp'gwealth_HH3 `tp'gwealth_HHall b_0) ///
	(mspline `tp'gwealth_HH1 b_0, lcolor(black)) ///
	(mspline `tp'gwealth_HH2 b_0, lcolor(gray)) ///
	(mspline `tp'gwealth_HH3 b_0, lcolor(ltblue)) ///
	(mspline `tp'gwealth_HHall b_0, lcolor(green)),  legend(order(1 2 3 4))
	gr export `file_aux', replace
}


// Erase output file
// erase output_$myexp.dta
