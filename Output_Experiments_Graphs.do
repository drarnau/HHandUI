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
global myexp = "UIfromSky"
global nexp = 8 // Number of experiments
global nbm = 1 // Number of benchmark experiment

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

save output, replace
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

	// Compute VF moments
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
	// Keep just one observation and relevant variables
	keep if _n == 1
	keep *_HH*

	// Add variables read in txt files
	foreach v in "b_0" "b_bar" "mu" "r" "tau" "KLratio" "average_z" {
		gen `v' = ${`v'}
		}

	// Append to output.dat
	local f = "$dir_work" + "output.dta"
	ap using `f'
	cd $dir_work
	save output, replace
	clear
	}

// Open output dataset
cd $dir_work
use output.dta

// Normalise all variables as percentage change with respect to benchmark
sort b_0
	foreach v of var *_HH* {
		local bm = `v'[$nbm]
		replace `v' = ((`v' - `bm')/`bm')*100
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
label var KLratio "K to L ratio"
label var average_z "Average z employed"
label var r "Interest rate"


// Plot all variables
set scheme plotplainblind

foreach v of global myvars {
	foreach s in "mean" "p50" {
		local file_aux = "$dir_output" + "$myexp" + "_" + "`v'" + "_" + "`s'" + ".eps"
		forval t = 1/3 {
			qui su `v'_`s'_HH`t'
			local max`t' = `r(max)'
			local min`t' = `r(min)'
			}
			local ub = ceil(max(`max1',`max2',`max3'))
			local lb = floor(min(`min1',`min2',`min3'))

		if abs(`ub') + abs(`lb') > 2  {
			twoway (scatter `v'_`s'_HH1 `v'_`s'_HH2 `v'_`s'_HH3 b_0) ///
			(mspline `v'_`s'_HH1 b_0, lcolor(black)) ///
			(mspline `v'_`s'_HH2 b_0, lcolor(gray)) ///
			(mspline `v'_`s'_HH3 b_0, lcolor(ltblue)), ymtick(`lb'(1)`ub') legend(order(1 2 3))
			}
		else {
			twoway (scatter `v'_`s'_HH1 `v'_`s'_HH2 `v'_`s'_HH3 b_0) ///
			(mspline `v'_`s'_HH1 b_0, lcolor(black)) ///
			(mspline `v'_`s'_HH2 b_0, lcolor(gray)) ///
			(mspline `v'_`s'_HH3 b_0, lcolor(ltblue)),  legend(order(1 2 3))
			}
		gr export `file_aux', replace
	}
	// conditional graphs
	foreach cvar in "income" "wealth" {
	forval myq = 1/$qs {
		local file_aux = "$dir_output" + "$myexp" + "_" + "`v'" + "_" + "`cvar'" + "q`myq'" + ".eps"
		forval t = 1/3 {
			qui su `v'_`cvar'q`myq'_HH`t'
			local max`t' = `r(max)'
			local min`t' = `r(min)'
			}
			local ub = ceil(max(`max1',`max2',`max3'))
			local lb = floor(min(`min1',`min2',`min3'))

		if abs(`ub') + abs(`lb') > 2  {
			twoway (scatter `v'_`cvar'q`myq'_HH1 `v'_`cvar'q`myq'_HH2 `v'_`cvar'q`myq'_HH3 b_0) ///
			(mspline `v'_`cvar'q`myq'_HH1 b_0, lcolor(black)) ///
			(mspline `v'_`cvar'q`myq'_HH2 b_0, lcolor(gray)) ///
			(mspline `v'_`cvar'q`myq'_HH3 b_0, lcolor(ltblue)), ymtick(`lb'(1)`ub') legend(order(1 2 3))
			}
		else {
			twoway (scatter `v'_`cvar'q`myq'_HH1 `v'_`cvar'q`myq'_HH2 `v'_`cvar'q`myq'_HH3 b_0) ///
			(mspline `v'_`cvar'q`myq'_HH1 b_0, lcolor(black)) ///
			(mspline `v'_`cvar'q`myq'_HH2 b_0, lcolor(gray)) ///
			(mspline `v'_`cvar'q`myq'_HH3 b_0, lcolor(ltblue)),  legend(order(1 2 3))
			}
		gr export `file_aux', replace
	}
	}
}


// Erase output file
// erase output.dta
