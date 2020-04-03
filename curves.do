clear all
run "c:\mijn documenten\projecten\stata\abm\covid\covid.mata"

set seed 123456
drop _all
graph drop _all
local n = 10000
mata:
model = covid()
model.N_nodes(`n')
model.degree(10)
model.pr_weak(0.05)
model.R0(2.5)
model.infect_dur(14)
model.N_index(10)

model.tdim(730)
model.run("quietly")
model.toStata("inf")
end

gen S = (`n' - cumul)/`n'*100
gen I = curent/`n'*100
gen R = (cumul-curent)/`n'*100

twoway line S R I t, ///
    legend(order(1 3 2) rows(1)) scheme(s1rcolor) ///
	ylab(,angle(0)) ///
	ytitle("percentage of population") ///
	xtitle("time") lw(*2 ..)           ///
	name(p5) title(5% long distance connections)

frame create r05
frame change r05
mata: model.toStata("R")
rename R R05
collapse R05, by(t_start)
frame change default
	
drop _all
mata
model.pr_weak(0.15)
model.run("quietly")
model.toStata("inf")
end

gen S = (`n' - cumul)/`n'*100
gen I = curent/`n'*100
gen R = (cumul-curent)/`n'*100

twoway line S R I t, ///
    legend(order(1 3 2) rows(1)) scheme(s1rcolor) ///
	ylab(,angle(0)) ///
	ytitle("percentage of population") ///
	xtitle("time") lw(*2 ..)           ///
	name(p15) title(15% long distance connections)

frame create r15
frame change r15
mata: model.toStata("R")
rename R R15
collapse R15, by(t_start)
frame change default
	
	
drop _all
mata
model.pr_weak(0.25)
model.run("quietly")
model.toStata("inf")
end

gen S = (`n' - cumul)/`n'*100
gen I = curent/`n'*100
gen R = (cumul-curent)/`n'*100

twoway line S R I t, ///
    legend(order(1 3 2) rows(1)) scheme(s1rcolor) ///
	ylab(,angle(0)) ///
	ytitle("percentage of population") ///
	xtitle("time") lw(*2 ..)           ///
	name(p25) title(25% long distance connections)	

frame create r25
frame change r25
mata: model.toStata("R")
rename R R25
collapse R25, by(t_start)
	
frame change r05	
frlink 1:1 t_start, frame(r15)
frget R15, from(r15)
frlink 1:1 t_start, frame(r25)
frget R25, from(r25)

twoway lowess R05 t_start, lw(*2) || ///
       lowess R15 t_start, lw(*2) || ///
	   lowess R25 t_start, lw(*2)    ///
	   scheme(s1rcolor)              ///
	   legend(order(- "% long distance" "connections" 1 "5%" 2 "15%" 3 "25%") rows(1)) ///
	   ytitle("average number of infections" "(smoothed)") ///
	   xtitle("day infected") name(R) ///
	   yline(1, lcolor(white) lpattern(dash) lw(*.5))
exit