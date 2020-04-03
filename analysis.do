clear all
run "c:\mijn documenten\projecten\stata\abm\covid\covid.mata"
set seed 123456

local n = 40
mata:
model = covid()
model.N_nodes(`n')
model.degree(4)
model.pr_weak(0.1)
model.R0(2.5)
model.infect_dur(5)
model.N_index(1)

model.tdim(10)
model.run("quietly")
model.toStata("nw")
end

gen x_ego = cos((ego-1)/`n'*2*_pi)*(1+mod(ego,2)*.2)
gen y_ego = sin((ego-1)/`n'*2*_pi)*(1+mod(ego,2)*.2)
gen x_alter = cos((alter-1)/`n'*2*_pi)*(1+mod(alter,2)*.2)
gen y_alter = sin((alter-1)/`n'*2*_pi)*(1+mod(alter,2)*.2)
    
bys ego (alter) : gen first = _n == 1
    
// the graphs
forvalues t = 1/10{
    twoway pcspike y_ego x_ego y_alter x_alter,                ///
                  lcolor(gs8) scheme(s1mono) ||                ///
        scatter y_ego x_ego if first & inf_t`t'==1 ,             ///
                   mcolor(red) msymbol(O) ||                   ///
        scatter y_ego x_ego if first & inf_t`t'==0,              ///
                   mcolor(black) msymbol(Oh)   ||                ///
        scatter y_ego x_ego if first & inf_t`t'==2,              ///
                   mcolor(gs8) msymbol(o)                   ///	
                   aspect(1) xscale(off) yscale(off)           ///
                   legend(order(3 "susceptible" 2 "infective" 4 "removed")) ///
                   name(t`t', replace) title(Time `t')     
}

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
gen R = cumul/`n'*100

twoway line S I R t, ///
    legend(rows(1)) scheme(s1rcolor) ///
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
gen R = cumul/`n'*100

twoway line S I R t, ///
    legend(rows(1)) scheme(s1rcolor) ///
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
gen R = cumul/`n'*100

twoway line S I R t, ///
    legend(rows(1)) scheme(s1rcolor) ///
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



mata:

result_weak = J(26,4,.)
trials = 10
res = J(trials,4,.)

model = covid()
model.N_nodes(10000)
model.degree(10)
model.pr_weak(0.1)
model.tdim(730)
model.R0(2)
model.infect_dur(14)
model.N_index(10)

for(j=0; j<=25; j++) {
	model.pr_weak(j/100)
	for(i=1; i<=trials; i++) {
		model.setup()
		model.run("quietly")
		res[i,.] = model.peak(), model.saved(), model.t_done(), model.R()
	}
	result_weak[j+1,.] = mean(res)   
}
result_weak
end

exit
// pr_weak = 0.05
mata:
model = covid()
model.N_nodes(10000)
model.degree(10)
model.pr_weak(0.05)
model.tdim(730)
model.R0(2)
model.infect_dur(14)
model.N_index(10)
model.setup()
for(i=1; i<=trials; i++) {
	model.run("quietly")
	res[i, .] = model.peak(), model.saved(), model.t_done()
}
mean(res)
end

// R0 = 2.5
mata:
model = covid()
model.N_nodes(10000)
model.degree(10)
model.pr_weak(0.1)
model.tdim(730)
model.R0(2.5)
model.infect_dur(14)
model.N_index(10)
model.setup()
for(i=1; i<=trials; i++) {
	model.run("quietly")
	res[i,.] = model.peak(), model.saved(), model.t_done()
}
mean(res)
end

// degree

mata:
model = covid()
model.N_nodes(10000)
model.degree(6)
model.pr_weak(0.1)
model.tdim(730)
model.R0(2)
model.infect_dur(14)
model.N_index(10)
model.setup()
for(i=1; i<=trials; i++) {
	model.run("quietly")
	res[i,.] = model.peak(), model.saved(), model.t_done()
}
mean(res)
end