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
                  lcolor(white) scheme(s1rcolor) ||            ///
        scatter y_ego x_ego if first & inf_t`t'==1             ///
                   mcolor(red) msymbol(O) ||                   ///
        scatter y_ego x_ego if first & inf_t`t'==0,            ///
                   mcolor(yellow) msymbol(O)   ||              ///
        scatter y_ego x_ego if first & inf_t`t'==              ///
                   mcolor(green) msymbol(O)                    ///	
                   aspect(1) xscale(off) yscale(off)           ///
                   legend(order(3 "susceptible"                ///
				                2 "infective"                  ///
								4 "removed") rows(1))          ///
                   name(t`t', replace) title(Time `t')     
}
