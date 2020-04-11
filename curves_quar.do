clear all
run "c:\mijn documenten\projecten\stata\abm\covid\covid.mata"

set seed 1234567

local n    = 10000
local tdim = 600

tempfile tofill
set obs `tdim'
gen t = _n
save `tofill'

mata:
model = covid()
model.N_nodes(`n')
model.degree(10)
model.pr_weak(0.25)
model.R0(2.5)
model.infect_traj((5,8, 11, 14))
model.q_effect(1)
model.N_index(10)
model.tdim(`tdim')
end

foreach i in 0.25 0.05 {
	foreach j in 1 0.75 0.5 0.25 0 {
		do ..\covid\curves_quar_sub.do `i' `j' `n' `tofill'		
	}
}


reshape long S5_ S25_ I5_ I25_ R5_ R25_, i(t) j(eff)
rename S5_ S5
rename S25_ S25
rename I5_ I5
rename I25_ I25
rename R5_ R5
rename R25_ R25
reshape long S I R, i(t eff) j(weak)

label variable t "time" 
label variable eff "% effectiveness in reducing infections of quarantine"
label variable weak "% weak ties"
twby weak eff, compact  :   ///
    twoway line S I R t,    ///
	sort legend(cols(3))    ///
	lwidth(*2 ..)           ///
	scheme(s1rcolor)        ///
	ytitle(% of population) ///
	ylab(,angle(0))
exit
