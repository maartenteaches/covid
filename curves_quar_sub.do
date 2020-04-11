args weak eff n tofill
drop _all
mata
model.q_effect(`eff')
model.pr_weak(`weak')
model.run("quietly")
model.toStata("inf")
end

local eff = round(`eff'*100)
local weak = round(`weak'*100)
gen S`weak'_`eff' = (`n' - c_total)/`n'*100
gen I`weak'_`eff' = total/`n'*100
gen R`weak'_`eff' = (c_total-total)/`n'*100

keep t S`weak'_`eff' I`weak'_`eff' R`weak'_`eff'
merge 1:1 t using `tofill'
drop _merge
save `tofill', replace