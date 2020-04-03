cscript
cd "c:\mijn documenten\projecten\stata\abm\abm_nw"
run abm_nw.mata

mata :

mata set matastrict on

class covid
{
	real                  scalar R0
    real                  scalar infect_dur
	real                  scalar N
	real                  scalar degree
	real                  scalar pr_weak
	real                  scalar N_index
	real                  scalar tdim

	real                  scalar w
	class abm_nw          scalar network
	
    real                  matrix susceptible
	pointer (real vector) vector infectives
	real                  matrix dur
	real                  vector R
	real                  vector t_start
	real                  vector N_infectives
	real                  vector cumul_infected

	transmorphic                 R0()	
    transmorphic                 infect_dur()
	transmorphic                 N_nodes()
	transmorphic                 degree()
	transmorphic                 pr_weak()
	transmorphic                 N_index()
	transmorphic                 tdim()

	void                         is_posint()
	void                         is_pr()
	
	void                         setup() 
	void                         run()
	void                         step()
	

	void                         infect()
	void                         infect_neighbours()
	void                         copy_props()
	                        

	
	void                         toStata()
	real                  scalar peak()
	real                  scalar saved()
	real                  scalar t_done()
}


// -------------------------------- settings

transmorphic covid::R0(| real scalar val)
{
    if (args()==1) {
	    if (val < 0) _error("reproduction number must be positive")
		R0 = val
	}
	else {
	    return(R0)
	}
}

transmorphic covid::infect_dur(| real scalar val)
{
    if (args()==1) {
		is_posint(val)
		infect_dur = val
	}
	else {
	    return(infect_dur)
	}
}

transmorphic covid::N_nodes(| real scalar val)
{
    if (args()==1) {
	    network.N_nodes(0,val)
		N = val
	}
	else {
	    return(N)
	}
}

transmorphic covid::degree(| real scalar val) 
{
    if (args()==1) {
		is_posint(val)
	    degree = val
	}
	else {
	    return(degree)
	}
}

transmorphic covid::pr_weak(| real scalar val)
{
    if (args() == 1){
	    is_pr(val)
		pr_weak = val
	}
	else {
	    return(pr_weak)
	}
}

transmorphic covid::N_index(| real scalar val)
{
    if (args()==1) {
		is_posint(val)
		N_index = val
	}
	else {
	    return(N_index)
	}
}

transmorphic covid::tdim(| real scalar val)
{
    if (args()==1) {
	    is_posint(val)
		tdim = val
	}
	else {
	    return(tdim)
	}
}

// ----------------------------------------- checks
void covid::is_posint(real scalar val, | string scalar zero_ok) 
{
	string scalar errmsg
	
	if (args() == 1) {
		if (val <= 0 | mod(val,1)!= 0 ){
			errmsg = "scalar must be a positive integer"
			_error(3300, errmsg)
		}
	}
	else {
		if (val < 0 | mod(val,1)!= 0 ){
			errmsg = "scalar must be zero or a positive integer"
			_error(3300, errmsg)
		}
	}
}
void covid::is_pr(real scalar val)
{
	if (val < 0 | val > 1) {
		_error(3300,"scalar is not a valid probability")
	}
}

// ------------------------------------------ run the model

void covid::setup()
{	
    real scalar i
	
	network.clear()
    
	if(N == .){
	    _error("number of nodes needs to be set")
	}
	if(tdim == .) tdim = 365
    network.tdim(0)
	if (infect_dur == .) infect_dur = 10
	if (R0 == .) {
	    _error("reproduction number needs to be set")
	}
	if (N_index==.) N_index = 1
	if (degree == .) degree = 4
	if (pr_weak == .) pr_weak = .1
	
	network.directed(0)
	network.weighted(0)
	network.sw(degree,pr_weak)
	w = R0/(infect_dur*degree)
	if ( w > 1 ) w = 1
	
	network.setup("fast")
	N_infectives   = J(tdim,1,0)
	cumul_infected = J(tdim,1,0)
	R              = J(N,1,0)
	t_start        = J(N,1,.)
	dur            = J(N,tdim,0)
	susceptible    = J(N,tdim,1)
	infectives     = J(N,1,NULL)
	//Recovered/Removed = not in infectives & susceptible == 0
	
	for(i=1; i<=N; i++) {
	    infectives[i] = &J(1,0,.)
	}
	
	for(i=1;i<=N_index;i++){
		infect(1,i)
	}
}

void covid::infect(real scalar t, real scalar id)
{
    if (susceptible[id,t]) {
		infectives[t] = &(*infectives[t], id)
		dur[id,t] = 1
		susceptible[id,t] = 0
		N_infectives[t] = N_infectives[t] + 1
		cumul_infected[t] = cumul_infected[t] + 1
		t_start[id] = t
	}
}

void covid::infect_neighbours(real scalar t, real scalar id)
{
    real scalar i
	real vector cols
	
	cols = network.neighbours(id, 0, "dropped_ok", "fast")
	for (i=1; i <= length(cols); i++) {
	    if(runiform(1,1) < w) {
			R[id] = R[id] + susceptible[cols[i],t]
		    infect(t,cols[i])
		}
	}
}

void covid::step(real scalar t, string scalar quietly)
{
    real vector inf
	real scalar i, j, k
	
	// infect
	inf = *(infectives[t])
	k = cols(inf)
	j = 1
	for(i=1; i<= k; i++) {
	    dur[inf[i],t] = dur[inf[i],t] + 1
		if (dur[inf[i],t] <= infect_dur){
		    infect_neighbours(t,inf[i])
		}
		else {
			j = j + 1
		}
	}
	// end infect
	
	// heal 
	inf = *(infectives[t])
	k = cols(inf)	
	if     (j > 1 & j <=k) infectives[t] = &inf[|1,j \ 1,k|]
	else if(j==k+1)        infectives[t] = &J(1,0,.)
	N_infectives[t] = N_infectives[t] - j + 1
	// end heal
	
	copy_props(t)
	
	if (quietly == "") {
		"day " + strofreal(t) + ", " + strofreal(N_infectives[t]) + " current infectives, " + strofreal(cumul_infected[t]) + " cumulative"
	}
}


void covid::copy_props(real scalar t)
{
	susceptible[.,t+1] = susceptible[.,t]
	infectives[t+1] = &(*infectives[t])
	dur[.,t+1] = dur[.,t]
	N_infectives[t+1] = N_infectives[t]
	cumul_infected[t+1] = cumul_infected[t]
}

void covid::run(| string scalar quietly)
{
    real scalar i
	
	setup()
	
    for(i=1;i<tdim;i++) {
	    step(i, quietly)
	}
}

void covid::toStata(string scalar what)
{
    real matrix res, temp
	real scalar k , ks, xid, i
	string scalar varnames

	if(what == "inf") {
		res = (1..tdim)', N_infectives, cumul_infected
		varnames = "t", "curent", "cumul"
	}
	else if (what == "R") {
		res = (1..N)', t_start, R
		varnames = "id", "t_start", "R"
	}
	else if (what == "nw") {
		res = network.export_edgelist(0, "ego_all") 
        res = res, J(rows(res),tdim,.)
		temp = J(N, tdim, 0)
		k = cols(res)
        varnames = J(1,k,"")
        varnames[1..3] = "ego", "alter","weight"

		for(i=1; i<=tdim; i++) {
			temp[(*infectives[i])', i] = J(length(*infectives[i]),1,1)
            varnames[i+3] = "inf_t" + strofreal(i)
        }
		temp = temp :+ 2*(temp:==0 :& susceptible:==0)
		for (i=1; i<=rows(res);i++) {
			res[|i,4 \ i, .|] = temp[res[i,1],.]
        }

	}
	else _error(what + " is not a valid what")
	
	ks = st_nobs()
	k = rows(res)
	if (k>ks) st_addobs(k-ks)
	xid = st_addvar("long", varnames)
	st_store(.,xid,res)
}

real scalar covid::peak()
{
    return(max(N_infectives))
}

real scalar covid::saved()
{
    return(N-cumul_infected[tdim])
}

real scalar covid::t_done()
{
	real scalar i
    for(i=1; i <= tdim; i++) {
	    if (N_infectives[i] == 0) break 
	}
	return(i)
}

end	