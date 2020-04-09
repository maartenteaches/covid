cscript
cd "c:\mijn documenten\projecten\stata\abm\abm_nw"
run abm_nw.mata

local no_symptoms    = 1
local mild_symptoms  = 2
local needs_hospital = 3
local needs_icu      = 4
local dead           = 5
local infective      = 6
local immune         = 7

mata :

mata set matastrict on

class covid
{
    real                  scalar R0
    real                  vector infect_traj // onset symptoms, need hospital, need ICU, over
    real                  vector traj_probs   
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
    real                  vector status
    real                  matrix current
    real                  matrix cumul

    transmorphic                 R0()   
    transmorphic                 infect_traj()
    transmorphic                 traj_probs()
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
	void                         disease_progression()
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

transmorphic covid::infect_traj(| real vector val)
{
    if (args()==1) {
        is_posint(val)
        if(cols(val)!=4) _error("4 columns need to be specified")
        if(val[2]<val[1] | val[3] < val[2] | val[4]<val[3]) {
            _error("the values need to be ordered")
        }
        infect_traj = val
    }
    else {
        return(infect_traj)
    }
}

transmorphic covid::traj_probs(| real vector val)
{
    if (args()==1) {
        is_pr(val)
        if(cols(val)!=4) _error("4 columns need to be specified")
        traj_probs = val
    }
    else {
        return(traj_probs)
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
void covid::is_posint(real matrix val, | string scalar zero_ok) 
{
    string scalar errmsg
    
    if (args() == 1) {
        if (any(val :<= 0) | any(mod(val,1))){
            errmsg = "argument must be a positive integer"
            _error(3300, errmsg)
        }
    }
    else {
        if (any(val :< 0) | any(mod(val,1)) ){
            errmsg = "argument must be zero or a positive integer"
            _error(3300, errmsg)
        }
    }
}
void covid::is_pr(real matrix val)
{
    if (any(val :< 0) | any(val :> 1)) {
        _error(3300,"argument is not a valid probability")
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
    if (infect_traj == J(1,0,.)) infect_traj = (5, 7, 9, 14)
    if (traj_probs == J(1,0,.)) traj_probs = (0.65, 0.4, 0.4, 0.1)
    if (R0 == .) {
        _error("reproduction number needs to be set")
    }
    if (N_index==.) N_index = 1
    if (degree == .) degree = 4
    if (pr_weak == .) pr_weak = .1
    
    network.directed(0)
    network.weighted(0)
    network.sw(degree,pr_weak)
    w = R0/(infect_traj[4]*degree)
    if ( w > 1 ) w = 1
    
    network.setup("fast")
    current        = J(tdim,6,0)
    cumul          = J(tdim,6,0)
    R              = J(N,1,0)
    t_start        = J(N,1,.)
    status         = J(N,1,.)
    dur            = J(N,tdim,0)
    susceptible    = J(N,tdim,1)
    infectives     = J(tdim,1,NULL)
    //Recovered/Removed = not in infectives & susceptible == 0

    for(i=1; i<=tdim; i++) {
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
        current[t,`infective'] = current[t, `infective'] + 1
        cumul[t, `infective'] = cumul[t, `infective'] + 1
		current[t,`no_symptoms'] = current[t, `no_symptoms'] + 1
        cumul[t, `no_symptoms'] = cumul[t, `no_symptoms'] + 1	
        t_start[id] = t
        status[id] = `no_symptoms'
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

void covid::disease_progression(real scalar id, real scalar t)
{
	if (dur[id,t] == infect_traj[1] & 
        runiform(1,1)<traj_probs[1]) {
            status[id] = `mild_symptoms'
            current[t,`mild_symptoms'] = current[t, `mild_symptoms'] + 1
            cumul[t, `mild_symptoms'] = cumul[t, `mild_symptoms'] + 1
            current[t,`no_symptoms'] = current[t, `no_symptoms'] - 1
            cumul[t, `no_symptoms'] = cumul[t, `no_symptoms'] - 1
    }
    if (dur[id,t] == infect_traj[2] & 
        status[id]==`mild_symptoms' &
        runiform(1,1) < traj_probs[2]) {
            status[id] = `needs_hospital'
            current[t,`needs_hospital'] = current[t, `needs_hospital'] + 1
            cumul[t, `needs_hospital'] = cumul[t, `needs_hospital'] + 1
            current[t,`mild_symptoms'] = current[t, `mild_symptoms'] - 1
            cumul[t, `mild_symptoms'] = cumul[t, `mild_symptoms'] - 1
    }
    if (dur[id,t] == infect_traj[3] & 
        status[id]==`needs_hospital' &
        runiform(1,1) < traj_probs[3]) {
            status[id] = `needs_icu'
            current[t,`needs_icu'] = current[t, `needs_icu'] + 1
            cumul[t, `needs_icu'] = cumul[t, `needs_icu'] + 1
            current[t,`needs_hospital'] = current[t, `needs_hospital'] - 1
            cumul[t, `needs_hospital'] = cumul[t, `needs_hospital'] - 1
    }
    if (dur[id,t] == infect_traj[4] & 
        status[id]==`needs_icu' &
        runiform(1,1) < traj_probs[4]) {
            status[id] = `dead'
            current[t,`dead'] = current[t, `dead'] + 1
            cumul[t, `dead'] = cumul[t, `dead'] + 1
            current[t,`needs_icu'] = current[t, `needs_icu'] - 1
            cumul[t, `needs_icu'] = cumul[t, `needs_icu'] - 1
    }
	
}

void covid::step(real scalar t, string scalar quietly)
{
    real vector inf, nhealed
    real scalar i, k
    
    // infect
    inf = *(infectives[t])
    k = cols(inf)
    nhealed = J(1,5,0), 1
    for(i=1; i<= k; i++) {
        dur[inf[i],t] = dur[inf[i],t] + 1
        if (dur[inf[i],t] <= infect_traj[4]){
            infect_neighbours(t,inf[i])
			disease_progression(inf[i],t)
        }
        else {
            nhealed[6] = nhealed[6] + 1
            nhealed[status[inf[i]]] = nhealed[status[inf[i]]] + 1
        }
    }
    // end infect
    
    // heal 
    inf = *(infectives[t])
    k = cols(inf)   
    if     (nhealed[6] > 1 & nhealed[6] <=k) infectives[t] = &inf[|1,nhealed[6] \ 1,k|]
    else if(nhealed[6]==k+1)                 infectives[t] = &J(1,0,.)
    nhealed[6] = nhealed[6] -1
    current[t,.] = current[t,.] :- nhealed 
    // end heal
    
    copy_props(t)
    
    if (quietly == "") {
        "day " + strofreal(t) + ", " + strofreal(current[t,6]) + " current infectives, " + strofreal(cumul[t,6]) + " cumulative"
    }
}


void covid::copy_props(real scalar t)
{
    susceptible[.,t+1] = susceptible[.,t]
    infectives[t+1] = &(*infectives[t])
    dur[.,t+1] = dur[.,t]
    current[t+1,.] = current[t,.]
    cumul[t+1,.] = cumul[t,.]
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
        res = (1..tdim)', current, cumul
        varnames = "t", "no_symptoms", "mild_symptoms", "needs_hospital", "needs_icu", "dead", "total", 
        "c_no_symptoms", "c_mild_symptoms", "c_needs_hospital", "c_needs_icu", "c_dead", "c_total"
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


end 