cscript
set rmsg on
cd "c:\mijn documenten\projecten\stata\abm\covid"
set seed 123456

run covid.mata

// -------------------------------- N_nodes
mata:
	model = covid()
	model.N_nodes(5)
	model.N_nodes()
end

// ---------------------------------- tdim
mata:
	model = covid()
	model.N_nodes(200)
	model.tdim(100)
	model.tdim()
end

// --------------------------------- repro
mata:
	model.R0(4)
	assert(model.R0() == 4)
end
rcof "mata: model.R0(-5)" == 3498

// ------------------------------- infect_dur
mata:
	model.infect_dur(10)
	assert(model.infect_dur() == 10)
end
rcof "mata: model.infect_dur(-2)" == 3300
rcof "mata: model.infect_dur(2.4)" == 3300


// ---------------------------------- degree
mata:
	model.degree(4)
	assert(model.degree() == 4)
end
rcof "mata: model.degree(2.5)" == 3300
rcof "mata: model.degree(-4)"  == 3300

mata:
	model.setup()
	assert(model.w == model.R0/(model.degree*model.infect_dur))
end


// --------------------------------- is_posint
mata:
	for(i=1; i <=100; i++) {
		model.is_posint(i)
	}
	model.is_posint(0, "zero_ok")
	model.is_posint(0, "maarten is great") // any (silly) string will make zero ok
end
rcof "mata: model.is_posint(0)" == 3300
rcof "mata: model.is_posint(1.3)" == 3300
rcof "mata: model.is_posint(-1)" == 3300

// --------------------------------- pr_weak
mata:
	model.pr_weak(.5)
	assert(model.pr_weak() == .5)
end
rcof "mata: model.pr_weak(-2)" == 3300
rcof "mata: model.pr_weak(50)" == 3300



// --------------------------------- setup
mata:
	model = covid()
	model.N_nodes(1000)
	model.R0(4)
	model.setup()
	assert(model.tdim() == 365)
	assert(model.infect_dur() == 10)
	assert(model.degree()==4)
	assert(model.pr_weak() == .1)
	assert(model.network.directed()==0)
	assert(model.network.N_edges() == 4000) // this does not have to be the case 
	                                        // when an edge is rewired to an existing edge the number of edges will decrease
	assert(model.w == .1)
end

// ------------------------------ infect
mata:
	model = covid()
	model.N_nodes(1000)
	model.R0(4)
	model.setup()
	model.infect(1,1)
	assert(*model.infectives[1]==1)
	assert(model.dur[1,1] == 1)
	assert(model.susceptible[1,1] == 0)
	assert(model.N_infectives[1] == 1)
	assert(model.cumul_infected[1] == 1)
	model.infect(1,100)
	assert(*model.infectives[1]==(1,100))
	assert(model.dur[100,1] == 1)
	assert(model.susceptible[100,1] == 0)
	assert(model.N_infectives[1] == 2)
	assert(model.cumul_infected[1] == 2)
end

// ------------------------------ copy_props
mata:
	model.copy_props(1)
	assert(*model.infectives[2]==(1,100))
	assert(model.dur[100,2] == 1)
	assert(model.susceptible[100,2] == 0)
	assert(model.N_infectives[2] == 2)
	assert(model.cumul_infected[2] == 2)
end

// ------------------------------- infect_neighbours
mata:
	model.infect_neighbours(2,506)
	assert(*model.infectives[2] == (1,100))
	model.infect_neighbours(2,506)
	assert(*model.infectives[2] == (1, 100, 508))
end

// ------------------------------- step
mata:
	model = covid()
	model.N_nodes(1000)
	model.R0(4)
	model.setup()
	model.step(1,"")
	assert(*model.infectives[2]==(1,999))
	model.step(2, "")
	assert(*model.infectives[3]==(1,999))
	model.step(3, "quietly")
	assert(*model.infectives[4]==(1,999,1000))
	// at t = 10 agent 1 should be healed
	for(i=4;i<=11;i++) {
		model.step(i, "")
	}
	assert((*model.infectives[9])[1] == 1)
	assert((*model.infectives[10])[1] != 1)
end

// -------------------------------- R
mata:
	model = covid()
	model.N_nodes(1000)
	model.R0(4)
	model.setup()
	model.run()
	model.toStata("inf")
	model.toStata("R")
	model.toStata("nw")
end
