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

// ------------------------------- infect_traj
mata:
	model.infect_traj((5,8, 10, 14))
	assert(model.infect_traj() == (5,8, 10, 14))
end

rcof "mata: model.infect_traj((-5,8, 10, 14))" == 3300
rcof "mata: model.infect_traj((5.2,8, 10, 14))" == 3300
rcof "mata: model.infect_traj((9,8, 10, 14))" == 3498
rcof "mata: model.infect_traj((9,8, 10))" == 3498
// ------------------------------- traj_probs
mata:
	model.traj_probs((.6,.3, .5, .1))
	assert(model.traj_probs() == (.6,.3, .5, .1))
end

rcof "mata: model.traj_probs((1.6,.3, .5, .1))" == 3300
rcof "mata: model.traj_probs((.6,.3, .5))" == 3498


// ---------------------------------- degree
mata:
	model.degree(4)
	assert(model.degree() == 4)
end
rcof "mata: model.degree(2.5)" == 3300
rcof "mata: model.degree(-4)"  == 3300

mata:
	model.setup()
	assert(model.w == model.R0/(model.degree*model.infect_traj[4]))
end


// --------------------------------- is_posint
mata:
	for(i=1; i <=100; i++) {
		model.is_posint(i)
	}
	model.is_posint(0, "zero_ok")
	model.is_posint(0, "maarten is great") // any (silly) string will make zero ok
	model.is_posint((1..10))
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
	assert(model.infect_traj() == (5, 7, 9, 14))
	assert(model.traj_probs() == (0.65, 0.4, 0.4, 0.1))
	assert(model.degree()==4)
	assert(model.pr_weak() == .1)
	assert(model.network.directed()==0)
	assert(model.network.N_edges() == 4000) // this does not have to be the case 
	                                        // when an edge is rewired to an existing edge the number of edges will decrease
	assert(model.w == 1/14)
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
	assert(model.current[1,.] == (1, J(1,4,0),1))
	assert(model.cumul[1,.] == (1, J(1,4,0),1))
	true = J(1000,1,.)
	true[1] = 1
	assert(model.status==true)
	model.infect(1,100)
	assert(*model.infectives[1]==(1,100))
	assert(model.dur[100,1] == 1)
	assert(model.susceptible[100,1] == 0)
	assert(model.current[1,.] == (2, J(1,4,0),2))
	assert(model.cumul[1,.] == (2, J(1,4,0),2))
	true[100] = 1
	assert(model.status== true)
end

// ------------------------------ copy_props
mata:
	model.copy_props(1)
	assert(*model.infectives[2]==(1,100))
	assert(model.dur[100,2] == 1)
	assert(model.susceptible[100,2] == 0)
	assert(model.current[2,.] == (2, J(1,4,0),2))
	assert(model.cumul[2,.] == (2, J(1,4,0),2))
end

// ------------------------- disease_progression
mata:
	true = J(1000,1,.)
	true[1] = 1
	model = covid()
	model.N_nodes(1000)
	model.R0(4)
	model.setup()
	model.infect(1,1)
	model.disease_progression(1,1)
	assert(model.status == true)
	assert(model.current[1,.] == (1, J(1,4,0),1))
	assert(model.cumul[1,.] == (1, J(1,4,0),1))
	for(i=2; i <=1000; i++) {
		model.infect(1,i)
	}
	for(i=1; i<5; i++) {
		model.copy_props(i)
	}
	model.dur[.,5] = J(1000,1,5)
	for(i=1; i<=1000; i++){
		model.disease_progression(i,5)
	}
	assert(model.current[5,.] == (350,650,0,0,0,1000))
	assert(model.cumul[5,.] == (350,650,0,0,0,1000))
	assert(sum(model.status:==2)==650) // exactly right
	for(i=5; i<7; i++) {
		model.copy_props(i)
	}
	model.dur[.,7] = J(1000,1,7)
	for(i=1; i<=1000; i++){
		model.disease_progression(i,7)
	}
	assert(model.current[7,.] == (350,380,270,0,0,1000))
	assert(model.cumul[7,.] == (350,380,270,0,0,1000))
	assert(sum(model.status:==3)== 270) // close enough, on average should be 260
	for(i=7; i<9; i++) {
		model.copy_props(i)
	}
	model.dur[.,9] = J(1000,1,9)
	for(i=1; i<=1000; i++){
		model.disease_progression(i,9)
	}
	assert(model.current[9,.] == (350, 380, 175,95,0,1000))
	assert(model.cumul[9,.] == (350, 380, 175,95,0,1000))
	assert(sum(model.status:==4)==95) // close enough, on average should be 104
	for(i=9; i<=13; i++) {
		model.copy_props(i)
	}
	model.dur[.,14] = J(1000,1,14)
	for(i=1; i<=1000; i++){
		model.disease_progression(i,14)
	}	
	assert(model.current[14,.] == (350,380, 175, 83,12,1000))
	assert(model.cumul[14,.] == (350,380, 175, 83,12,1000))
	assert(sum(model.status:==5)==12) // close enough, on average should be 10
	
end

// ------------------------------- infect_neighbours
mata:
	model = covid()
	model.N_nodes(1000)
	model.R0(4)
	model.tdim(2)
	model.setup()
	for(i=1; i<1000; i = i + 5) {
		model.infect_neighbours(2,i)	
	}
	// 200 "origins" with 4 neighbours each with a chance of 1/14
	// on average should be 57 infections
	assert(model.current[2,6] == 58) // close enough
	index = (mod((1..1000),5):==1)'
	
	//each origin has 4 chances to infect 1/14 each, 
	// so should infect 4/14 ~= .2857 neighbours
	assert(mean(select(model.R,index))==.29) // close enough
end


// ------------------------------- step
mata:
	model = covid()
	model.N_nodes(1000)
	model.R0(4)
	model.setup()
	model.step(1,"")
	assert(*model.infectives[2]==(1))
	model.step(2, "")
	assert(*model.infectives[3]==(1))
	model.step(3, "quietly")
	assert(*model.infectives[4]==(1))
	// at t = 14 agent 1 should be healed
	for(i=4;i<=14;i++) {
		model.step(i, "")
	}
	assert((*model.infectives[13])[1] == 1)
	assert((*model.infectives[14])[1] != 1)
	assert(model.current[13,.] == (3,2,0,0,0,5))
	assert(model.current[14,.] == (3,1,0,0,0,4))
	assert(model.cumul[13,.] == (3,2,0,0,0,5))
	assert(model.cumul[14,.] == (3,2,0,0,0,5))
	
end

// -------------------------------- toStata
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
