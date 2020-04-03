cscript

mata:

mata set matastrict on

class sir
{
	real vector  S
	real vector  I
	real vector  R
	
	real scalar  n
	real scalar  a
	real scalar  b
	real scalar  t
	
	transmorphic initial()
	transmorphic a()
	transmorphic b()
	transmorphic t()
	
	void         setup()
	void         step()
	void         run()
	
	void         toStata()
}

transmorphic sir::a(|real scalar val)
{
	if (args()==1) {
		if (val <= 0) _error("a needs to be positive")
		a = val
	}
	else {
		return(a)
	}
}
transmorphic sir::b(|real scalar val)
{
	if (args()==1) {
		if (val <= 0 ) _error("b needs to be a postive number")
		b = val
	}	
	else {
		return(b)
	}
}
transmorphic sir::t(|real scalar val)
{
	if (args()==1) {
		if (val <= 0 | mod(val,1) != 0){
			 _error("t needs to be a positive integer")
		}
		t = val
	}	
	else {
		return(b)
	}
}


transmorphic sir::initial(| real scalar S0, real scalar I0)
{
	if (args()==1) {
		_error("initial conditions need to be specified for both S and I")
	} 
	else if (args()==2) {
		if (S0 < 0 | I0 < 0) {
				_error("initial conditions need to be positive")
		}
		S = S0
		I = I0
	}
	else {
		return((S[1],I[1]))
	}
}

void sir::setup()
{
	if(S==J(1,0,.) | I==J(1,0,.)) _error("initial conditions need to be specified")
	if(a==. | b == .) _error("parameters a and b need to be specified")
	if(t==.) _error("time needs to be specified")
	
	S = S \ J(t-1,1,.)
	I = I \ J(t-1,1,.)
	R = 0 \ J(t-1,1,.)
	n = S[1]+I[1]
}

void sir::step(real scalar now)
{
	S[now] = S[now-1] - a*S[now-1]*I[now-1]/n
	I[now] = I[now-1] + a*S[now-1]*I[now-1]/n - b*I[now-1]
	R[now] = R[now-1] + b*I[now-1]
}

void sir::run()
{
	real scalar i
	
	setup()
	for(i=2; i<=t; i++) {
		step(i)
	}
}
end


mata:
model = sir()
model.initial(1000000,3)
model.a(2)
model.b(0.005)
model.t(20)
model.run()
model.S, model.I, model.R
end