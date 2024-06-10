# Cpp-version-MATLAB-odesolver
The repository achieve the major functions  of the MATLAB ode, like ode45 and ode113.

To using `ode45`, run like

```c++
ode45<mode>(odefun, t0, tf, y0, sol, options, pointer);
```

`odefun` is a function-pointer of the ode to be solved, callback as:

```
typedef std::vector<double> state_type;
void (*odefun)(double t, const state_type& x, state_type& dxdt, void* p);
```

`t0` is the starting time of the ode.

`tf` is the final time of the ode.

`y0` is the initial condition, of which the type is `std::vector<double>`.

`sol` is a pointer of the solution struct `OdeSol`, defined as below:

```c++
struct OdeSol
{
	state_type tout;
	std::vector<state_type> yout;
};
```

If `mode` is set `1`, only `yout` is used, to output the exact solution at the final time `tf`. If `mode` is set `2`, `tout` and `yout` are both used to store the whole solution.

`options` is used to set ode settings, defined as below:

```c++
struct OdeSet 
{
	double AbsTol;
	double RelTol;
	double InitialStep;
	double MaxStep;
	bool NormControl;
	event_fun event_f;

	// 默认构造函数
	OdeSet() : AbsTol(1.0e-6), RelTol(1.0e-3), InitialStep(0.0),
		MaxStep(0.0), NormControl(false), event_f(nullptr) { /*不做任何事*/};

	// 带参数的构造函数
	OdeSet(double atol = 1.0e-6, double rtol = 1.0e-3, double h1 = 0.0,
           double hmax = 0.0, bool isnormcontrol = false, event_fun f = nullptr)
	{
		this->AbsTol = atol;
		this->RelTol = rtol;
		this->InitialStep = h1;
		this->MaxStep = hmax;
		this->NormControl = isnormcontrol;
		this->event_f = f;
	}

};
```

`pointer` is used to pass additional parameters to the `odefun`.



To using `ode113`, run like

```
ode113(odefun, t0, tf, y0, sol, options, pointer);
```

`mode` is not used.

`sol` is always a `std::vector<double>` type variable.

The other arguments are the same as those in `ode45`.