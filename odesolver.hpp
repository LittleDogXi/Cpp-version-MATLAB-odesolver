#pragma once
#include<limits>
#include<cmath>
#include<vector>
#include<string>
#include<iostream>
#include<iomanip>
/*
	Builded in UTF-8
	Done by Chenxi Ding
	2024/06/01
	Harbin, China

	请使用 UTF-8 重构！

	用于求解微分方程组的工具包，实现了MATLAB中的 ode solver

	a solver suite for ordinary differential equations, to achieve
	the major functions of MATLAB ode solver.

*/


namespace odesolver {

	// 实现MATLAB的realmin常数
	// 经过试验，这是与MATLAB函数等效的
#define realmin ((std::numeric_limits<double>::min)( ))

// 实现阉割版本的max和min函数
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)>(b))?(b):(a))

// 实现MATLAB的阉割版eps函数，缺省输入为1.0
	double eps(double x = 1.0)
	{
		 //double epsd = std::numeric_limits<double>::epsilon();
		 //return double(x) - std::nextafter(double(x), epsd);

		double xp = std::abs(x);
		double x1 = std::nextafter(xp, xp + 1.0);
		return x1 - xp;

	}
	// 实现MATLAB的sign函数，但是返回值为int
	int sign(double x)
	{
		if (x > 0.0)
		{
			return +1;
		}
		if (x < 0.0)
		{
			return -1;
		}
		if (x == 0.0)
		{
			return 0;
		}

	}

	typedef std::vector<double> state_type;
	typedef std::vector<bool> bool_vec;
	typedef std::vector<int> int_vec;
	// 事件函数类型
	typedef void(*event_fun) (double t, const state_type& y, state_type& value,
		bool_vec& isterminal, int_vec& direction, void* auxdata);

	// 打印向量到同一行
	template<typename T = state_type>
	void printVec(const T& v)
	{
		for (auto& value : v)
		{
			std::cout << value << '\t';
		}
		std::cout << '\n';
		// 如果这个向量为空就只会打印一个换行符
	}

	// 重载向量加法运算符
	template<typename T = double>
	std::vector<T> operator+(const std::vector<T>& lhs, const std::vector<T>& rhs)
	{
		int length = lhs.size();

		std::vector<T> result(length);

		for (int i = 0;i < length; i++)
		{
			result[i] = lhs[i] + rhs[i];
		}

		return result;
	}
	template<typename T = double>
	std::vector<T> operator-(const std::vector<T>& lhs, const std::vector<T>& rhs)
	{
		int length = lhs.size();

		std::vector<T> result(length);

		for (int i = 0;i < length; i++)
		{
			result[i] = lhs[i] - rhs[i];
		}

		return result;
	}
	// 重载向量数乘运算符
	template<typename T = double>
	std::vector<T> operator*(const std::vector<T>& vec, const T& scalar)
	{
		std::vector<T> result(vec);
		for (auto& value : result)
		{
			value = value * scalar;
		}
		return result;
	}
	// 左乘和右乘都定义一下
	template<typename T = double>
	std::vector<T> operator*(const T& scalar, const std::vector<T>& vec)
	{
		std::vector<T> result(vec);
		for (auto& value : result)
		{
			value = value * scalar;
		}
		return result;
	}
	template<typename T = double>
	std::vector<T> operator/(const std::vector<T>& vec, const T& scalar)
	{
		std::vector<T> result(vec);
		for (auto& value : result)
		{
			value = value / scalar;
		}
		return result;
	}

	// odeset类，指得是微分方程的一些设置
	struct OdeSet {
		double AbsTol;
		double RelTol;
		double InitialStep;
		double MaxStep;
		bool NormControl;
		event_fun event_f;

		// 默认构造函数
		OdeSet() : AbsTol(1.0e-6), RelTol(1.0e-3), InitialStep(0.0),
			MaxStep(0.0), NormControl(false), event_f(nullptr) { /*不做任何事*/
		};

		// 带参数的构造函数
		OdeSet(double atol = 1.0e-6, double rtol = 1.0e-3, double h1 = 0.0, double hmax = 0.0,
			bool isnormcontrol = false, event_fun f = nullptr)
		{
			this->AbsTol = atol;
			this->RelTol = rtol;
			this->InitialStep = h1;
			this->MaxStep = hmax;
			this->NormControl = isnormcontrol;
			this->event_f = f;
		}

	};

	// odesol类，指的是微分方程的输出
	struct OdeSol
	{
		state_type tout;
		std::vector<state_type> yout;
	};

	// ode_fun函数类型，可以被ode45函数调用，注意到这里给出了一个额外接口 void* auxdata
	// 出于简化的目的，这里的 auxdata 和事件函数中的额外接口共用
	typedef void(*ode_fun) (double t, const state_type& x, state_type& dxdt, void* auxdata);

	// 求范数， 1、2、无穷范数均可以求
	double norm(const state_type& vec, int mode = 2)
	{
		double output = 0.0;
		if (mode == 1)
		{
			for (auto& value : vec)
			{
				output += fabs(value);
			}
			return output;
		}
		if (mode == 2)
		{
			for (auto& value : vec)
			{
				output += value * value;
			}
			return sqrt(output);
		}
		if (mode == 3) // mode=3就是无穷范数
		{
			double temp;
			for (auto& value : vec)
			{
				temp = fabs(value);
				if (temp > output) { output = temp; }
			}
			return output;
		}
	}

	void ntrp45(double tinterp, double t, state_type& y, double h, state_type f[7], state_type& yinterp)
	{
		// BI 是 7×4 矩阵
		// BI = [
		// 1       -183/64      37/12       -145/128
		// 0          0           0            0
		// 0       1500/371    -1000/159    1000/371
		// 0       -125/32       125/12     -375/64 
		// 0       9477/3392   -729/106    25515/6784
		// 0        -11/7        11/3        -55/28
		// 0         3/2         -4            5/2
		// ];
		int dim = f[0].size();
		double BI[7][4] =
		{
			{1.0 * h,       -183.0 / 64.0 * h,      37.0 / 12.0 * h,       -145.0 / 128.0 * h},
			{0.0,          0.0,           0.0,            0.0},
			{0.0,       1500.0 / 371.0 * h,    -1000.0 / 159.0 * h,    1000.0 / 371.0 * h},
			{0.0,       -125.0 / 32.0 * h,       125.0 / 12.0 * h,     -375.0 / 64.0 * h},
			{0.0,       9477 / 3392.0 * h,   -729 / 106.0 * h,    25515 / 6784.0 * h},
			{0.0,        -11.0 / 7.0 * h,        11.0 / 3.0 * h,      -55.0 / 28.0 * h},
			{0.0,        3.0 / 2.0 * h,         -4.0 * h,            5.0 / 2.0 * h},
		};

		double s = (tinterp - t) / h;
		state_type cumprod1({ s, s * s, s * s * s, s * s * s * s });

		state_type temp[4];
		for (int jcol = 0;jcol < 4;jcol++)
		{
			temp[jcol].reserve(dim);
			for (int irow = 0;irow < dim;irow++)
			{
				temp[jcol][irow] = 0.0;
				for (int k = 0;k < 7;k++)
				{
					temp[jcol][irow] += f[k][irow] * BI[k][jcol];
				}
			}
		}
		yinterp = y + temp[0] * cumprod1[0]
			+ temp[1] * cumprod1[1] + temp[2] * cumprod1[2] + temp[3] * cumprod1[3];

	}

	// 找零点
	void odezero(event_fun f, state_type& v, double t, state_type& y, double tnew, state_type& ynew,
		double t0, double h, state_type f7[7], state_type& tout, std::vector<state_type>& yout,
		std::vector<bool_vec>& iout, state_type& vnew, bool stop, void* auxdata = NULL)
	{
		double tol = 128.0 * max(eps(t), eps(tnew));
		tol = min(tol, fabs(tnew - t));
		// state_type tout;
		// std::vector<state_type> yout;
		// std::vector<bool_vec> iout;
		int tdir = tnew - t > 0 ? +1 : -1;
		stop = false;
		double constexpr rmin = realmin;

		double tL = t;
		state_type yL = y;
		state_type vL = v;
		int eventdim = v.size(), odedim = y.size();
		// state_type vnew(eventdim);
		bool_vec isterminal(eventdim);
		int_vec direction(eventdim);
		f(tnew, ynew, vnew, isterminal, direction, auxdata);

		double tR = tnew, tswap;
		state_type yR = ynew, yswap;
		state_type vR = vnew, vswap;

		bool_vec useless1;
		int_vec useless2;
		double ttry = tR;
		double delta;
		state_type vtry(eventdim), ytry(odedim);
		while (true)
		{
			int lastmoved = 0;
			bool_vec indzc(eventdim);
			while (true)
			{
				bool isempty = false;
				for (int i = 0;i < eventdim;i++)
				{
					indzc[i] = (sign(vL[i]) != sign(vR[i])) && (direction[i] * (vR[i] - vL[i]) >= 0.0);
					isempty = isempty || indzc[i];
				}
				if (isempty)
				{
					if (lastmoved != 0)
					{
						// 这里要报一下错，后面再添加吧
					}
					return;
				}

				delta = tR - tL;
				if (fabs(delta) <= tol)
				{
					break;
				}
				double maybe;

				bool any = false;
				for (int i = 0;i < eventdim;i++)
				{
					if (indzc[i] == false) { continue; }
					else
					{
						if (vL[i] == 0.0 && vR[i] != 0.0)
						{
							any = true;
							break;
						}
					}

				}
				if (tL == t && any)
				{
					ttry = tL + tdir * 0.5 * tol;
				}
				else // Compute Regula Falsi change, using leftmost possibility.
				{
					double change = 1.0;
					for (int j = 0;j < eventdim;j++)
					{
						if (indzc[j] == false) { continue; }

						if (vL[j] == 0.0)
						{
							if ((tdir * ttry > tdir * tR) && (vtry[j] != vR[j]))
							{
								maybe = 1.0 - vR[j] * (ttry - tR) / ((vtry[j] - vR[j]) * delta);
								if ((maybe < 0.0) || (maybe > 1.0))
								{
									maybe = 0.5;
								}
							}
							else
							{
								maybe = 0.5;
							}
						}
						else if (vR[j] == 0.0)
						{
							if ((tdir * ttry < tdir * tL) && (vtry[j] != vL[j]))
							{
								maybe = vL[j] * (tL - ttry) / ((vtry[j] - vL[j]) * delta);
								if ((maybe < 0.0) || (maybe > 1.0))
								{
									maybe = 0.5;
								}
							}
							else
							{
								maybe = 0.5;
							}
						}
						else
						{
							maybe = -vL[j] / (vR[j] - vL[j]);
						}

						if (maybe < change)
						{
							change = maybe;
						}
					}

					change = change * fabs(delta);
					change = max(0.5 * tol, min(change, fabs(delta) - 0.5 * tol));
					ttry = tL + tdir * change;
				}

				ntrp45(ttry, t, y, h, f7, ytry);
				f(ttry, ytry, vtry, useless1, useless2, auxdata); // 只是为了语法不得不加一个

				isempty = false;
				for (int i = 0;i < eventdim;i++)
				{
					indzc[i] = (sign(vL[i]) != sign(vtry[i])) && (direction[i] * (vtry[i] - vL[i]) >= 0.0);
					isempty = isempty || indzc[i];
				}

				if (!isempty)
				{
					tswap = tR; tR = ttry; ttry = tswap;
					yswap = yR; yR = ytry; ytry = yswap;
					vswap = vR; vR = vtry; vtry = vswap;

					if (lastmoved == 2)
					{
						double temp;
						for (int i = 0;i < eventdim;i++)
						{
							temp = fabs(0.5 * vL[i]);
							if (temp > rmin)
							{
								vL[i] = temp;
							}
						}
					}
					lastmoved = 2;

				}
				else
				{
					tswap = tL; tL = ttry; ttry = tswap;
					yswap = yL; yL = ytry; ytry = yswap;
					vswap = vL; vL = vtry; vtry = vswap;
					if (lastmoved == 1)
					{
						double temp;
						for (int i = 0;i < eventdim;i++)
						{
							temp = fabs(0.5 * vR[i]);
							if (temp > rmin)
							{
								vR[i] = temp;
							}
						}
					}
					lastmoved = 1;
				}

			}

			tout.push_back(tR);
			yout.push_back(yR);
			iout.push_back(indzc);
			// isterminal
			bool any = false;
			for (int i = 0;i < eventdim;i++)
			{
				if (indzc[i] == false) { continue; }
				else
				{
					if (isterminal[i] == true)
					{
						any = true;
						break;
					}
				}
			}
			if (any == true)
			{
				if (tL != t0)
				{
					stop = true;
				}
				break;
			}
			else if (fabs(tnew - tR) <= tol)
			{
				break;
			}
			else
			{
				ttry = tR; ytry = yR; vtry = vR;
				tL = tR + tdir * 0.5 * tol;
				ntrp45(tL, t, y, h, f7, yL);
				f(tL, yL, vL, useless1, useless2, auxdata);
				tR = tnew; yR = ynew; vR = vnew;
			}
			// 明天检测一下
		}
	}

	// ode45函数
	template<int outputmode = 1> // 1代表只输出最后一个时刻的状态，其他表示输出中间时刻的状态
	void ode45(ode_fun f, double t0, double tf, const state_type& y0, OdeSol* sol, OdeSet* odeset,
		void* odeargs = nullptr)
	{
		// OdeSet赋值
		double atol = odeset->AbsTol, rtol = odeset->RelTol;
		double hmax = odeset->MaxStep, htry = odeset->InitialStep;
		bool isnormcontrol = odeset->NormControl;
		event_fun event_f = odeset->event_f;

		int refine = 1; // 这个功能还没添加，暂时设置为定值，即1
		// 一些别的参数
		double threshold = atol / rtol;
		int tdir;
		double htspan = tf - t0;
		if (htspan > 0) { tdir = +1; }
		else { htspan = -htspan; tdir = -1; } // 这里要保证htspan是大于0的
		int dim = y0.size();

		// 这里要矫正一下，因为默认hmax = 0，这是不能用的
		if (hmax == 0.0) { hmax = 0.1 * htspan; }
		double t = t0;
		state_type y(y0);
		state_type f0(dim);
		f(t0, y0, f0, odeargs); // 计算 f0 
		double normy = norm(y0);

		// 分配内存
		int nout = 0;
		int chunk;

		if (outputmode == 1)
		{
			sol->yout.reserve(1);
		}
		else
		{
			// std::cout << "要Output输出\n";
			chunk = min(max(100, 50 * refine), refine + floor(2048.0 / dim));
			std::cout << "chunk = " << chunk << std::endl;
			sol->tout.reserve(chunk);
			sol->yout.reserve(chunk);
			// std::cout << sol.tout.size() << '\n';
			// std::cout << sol.yout.size() << '\n';
			nout = 1;
			sol->tout.push_back(t);
			sol->yout.push_back(y);
		}

		/* --------定义一些常数------------*/
		// 用于
		double power = 1.0 / 5.0;
		// state_type A = ({1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0});
		double a2 = 1.0 / 5.0,
			a3 = 3.0 / 10.0,
			a4 = 4.0 / 5.0,
			a5 = 8.0 / 9.0,
			// 
			b11 = 1.0 / 5.0,
			b21 = 3.0 / 40.0,
			b31 = 44.0 / 45.0,
			b41 = 19372.0 / 6561.0,
			b51 = 9017.0 / 3168.0,
			b61 = 35.0 / 384.0,
			b22 = 9.0 / 40.0,
			b32 = -56.0 / 15.0,
			b42 = -25360.0 / 2187.0,
			b52 = -355.0 / 33.0,
			b33 = 32.0 / 9.0,
			b43 = 64448.0 / 6561.0,
			b53 = 46732.0 / 5247.0,
			b63 = 500.0 / 1113.0,
			b44 = -212.0 / 729.0,
			b54 = 49.0 / 176.0,
			b64 = 125.0 / 192.0,
			b55 = -5103.0 / 18656.0,
			b65 = -2187.0 / 6784.0,
			b66 = 11.0 / 84.0,
			// 一些常数
			e1 = 71.0 / 57600.0,
			e3 = -71.0 / 16695.0,
			e4 = 71.0 / 1920.0,
			e5 = -17253.0 / 339200.0,
			e6 = 22.0 / 525.0,
			e7 = -1.0 / 40.0;
		/* --------定义一些常数------------*/

		double absh;
		double hmin = 16.0 * eps(t); // 最小步长不能低于浮点精度的16倍
		// 计算初始步长
		if (htry == 0.0) // htry=0就是指没指定初始步长
		{
			absh = min(hmax, htspan);
			double rh;
			if (isnormcontrol) // 如果是范数控制
			{
				double normf0 = norm(f0);
				rh = (normf0 / max(normy, threshold)) / (0.8 * pow(rtol, power));
			}
			else // 如果不是范数控制
			{
				double norminf = 0.0, temp;
				for (int i = 0; i < dim;i++)
				{
					temp = fabs(f0[i] / max(fabs(y0[i]), threshold));
					if (temp > norminf)
					{
						norminf = temp;
					}
				}
				rh = norminf / (0.8 * pow(rtol, power));
			}
			if (absh * rh > 1.0)
			{
				absh = 1.0 / rh;
			}

			absh = max(absh, hmin);
		}
		else // 指定步长就修正一下即可
		{
			absh = min(hmax, max(hmin, htry));
		}

		state_type f1(f0), f2(dim), f3(dim), f4(dim), f5(dim), f6(dim), f7(dim), fE(dim);
		state_type y1(dim), y2(dim), y3(dim), y4(dim), y5(dim), y6(dim), ynew(dim); // 在这里定义一下，这样可以在循环外面赋值
		double t1, t2, t3, t4, t5, t6, tnew;
		double err, normynew;
		bool done = false;
		while (!done)
		{
			std::cout << std::setprecision(20) << "Now t = " << t << '\n';

			// By default, hmin is a small number such that t+hmin is only slightly
			// different than t.It might be 0 if t is 0.
			hmin = 16.0 * eps(t);
			absh = min(hmax, max(hmin, absh));
			double h = tdir * absh;
			if (1.1 * absh >= fabs(tf - t))
			{
				h = tf - t;
				absh = fabs(h);
				done = true;
			}

			bool nofailed = true;
			bool NNrejectStep, NNreset_f7;
			// state_type f1(dim), f2(dim), f3(dim), f4(dim), f5(dim), f6(dim), f7(dim), fE(dim);
			// state_type y1(dim), y2(dim), y3(dim), y4(dim), y5(dim), y6(dim);
			//double t1, t2, t3, t4, t5, t6, tnew;
			//double err, normynew;
			while (true)
			{
				y2 = y + h * (b11 * f1);
				t2 = t + h * a2;
				f(t2, y2, f2, odeargs);

				 std::cout << "t2 = " << t2 << std::endl;
				 std::cout << "y2 = ";
				 printVec(y2);

				y3 = y + h * (b21 * f1 + b22 * f2);
				t3 = t + h * a3;
				f(t3, y3, f3, odeargs);

				 std::cout << "t3 = " << t3 << std::endl;
				 std::cout << "y3 = ";
				 printVec(y3);

				y4 = y + h * (b31 * f1 + b32 * f2 + b33 * f3);
				t4 = t + h * a4;
				f(t4, y4, f4, odeargs);

				 std::cout << "t4 = " << t4 << std::endl;
				 std::cout << "y4 = ";
				 printVec(y4);

				y5 = y + h * (b41 * f1 + b42 * f2 + b43 * f3 + b44 * f4);
				t5 = t + h * a5;
				f(t5, y5, f5, odeargs);

				 std::cout << "t5 = " << t5 << std::endl;
				 std::cout << "y5 = ";
				 printVec(y5);

				y6 = y + h * (b51 * f1 + b52 * f2 + b53 * f3 + b54 * f4 + b55 * f5);
				t6 = t + h;
				f(t6, y6, f6, odeargs);

				 std::cout << "t6 = " << t6 << std::endl;
				 std::cout << "y6 = ";
				 printVec(y6);

				tnew = t + h;
				if (done)
				{
					tnew = tf;
				}
				h = tnew - t;

				ynew = y + h * (b61 * f1 + b63 * f3 + b64 * f4 + b65 * f5 + b66 * f6);
				f(tnew, ynew, f7, odeargs);

				 std::cout << "f7 = ";
				 printVec(f7);
				 std::cout << "ynew = ";
				 printVec(ynew);

				// nfevals += 6;

				// 估计误差
				NNrejectStep = false;
				fE = f1 * e1 + f3 * e3 + f4 * e4 + f5 * e5 + f6 * e6 + f7 * e7;
				if (isnormcontrol) // 如果是范数控制
				{
					normynew = norm(ynew);
					double errwt = max(max(normy, normynew), threshold);
					err = absh * (norm(fE) / errwt);
					// 这里没有使用nonNegative的功能
				}
				else // 如果是非范数控制
				{
					double norminf = 0.0, temp;
					for (int i = 0;i < dim;i++)
					{
						temp = fabs(fE[i]) / max(
							max( fabs(y[i]), fabs(ynew[i]) ), threshold
						);
						if (temp > norminf)
						{
							norminf = temp;
						}
					}
					err = absh * norminf;
					 std::cout << "err = " << err << '\n'; 
					// 这里没有使用nonNegative的功能
				}

				// Accept the solution only if the weighted error is no more than the
				// tolerance rtol.  Estimate an h that will yield an error of rtol on
				// the next step or the next try at taking this step, as the case may be,
				// and use 0.8 of this value to avoid failures.
				if (err > rtol)   // 失败的一步
				{
					// nfailed ++;
					if (absh <= hmin)
					{
						std::cout << "ode45在时刻 t = " << t
							<< " 的积分失败，为达到所需精度会使步长低于浮点相对精度！";
						return;
					}

					if (nofailed)
					{
						nofailed = false;
						if (NNrejectStep)
						{
							absh = max(hmin, 0.5 * absh);
						}
						else
						{
							absh = max(hmin, absh * max(0.1, 0.8 * pow((rtol / err), power)));
						}
					}
					else
					{
						absh = max(hmin, 0.5 * absh);
					}
					h = tdir * absh;
					done = false;
				}
				else // 成功的一步
				{
					NNreset_f7 = false;
					break;
				}

			}
			// nsteps++;

			// 触发事件
			if (event_f != nullptr)
			{
				std::vector<state_type> ff({ f1, f2, f3, f4, f5, f6, f7 });
				// 等后面补充 odezero(event_f, valt, t, y, tnew, ynew, t0, h, ff, );
				;
			}

			// 输出
			if (outputmode != 1)
			{
				nout++;
				int len = sol->tout.size(); // 一般来说nout == len
				int capacity = sol->tout.capacity();
				// std::cout << "size = " << len << '\n';
				// std::cout << "cap = " << capacity << '\n';
				if (len >= capacity)
				{
					// 当超过内存时，直接给再多赋值chunk
					sol->tout.reserve(capacity + chunk);
					sol->tout.reserve(capacity + chunk);
					// std::cout << "内存不够了" << sol.tout.capacity() <<'\n';
				}
				sol->tout.push_back(tnew);
				sol->yout.push_back(ynew);
			}


			if (done) { break; }

			if (nofailed)
			{
				double temp = 1.25 * pow((err / rtol), power);
				if (temp > 0.2)
				{
					absh = absh / temp;
				}
				else
				{
					absh = 5.0 * absh;
				}
			}

			t = tnew;
			y = ynew;
			if (isnormcontrol)
			{
				normy = normynew;
			}
			// if (NNreset_f7)
			// {
			//     f(tnew, ynew, f7, odeargs);
			//     // nfevals ++;
			// }

			f1 = f7;
		}

		if (outputmode == 1) // 如果返回的类型就是一个向量，直接赋值
		{
			sol->yout.push_back(ynew);
		}
		else // 如果返回的是一个那啥
		{
			sol->tout.shrink_to_fit();
			sol->yout.shrink_to_fit();
		}


	}

	// ode113函数
	void ode113(ode_fun odeFcn, double t0, double tf, const state_type& y0, state_type& sol, OdeSet* odeset,
		void* odeargs = nullptr)
	{
		// OdeSet赋值
		double atol = odeset->AbsTol, rtol = odeset->RelTol;
		double hmax = odeset->MaxStep, htry = odeset->InitialStep;
		bool isnormcontrol = odeset->NormControl;

		int refine = 1; // 这个功能还没添加，暂时设置为定值，即1
		// 一些别的参数
		double threshold = atol / rtol;
		int tdir;
		double htspan = tf - t0;
		if (htspan > 0.0) { tdir = +1; }
		else { htspan = -htspan; tdir = -1; } // 这里要保证htspan是大于0的
		int dim = y0.size();

		// 这里要矫正一下，因为默认hmax = 0，这是不能用的
		if (hmax == 0.0) { hmax = 0.1 * htspan; }
		double t = t0;
		state_type y(y0), yp;
		state_type f0(dim);
		odeFcn(t0, y0, f0, odeargs);
		yp = f0;
		double normy = norm(y0);

		int maxk = 12;
		state_type two({ 2., 4., 8., 16., 32., 64., 128., 256.,
			512., 1024., 2048., 4096., 8192. });
		state_type gstar({ 0.5000,  0.0833,  0.0417,  0.0264,
			  0.0188,  0.0143,  0.0114,  0.00936,
			  0.00789,  0.00679, 0.00592, 0.00524, 0.00468 });
		double hmin = 16.0 * eps(t), h;

		double absh;
		if (htry == 0.0)
		{
			absh = min(hmax, htspan);
			double rh;
			if (isnormcontrol)
			{
				rh = (norm(yp) / max(normy, threshold)) / (0.25 * sqrt(rtol));
			}
			else
			{
				double norminf = 0.0, temp;
				for (int i = 0;i < dim;i++)
				{
					temp = fabs(yp[i]) / max(fabs(y[i]), threshold);
					if (temp > norminf) { norminf = temp; }
				}
				rh = norminf / (0.25 * sqrt(rtol));
			}
			if (absh * rh > 1.0)
			{
				absh = 1.0 / rh;
			}
			absh = max(absh, hmin);
		}
		else
		{
			absh = min(hmax, max(hmin, htry));
		}

		int k = 1;
		int_vec K({ 1 });
		// std::vector<state_type> phi(14);for(int i = 0;i<14;i++){phi[i].reserve(dim);} // 分配内存
		std::vector<state_type> phi(14, state_type(dim, 0));
		phi[0] = yp;
		state_type psi(12, 0.0), alpha(12, 0.0), beta(12, 0.0), sig(13, 0.0);
		sig[0] = 1.0;
		state_type w(12), v(12), g(13);
		g[0] = 1.0;g[1] = 0.5;

		double hlast = 0.0;
		int klast = 0;
		bool phase1 = true;

		bool done = false;
		int failed, ns, knew, kold;
		double temp1, temp2, temp3, tlast;
		double err, erk, erkm1, erkm2, erkp1;
		double reduce;
		state_type invwt, phikp1(dim);
		/*这里预分配空间，否则会后续无法直接赋值*/
		// 注意要用resize，否则不能用 invwt[0] = ... 这样的方法赋值
		if (isnormcontrol) { invwt.resize(1); }
		else { invwt.resize(dim); }
		/*----------------------------------*/
		state_type p(dim), ylast(dim);
		while (!done)
		{
			hmin = 16.0 * eps(t);
			absh = min(hmax, max(hmin, absh));
			h = tdir * absh;

			if (1.1 * absh > fabs(tf - t))
			{
				h = tf - t;
				absh = fabs(h);
				done = true;
			}

			failed = 0;
			if (isnormcontrol)
			{
				invwt[0] = 1.0 / max(norm(y), threshold);
			}
			else
			{
				for (int i = 0;i < dim;i++)
				{
					invwt[i] = 1.0 / max(fabs(y[i]), threshold);
				}
			}

			while (true)
			{
				if (h != hlast) { ns = 0; }
				if (ns <= klast) { ns++; }
				if (k >= ns)
				{
					// 注意这里的 index 要比MATLAB的 index 小1
					beta[ns - 1] = 1.0;
					alpha[ns - 1] = 1.0 / double(ns);
					temp1 = h * double(ns);
					sig[ns] = 1.0;
					for (int i = ns + 1;i <= k;i++)
					{
						temp2 = psi[i - 2];
						psi[i - 2] = temp1;
						temp1 = temp2 + h;

						beta[i - 1] = beta[i - 2] * psi[i - 2] / temp2;
						alpha[i - 1] = h / temp1;
						sig[i] = double(i) * alpha[i - 1] * sig[i - 1];
					}
					psi[k - 1] = temp1;

					if (ns == 1)
					{
						int Klength = K.size();
						v.clear();
						w.clear();
						for (int i = 0;i < Klength;i++)
						{
							v.push_back(1.0 / double(K[i] * (K[i] + 1)));
							w.push_back(v[i]);
						}
					}
					else
					{
						// If order was raised, update diagonal part of v.
						if (k > klast)
						{
							// 注意，这里需要push_back，而不能用下标赋值
							// 原因在于 idx = k - 1 位置的可能越过了 size 的终点
							// 然而 k 通常只比 klast 大 1 ，因此直接push_back即可
							// v[k - 1] = 1.0 / double(k * (k + 1));
							v.push_back(1.0 / double(k * (k + 1)));
							for (int j = 1;j <= ns - 2;j++)
							{
								v[k - j - 1] = v[k - j - 1] - alpha[j] * v[k - j];
							}
						}
						// Update v and set w.
						for (int iq = 1;iq <= k + 1 - ns;iq++)
						{
							v[iq - 1] = v[iq - 1] - alpha[ns - 1] * v[iq];
							w[iq - 1] = v[iq - 1];
						}
						g[ns] = w[0];
					}

					for (int i = ns + 2;i <= k + 1;i++)
					{
						for (int iq = 1;iq <= k + 2 - i;iq++)
						{
							w[iq - 1] = w[iq - 1] - alpha[i - 2] * w[iq];
						}
						g[i - 1] = w[0];
					}

				}

				// Change phi to phi star.
				for (int i = ns + 1;i <= k;i++) // 这是列数，但是要减1
				{
					for (int j = 0;j < dim; j++) // 这是行数，不用减1
					{
						phi[i - 1][j] = phi[i - 1][j] * beta[i - 1];
					}
				}

				// Predict solution and differences.
				phi[k + 1] = phi[k];
				std::fill(phi[k].begin(), phi[k].end(), 0.0);
				// for(int i = 0;i < dim;i++)
				// {
				//     phi[k][i] = 0.0;
				// }

				std::fill(p.begin(), p.end(), 0.0);
				for (int i = k;i >= 1;i--)
				{
					p = p + g[i - 1] * phi[i - 1];
					phi[i - 1] = phi[i - 1] + phi[i];
				}

				p = y + h * p;
				tlast = t;
				t = tlast + h;
				if (done)
				{
					t = tf;
				}

				odeFcn(t, p, yp, odeargs);
				// nfevals ++;

				// Estimate errors at orders k, k-1, k-2.
				phikp1 = yp - phi[0];
				if (isnormcontrol)
				{
					temp3 = norm(phikp1) * invwt[0];
					err = absh * (g[k - 1] - g[k]) * temp3;
					erk = absh * sig[k] * gstar[k - 1] * temp3;
					if (k >= 2)
					{
						erkm1 = absh * sig[k - 1] * gstar[k - 2] *
							(norm(phi[k - 1] + phikp1) * invwt[0]);
					}
					else { erkm1 = 0.0; }
					if (k >= 3)
					{
						erkm2 = absh * sig[k - 2] * gstar[k - 3] *
							(norm(phi[k - 2] + phikp1) * invwt[0]);
					}
					else { erkm2 = 0.0; }

				}
				else
				{
					double norminf = 0.0, temptemp;
					for (int i = 0;i < dim;i++)
					{
						temptemp = fabs(phikp1[i] * invwt[i]);
						if (temptemp > norminf)
						{
							norminf = temptemp;
						}
					}
					temp3 = norminf;
					err = absh * (g[k - 1] - g[k]) * temp3;
					erk = absh * sig[k] * gstar[k - 1] * temp3;
					if (k >= 2)
					{
						double norminf = 0.0, temptemp;
						for (int i = 0;i < dim;i++)
						{
							temptemp = fabs((phi[k - 1][i] + phikp1[i]) * invwt[i]);
							if (temptemp > norminf)
							{
								norminf = temptemp;
							}
						}
						erkm1 = absh * sig[k - 1] * gstar[k - 2] * norminf;
					}
					else { erkm1 = 0.0; }
					if (k >= 3)
					{
						double norminf = 0.0, temptemp;
						for (int i = 0;i < dim;i++)
						{
							temptemp = fabs((phi[k - 2][i] + phikp1[i]) * invwt[i]);
							if (temptemp > norminf)
							{
								norminf = temptemp;
							}
						}
						erkm2 = absh * sig[k - 2] * gstar[k - 3] * norminf;
					}
					else { erkm2 = 0.0; }
				}

				// Test if order should be lowered
				knew = k;
				if ((k == 2) && (erkm1 <= 0.5 * erk))
				{
					knew = k - 1;
				}
				if ((k > 2) && (max(erkm1, erkm2) <= erk))
				{
					knew = k - 1;
				}

				// Test if step successful
				if (err > rtol) // Failed step
				{
					// nfailed ++;
					if (absh <= hmin)
					{
						std::cout << "ode113在时刻 t = " << t
							<< " 的积分失败，为达到所需精度会使步长低于浮点相对精度！";
						return;
					}

					// Restore t, phi, and psi.
					phase1 = false;
					t = tlast;

					int Klength = K.size();
					for (int i = 1;i <= Klength;i++)
					{
						phi[i - 1] = (phi[i - 1] - phi[i]) / beta[i - 1];
					}
					for (int i = 2;i <= k;i++)
					{
						psi[i - 2] = psi[i - 1] - h;
					}

					failed++;
					reduce = 0.5;
					if (failed == 3) { knew = 1; }
					else if (failed > 3) { reduce = min(0.5, sqrt(0.5 * rtol / erk)); }
					absh = max(reduce * absh, hmin);
					h = tdir * absh;
					k = knew;

					/*这一步对应ode113的458行，不知道对不对*/
					K.clear();
					for (int i = 1;i <= k;i++) { K.push_back(i); }
					/*----------------------------------*/

					done = false;
				}
				else // Successful step
				{
					break;
				}

			}
			// nsteps ++;

			klast = k;
			hlast = h;

			// Correct and evaluate.
			ylast = y;
			y = p + h * g[k] * phikp1;
			odeFcn(t, y, yp, odeargs);
			// nfevals ++;

			// Update differences for next step.
			phi[k] = yp - phi[0];
			phi[k + 1] = phi[k] - phi[k + 1];
			for (int i = 1;i <= K.size();i++)
			{
				phi[i - 1] = phi[i - 1] + phi[k];
			}

			if ((knew == k - 1) || (k == maxk))
			{
				phase1 = false;
			}

			// Select a new order.
			kold = k;
			if (phase1) // Always raise the order in phase1
			{
				k++;
			}
			else if (knew == k - 1) // Already decided to lower the order
			{
				k--;
				erk = erkm1;
			}
			else if (k + 1 <= ns) // Estimate error at higher order
			{
				if (isnormcontrol)
				{
					erkp1 = absh * gstar[k] * (norm(phi[k + 1]) * invwt[0]);
				}
				else
				{
					double norminf = 0.0, temptemp;
					for (int i = 0;i < dim;i++)
					{
						temptemp = fabs(phi[k + 1][i] * invwt[i]);
						if (temptemp > norminf)
						{
							norminf = temptemp;
						}
					}
					erkp1 = absh * gstar[k] * norminf;
				}
				if (k == 1)
				{
					if (erkp1 < 0.5 * erk)
					{
						k++;
						erk = erkp1;
					}
				}
				else
				{
					if (erkm1 <= min(erk, erkp1))
					{
						k--;
						erk = erkm1;
					}
					else if ((k < maxk) && (erkp1 < erk))
					{
						k++;
						erk = erkp1;
					}
				}

			}

			if (k != kold)
			{
				/*这一步对应ode113的517行，不知道对不对*/
				K.clear();
				for (int i = 1;i <= k;i++) { K.push_back(i); }
				/*----------------------------------*/
			}

			// NNreset_phi = false;

			if (done) { break; }

			if (phase1)
			{
				absh = 2.0 * absh;
			}
			else if (0.5 * rtol >= erk * two[k])
			{
				absh = 2.0 * absh;
			}
			else if (0.5 * rtol < erk)
			{
				reduce = pow((0.5 * rtol / erk), (1.0 / double(k + 1)));
				absh = absh * max(0.5, min(0.9, reduce));
			}
		}
		sol.resize(dim);
		sol = y;
	}

}

