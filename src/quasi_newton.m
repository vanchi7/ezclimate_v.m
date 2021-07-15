function [f , x , g, iter] = quasi_newton(fcn, x, param, opt, verbose, ftol, gtol)
	if nargin < 4
		error('missing argument for admat opt');
	end
	if nargin < 5
		verbose = true;
	end
	if nargin < 6
		ftol = 1e-5;
	end
	if nargin < 7
		gtol = 1e-7;
	end
	max_iter = 500;
	n = numel(x);
	H = eye(n);
	[f, g] = feval(fcn, x, param, opt);
	iter = 0;
	norm_g = norm(g);
	df = 1;
	while norm_g > gtol && abs(df) > ftol && iter < max_iter
		p = -g*H;
		if g'*p > 0
			H = eye(n);
			p = -g;
		end
		a = line_search(fcn, f, g, x, p, param,opt);
		s = a*p;
		x = x + s;
		[f_new, g_new] = feval(fcn, x, param, opt);
		df = f_new - f;
		dg = g_new - g;
		f = f_new;
		g = g_new;
		if iter == 0
			H = dg*s'/(dg*dg')*eye(n);
		end
		iter = iter + 1;
		rho = 1/(dg*s');
		H = (eye(n) - rho*s'*dg) * H * (eye(n) - rho*dg'*s) + rho*(s'*s);
		norm_g = norm(g);
		if(verbose)
			fprintf('iter #%d, |g| = %f, f = %f\n', iter, norm_g, f);
		end
	end
end

function a = line_search(fcn , f , g, x , p, param,opt)
	a = 1;
	rho = 0.5;
	c = 1e-4;
	fval = feval(fcn, x + a*p, param, opt);
	while (fval > f + c*a*p*g' && a > 1e-5)
		a = a*rho;
		fval = feval(fcn, x + a*p, param, opt);
	end
end
