function solve_demo
	addpath ..\data % poor man's import management
    addpath ..\src
    addpath ..\EZClimate-MATLAB
	addpath(genpath('..\dep\admat-2.0'))
	close all

	admat_startup
	carbon_price_by_jump(10)
	
	fprintf('solve demo finished\n')
    rmpath ..\data
    rmpath ..\src
    rmpath ..\EZClimate-MATLAB
    rmpath(genpath('..\dep\admat-2.0'))
end

function carbon_price_by_jump(ndecisions)
	%dtimes = [0, 15, 45, 85, 185, 285, 385];
	dtimes = [0:15:(ndecisions-2)*15, 285, 385];
	p_ups = [.5, .75];
	p_jumps = [.5, .6, .7, .8, .9];
	fname = sprintf('ndecisions_%d.md', ndecisions);
	diary off
	diary(fname)
	fprintf('solving %d-decision carbon price for %d jump and %d up probabilities\n', ndecisions, numel(p_jumps), numel(p_ups));
	mits = zeros(numel(p_ups) * numel(p_jumps), ndecisions);
	prices = zeros(size(mits));
	names = cell(size(prices, 1), 1);
	names(:) = {''};
	k = 0;
	for p_up = p_ups
	for p_jump = p_jumps
		k = k + 1;
		if p_up ~= p_ups(1) && p_jump == 0
			continue
		end
		[u, x, iter, param] = solve_carbon_price(dtimes, p_jump, p_up);
		[mits(k, :), prices(k, :)] = mean_path(x, param);
		fprintf('price_0 = %f, utility = %f, iter = %d, p_jump = %.2f, p_up = %.2f\n', ...
			prices(k, 1), u, iter, p_jump, p_up);
		names{k} = sprintf('p_{jump}=%.2f, p_{up}=%.2f', p_jump, p_up);
	end
	end
	diary off
	figure; plot(2015 + dtimes(1:end-1), mits'); legend(names);
	figure; plot(2015 + dtimes(1:end-1), prices'); legend(names);
	save(fname, 'mits', 'prices', 'names');
end

function [u, x, iter, param] = solve_carbon_price(dtimes, jump_prob, upjump_prob)
	mit = ones(1, 2^(numel(dtimes)-1) - 1);
	fcn = ADfun(@utility_v,1);
	opt = setgradopt('revprod', numel(mit));
	subinterval = 5;
	param = utility_v_config(jump_prob, dtimes, subinterval, upjump_prob);
	[u, x, ~, iter] = quasi_newton(fcn, mit, param, opt, false);
    utility_v(x, param);
end	

function [ms, ps] = mean_path(mits, param)
	prices = tax(mits, param);
	ms = zeros(1, param.ndecisions);
	ps = zeros(size(ms));
	for n = 1:param.ndecisions
		lb = 2 ^ (n-1);
        ub = 2 ^ n - 1;
		ms(n) = mits(lb:ub) * param.node_probs{n};
		ps(n) = prices(lb:ub) * param.node_probs{n};
	end
end
