function cs = tax(mits, param)
	% This function converts mitigation levels at all nodes to carbon tax rate
	% at all nodes. It is similar to the costs() function except the paramters
	% g and a are slightly different. It is exactly like the 
	% price(obj, years, mitigation, ave_mitigation) in DLWCost but vectorized
	g = param.g*param.a;
	a = param.a - 1;
	mitc = arrayfun(@(n) mits(2^(n-1) : 2^n - 1), 1:param.ndecisions, 'UniformOutput', false);
	avg_mits  = average_mitigation(mitc, param);
	avg_mits  = [0, repelem(avg_mits(1:floor(end/2)), 2), avg_mits(ceil(end/2):end)];
	t         = repelem(param.dtimes, 2.^[0:param.ndecisions-1, param.ndecisions-1]);
	tech_term = (1 - ((param.tech_const + param.tech_scale .* avg_mits)./100)) .^ t;
	mits      = max(mits, 0);
	mits      = [mits, mits(ceil(end/2):end)];
	cbs       = g .* (mits .^ a);
	backstop  = 1:numel(getval(mits));
	backstop  = backstop(mits > param.cbs_level);
	if ~isempty(backstop)
		m = mits(backstop);
		cbs(backstop) = param.max_price - (param.cbs_k./m).^(1/param.cbs_b);
	end
	cs = cbs .* tech_term;
end
	
% From utility_v:
function avg_mits = average_mitigation(mits, param)
	cum_emits    = cell(1, numel(mits));
	cum_emits{1} = mits{1} .* param.bau_emits(1);
	avg_mits     = cell(1, numel(mits));
	avg_mits{1}  = mits{1};
	for p = 2:param.ndecisions
		cum_emits{p}  = repelem2(cum_emits{p-1}) + mits{p} .* param.bau_emits(p);
		avg_mits{p}   = cum_emits{p} ./ param.bau_cum_emits(p);
	end
	avg_mits = [avg_mits{:}];
end
function y = repelem2(x)
	y = cons(zeros(1, 2 * numel(getval(x))), x);
	y(1:2:end) = x;
	y(2:2:end) = x;
end
	