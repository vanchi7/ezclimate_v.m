function param = utility_v_config (jump_prob, dtimes, subinterval, upjump_prob, filename)
	if nargin < 1
        jump_prob = 0.05;
    end
    if nargin < 2
        dtimes = [0, 15, 45, 85, 185, 285, 385];
	end
	if nargin < 3
		subinterval = 5;
	end
	if nargin < 4
		upjump_prob = .5;
	end
	if nargin < 5
		filename = ['simulated_damages_' num2str(numel(dtimes) - 1) '_sub_' num2str(subinterval) '.mat'];
	end
    % dtimes     - Decision times in years. dtimes(end) is the effective limit of the
    %                last decision. No decision is made at dtimes(end).
    % ndecisions - Number of decision times.
    % probs      - Probabilities of nodes reaching each ending nodes.
    % bau_emit-levels - Business-as-usual scenario emission levels during each period.
    % bau_emits       - BAU emissions at the end of each period.
    % bau_cumemits    - BAU cumulative emissions at the end of each period.
    param.dtimes      = dtimes;
    param.elapses     = diff(param.dtimes);
    param.ndecisions  = numel(dtimes) - 1;
    param.jump_prob  = jump_prob;
	param.upjump_prob = upjump_prob;
	[param.probs, param.node_probs, param.adj_probs] = probabilities(param);
    param.bau_emit_levels = ...
        interp1([0, 30, 60], [52.0, 70.0, 81.4], param.dtimes, 'linear', 81.4);
    param.bau_emits      = (param.bau_emit_levels(1:end-1)+param.bau_emit_levels(2:end)) ./ 2 .* param.elapses;
    param.bau_co2_levels = param.bau_emit_levels .* 0.71 ./ 3.67 ./ 2.13;
    param.bau_cum_emits  = cumsum(param.bau_emits);
    param.subinterval    = subinterval;  % in years
    % Forcing and damage parameters
    param.ghg_start      = 400;
    param.ghg_end        = 1000;
    param.emit_pcts      = 1 - ([450, 650, 1000] - param.ghg_start) ./ (param.ghg_end - param.ghg_start);
    param.sink_start     = 35.596;
    param.forcing_start  = 4.926;
    param.forcing_p1     = 0.13173;
    param.forcing_p2     = 0.607773;
    param.forcing_p3     = 315.3785;
    param.forcing_log_p1 = 5.35067129;
    param.forcing_log_p2 = log(278.06340701);
    param.absorbtion_p1  = 0.94835;
    param.absorbtion_p2  = 0.741547;
    param.lsc_p1         = 285.6268;
    param.lsc_p2         = 0.88414;
    param.steps          = floor(param.elapses/param.subinterval);
    param.cum_forcings   = cum_forcing_table(param);
    dsim = load(filename);
	param.damage_coefs   = permute(dsim.damage_coefs, [2 1 3 4]);
    param.damage_rcomb   = permute(dsim.damage_rcomb, [3 2 1]);
    % Cost parameters
    x60 = 0.543; % McKinsey estimate of % mitigation for 60 euros
    x100 = 0.671;% McKinsey estimate of % mitigation for 100 euros
    k = log(x60./x100)./log(0.6);
    c = x60./(60*1.206).^k;  % 1.206 is euro to USD exchange rate
    param.g = c.*(1- 1./(k+1)).*(1./c).^((k+1)./k);
    param.a = ((k+1)./k);
    param.max_price      = 2500;
    param.tech_const     = 1.5;
    param.tech_scale     = 0;
    join_price           = 2000;
    param.cbs_level      = (join_price ./ param.g ./ param.a) .^ (1 ./ (param.a-1));
    param.cbs_deriv      = param.cbs_level / (join_price .* (param.a - 1));
    param.cbs_b          = param.cbs_deriv .* (param.max_price - join_price) ./ param.cbs_level;
    param.cbs_k          = param.cbs_level .* (param.max_price - join_price) .^ param.cbs_b;
    param.cons_per_ton   = 30460 / param.bau_emit_levels(1);
    % Utility parameters
    param.growth         = 1.015;
    param.rho            = 1 - 1/0.9;
    param.alpha          = 1 - 7.0;
    param.beta           = (1 - 0.005) .^ param.subinterval;
    param.potential_cons = repmat(param.growth, 1, numel(param.dtimes)) .^ param.dtimes;
end

% Indices of binomial tree. The last non-leaf layer has single child.
% %
% p    - Layer number of the tree. p=1 is head of the tree.
% curr - Indices of nodes at layer p.
% up   - Indices of up children.
% down - Indices of down children.
% %
% curr, up, down have the same number of elements, except for last two layers.
function [curr, up, down] = indices(p, param)
    if p <= param.ndecisions
        curr = 2^(p-1) : 2^p-1;
    else
        curr = (0:2^(param.ndecisions-1)-1) + 2^param.ndecisions;
    end
    if nargout > 1
        if p < param.ndecisions
            up   = 2 * curr;
            down = up + 1;
        else
            up   = curr + 2.^(p-1);
            down = [];
        end
    end
end

% Forcing and GHG level at the end of each period.
function [cum_forcing, cum_ghg] = forcing_and_ghg(mits, param)
    cum_sink    = zeros(size(mits)); cum_sink(1)    = param.sink_start;
    cum_forcing = zeros(size(mits)); cum_forcing(1) = param.forcing_start;
    cum_ghg     = zeros(size(mits)); cum_ghg(1)     = param.ghg_start;
    for p = 1:param.ndecisions
        [curr, up, down] = indices(p, param);
        nsteps = param.steps(p);
        bau_co2 = interp1([1, nsteps+1], param.bau_co2_levels(p:p+1) .* param.subinterval, 1:nsteps);
        for n = 1:nsteps
            grt = curr(cum_ghg(curr) > 260); lst = curr(cum_ghg(curr) <= 260);
            cum_forcing(grt) = cum_forcing(grt) ... 
                + param.forcing_log_p1 .* (log(cum_ghg(grt))-param.forcing_log_p2);
            cum_forcing(lst) = cum_forcing(lst) ...
                + param.forcing_log_p1 .* (log(260)-param.forcing_log_p2) ...
                + param.forcing_log_p1 ./ 260 .* (cum_ghg(lst)-260);
            
            ppm = bau_co2(n) .* (1-mits(curr));
            ghg_lsc = cum_ghg(curr) - (param.lsc_p1 + param.lsc_p2 .* cum_sink(curr));
            absorb = 0.5 .* param.absorbtion_p1 .* sign(ghg_lsc) .* abs(ghg_lsc) .^ param.absorbtion_p2;
            cum_sink(curr) = cum_sink(curr) + absorb;
            cum_ghg(curr)  = cum_ghg(curr) + ppm - absorb;
        end
        if p < param.ndecisions
            cum_sink(up)    = cum_sink(curr);    cum_sink(down)    = cum_sink(curr);
            cum_forcing(up) = cum_forcing(curr); cum_forcing(down) = cum_forcing(curr);
            cum_ghg(up)     = cum_ghg(curr);     cum_ghg(down)     = cum_ghg(curr);
        end
    end
end

function cum_forcings = cum_forcing_table(param)
    cum_forcings = zeros(param.ndecisions, numel(param.emit_pcts));
    for k = 1 : numel(param.emit_pcts)
        mits = repmat(param.emit_pcts(k), 1, 2.^param.ndecisions-1);
        cum_forcing = forcing_and_ghg(mits, param);
        cum_forcings(:, k) = cum_forcing(2.^(0:param.ndecisions-1));
    end
end

function [ps, qs, adj] = probabilities(param)
    dnodes = 2 ^ param.ndecisions - 1;
    nodes  = dnodes + 2 ^ (param.ndecisions-1);
    adj    = sparse(nodes, nodes);
    for p = 1:param.ndecisions
        lb = 2 ^ (p-1);
        ub = 2 ^ p - 1;
        for crr = lb:ub
            if p == param.ndecisions
                chd = crr - lb + 1 + dnodes;
                jmp = [];
            else
                chd = 2*crr + [0  1];
                jmp = 2*crr + [-1 2];
                jmp = jmp (2*lb <= jmp & jmp <= 2*ub+1);
            end
            if param.jump_prob == 0 || isempty(jmp)
                adj(chd, crr) = 1 / numel(chd);
            else
				adj(chd, crr) = (1-param.jump_prob) / numel(chd);
				if numel(jmp) == 1
					adj(jmp, crr) = param.jump_prob;
				else
					adj(jmp, crr) = param.jump_prob .* [param.upjump_prob; 1-param.upjump_prob];
				end
            end
        end
    end
    endings = nodes - 2 ^ (param.ndecisions - 1) + 1 : nodes;
    adj(endings, endings) = speye(numel(endings));
    s = speye(nodes);
    s = s(:, 1 : nodes - 2 ^ (param.ndecisions - 1));
    s = (adj ^ param.ndecisions) * s;
    assert(0 == nnz(s(1:dnodes, :)));
    assert(all(abs(1-sum(s(dnodes+1:end, :))) < 1e-14));
	ps = cell(1, param.ndecisions);
	qs = cell(1, param.ndecisions);
    for p = 1:param.ndecisions
        lb = 2 ^ (p-1);
        ub = 2 ^ p - 1;
		ps{p} = s(dnodes+1:end, lb:ub);
		q = adj^(p-1);
		qs{p} = full(q(lb:ub, 1));
    end
end
    