% Vectorized utility function - trees are stored as cell arrays, with each layer stored in a cell.
function u = utility_v(mits, param)
    mitc = arrayfun(@(n) mits(2^(n-1) : 2^n - 1), 1:param.ndecisions, 'UniformOutput', false);
    cs   = costs(mits, mitc, param);
    cs   = arrayfun(@(n) cs(2^(n-1) : min(2^n - 1, numel(getval(cs)))), 1:param.ndecisions+1, 'UniformOutput', false);
    dt   = param.subinterval;
    ds   = damages(mitc, param);
    cons = max(param.potential_cons(end) .* (1 - ds{end}), 1e-18);
    u    = ((1 - param.beta) ./ (1 - param.beta .* (param.growth .^ param.rho))) .^ (1 / param.rho) .* cons;
    for p = param.ndecisions : -1 : 1
        pcons = max(param.potential_cons(p) .* (1 - ds{p}) .* (1 - cs{p}), 1e-18);
        if p < param.ndecisions
            cons  = max(cons .* (1 - repelem2(cs{p})) ./ (1 - cs{p+1}), 1.1e-18);
            pcons = repelem2(pcons);
        end
        for t = param.dtimes(p+1)-dt : -dt : param.dtimes(p)+dt
            seg  = t - param.dtimes(p);
            cons = cons ./ pcons;
            cons = sign(cons) .* abs(cons) .^ (seg / (seg + param.subinterval)) .* pcons;
            u    = ((1 - param.beta) .* cons .^ param.rho + param.beta .* u .^ param.rho) .^ (1 / param.rho);
        end
        if p < param.ndecisions
            pcons = pcons(1:2:end);
			u     = u .^ param.alpha;
			prob = param.adj_probs(2^p:2^(p+1)-1, 2^(p-1):2^p-1);
            u     = (u * prob) .^ (1/param.alpha);
        end
        cons = pcons;
        u = ((1 - param.beta) .* cons .^ param.rho + param.beta .* u .^ param.rho) .^ (1 / param.rho);
    end
    u = -u;
end

% Workaround of ADMAT bug - avoid duplicate indices.  
function y = repelem2(x)
    y = cons(zeros(1, 2 * numel(getval(x))), x);
    y(1:2:end) = x;
    y(2:2:end) = x;
end

% Average mitigation at the end of each period.
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

% Forcing and GHG level at the end of each period.
function [cum_forcing, cum_ghg] = forcing_and_ghg(mits, param)
    cum_sink    = cell(1, numel(mits)); cum_sink{1}    = cons(param.sink_start, mits{1});
    cum_forcing = cell(1, numel(mits)); cum_forcing{1} = cons(param.forcing_start, mits{1});
    cum_ghg     = cell(1, numel(mits)); cum_ghg{1}     = cons(param.ghg_start, mits{1});
    for p = 1:param.ndecisions
        nsteps = param.steps(p);
        bau_co2 = interp1([1, nsteps+1], param.bau_co2_levels(p:p+1) .* param.subinterval, 1:nsteps);
        for n = 1:nsteps
            grt = cum_ghg{p} > 260;
            idx = 1:numel(getval(cum_forcing{p}));
            g1  = idx(grt); g2 = idx(~grt);
            cum_forcing{p}(g1) = cum_forcing{p}(g1) ... 
                + param.forcing_log_p1 .* (log(cum_ghg{p}(g1))-param.forcing_log_p2);
            cum_forcing{p}(g2) = cum_forcing{p}(g2) ...
                + param.forcing_log_p1 .* (log(260)-param.forcing_log_p2) ...
                + param.forcing_log_p1 ./ 260 .* (cum_ghg{p}(g2)-260);
            
            ppm = bau_co2(n) .* (1-mits{p});
            ghg_lsc = cum_ghg{p} - (param.lsc_p1 + param.lsc_p2 .* cum_sink{p});
            absorb = 0.5 .* param.absorbtion_p1 .* sign(ghg_lsc) .* abs(ghg_lsc) .^ param.absorbtion_p2;
            cum_sink{p} = cum_sink{p} + absorb;
            cum_ghg{p}  = cum_ghg{p} + ppm - absorb;
        end
        if p < param.ndecisions
            cum_sink{p+1}    = repelem2(cum_sink{p});
            cum_forcing{p+1} = repelem2(cum_forcing{p});
            cum_ghg{p+1}     = repelem2(cum_ghg{p});
        end
    end
end

% Forcing based mitigation at the end of each period.
function fmits = forcing_based_mitigation(cum_forcing, param)
    fmits = cell(1, numel(cum_forcing));
    for p = 1 : param.ndecisions
        weight2 = cons(zeros(size(cum_forcing{p})), cum_forcing{1});
        weight3 = cons(zeros(size(cum_forcing{p})), cum_forcing{1});
        c1  = cum_forcing{p} > param.cum_forcings(p,1);
        c2  = cum_forcing{p} > param.cum_forcings(p,2);
        idx = 1:numel(getval(weight2));
        g1  = idx(c2); g2 = idx(c1 & ~c2); g3 = idx(~c1);
        weight2(g1) = (param.cum_forcings(p,3) - cum_forcing{p}(g1)) ...
            / (param.cum_forcings(p,3) - param.cum_forcings(p,2));
        weight2(g2) = (cum_forcing{p}(g2) - param.cum_forcings(p,1)) ...
            / (param.cum_forcings(p,2) - param.cum_forcings(p,1));
        weight3(g2) = (param.cum_forcings(p,2) - cum_forcing{p}(g2)) ...
            / (param.cum_forcings(p,2) - param.cum_forcings(p,1));
        weight3(g3) = 1 + (param.cum_forcings(p,1) - cum_forcing{p}(g3)) / param.cum_forcings(p,1);
        fmits{p} = weight2 .* param.emit_pcts(2) + weight3 .* param.emit_pcts(1);
    end
end

% Damages at the beginning of each period.
function ds = damages(mits, param)
    [cum_forcing, cum_ghg] = forcing_and_ghg(mits, param);
    fmits = forcing_based_mitigation(cum_forcing, param);
    ghg_extension = cellfun(@(cg) 1 ./ (1 + exp(0.05 .* (cg - 200))), cum_ghg, 'UniformOutput', false);
    ds = cell(1, numel(mits) + 1); ds{1} = 0;
    emit_pct1 = param.emit_pcts(1);
    emit_pct2 = param.emit_pcts(2);
    for p = 1:param.ndecisions
        if p == param.ndecisions
            prob = param.probs{p};
            fm   = fmits{p};
            ghg  = ghg_extension{p};
        else
            prob = param.probs{p+1};
            fm   = repelem2(fmits{p});
            ghg  = repelem2(ghg_extension{p});
        end
        dc11  = param.damage_coefs(p, :, 1, 1);
        dc12  = param.damage_coefs(p, :, 1, 2);
        dc13  = param.damage_coefs(p, :, 1, 3);
        dc22  = param.damage_coefs(p, :, 2, 2);
        dc23  = param.damage_coefs(p, :, 2, 3);
        v1    = dc22 * prob .* fm + dc23 * prob;
        v2    = dc11 * prob .* fm .^ 2 + dc12 * prob .* fm + dc13 * prob;
        idx1  = fm < emit_pct2;
        idx2  = fm < emit_pct1;
        v     = cons(zeros(1, size(prob, 2)), fm);
        dr1   = param.damage_rcomb(p, :, 1);
        idz   = dr1 > 1e-5;
        deriv = dc11 * 2 * emit_pct1 + dc12;
        decay = deriv ./ (dr1 * log(0.5));
        for k = 1:size(prob, 2)
            endings = prob(:, k)' ~= 0 & idz;
            if idx1(k)
                v(k) = v1(k);
            elseif idx2(k)
                v(k) = v2(k);
            elseif any(endings)
                dist = decay(endings) * (fm(k) - emit_pct1) + log(dr1(endings)) / log(0.5);
                 v(k) = (0.5 .^ dist) * prob(endings, k) ./ exp((fm(k) - emit_pct1) .^ 2 ./ 60);
            end
        end
        ds{p+1} = v + ghg;
    end
end

% Cost at the beginning of each period.
function cs = costs(mits, mitc, param)
    avg_mits  = average_mitigation(mitc, param);
    avg_mits  = [0, repelem2(avg_mits(1:floor(end/2))), avg_mits(ceil(end/2):end)];
    t         = repelem(param.dtimes, 2.^[0:param.ndecisions-1, param.ndecisions-1]);
    tech_term = (1 - ((param.tech_const + param.tech_scale .* avg_mits)./100)) .^ t;
    mits      = max(mits, 1e-18);
    mits      = [mits, mits(ceil(end/2):end)];
    cbs       = param.g .* (mits .^ param.a);
    backstop  = 1:numel(getval(mits));
    backstop  = backstop(mits > param.cbs_level);
    if ~isempty(backstop)
        m = mits(backstop);
        cbs(backstop) = (param.g .* param.cbs_level .^ param.a) + ((m - param.cbs_level) .* param.max_price ...
            - (param.cbs_b ./ (param.cbs_b-1)) .* m .* (param.cbs_k./m) .^ (1./ param.cbs_b) ...
            + param.cbs_b .* param.cbs_level .* (param.cbs_k./param.cbs_level) .^ (1./param.cbs_b) ./ (param.cbs_b-1));
    end
    cs = cbs .* tech_term ./ param.cons_per_ton;
end
