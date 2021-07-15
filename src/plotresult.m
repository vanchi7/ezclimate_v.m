function plotresult
	close all
	load('ndecisions_10.mat');
	ndecisions = 10;
	dtimes = 2015 + [0:15:(ndecisions-2)*15, 285];
	qs = [.5, .6, .7, .8, .9, .5, .6, .7, .8, .9,];
    names = arrayfun(@(q) ['$q =' num2str(q) '$'], qs, 'UniformOutput', false);
	names(1:5) = strcat(names(1:5), ' $(1/2$ up$)$');
	names(6:10) = strcat(names(6:10), ' $(3/4$ up$)$');
	clrs = rand(10, 3);
	syms = ['x';'+'; '.';'*';'^';'>';'<';'s';'+'; '.';'*';'^';'>';'<';'s'];
	plotpaths(dtimes, 100*mits, clrs, syms, names, 'Year', 'Percent of Mitigation')
	plotpaths(dtimes, prices, clrs, syms, names, 'Year', 'GHG Taxation Rate ($USD)')
end

function plotpaths(x, ys, clrs, syms, names, xname, yname)
	figure
	set(gcf, 'position', [200, 200, 500, 500])
	set(gca, 'FontName', 'Times New Roman');
	hold on
	for n = 1:size(ys, 1)
		plot(x, ys(n, :), [syms(n) '-'], 'color', clrs(n, :), 'LineWidth', 1.5, 'MarkerSize', 10);
	end
	hold off
	legend(names, 'Interpreter', 'latex')
	legend('boxoff')
	xlim([2000 2315])
	xlabel(xname)
	ylabel(yname)
end
