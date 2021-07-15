function regression_test
    addpath ..\data % poor man's import management
    addpath ..\src
    addpath(genpath('..\dep\admat-2.0'))
	
	jmp_prob = 0.05;
	dtimes = [0, 15, 45, 85, 185, 285, 385];
	subinterval = 5;
    param = utility_v_config(jmp_prob, dtimes, subinterval);
	mits = [0.697761696466805, 0.853034558526121, 0.725151474155350, 1.041134384715706, 0.956013916061470, 0.972747784169145, 0.761556319227664, 1.193478341960530, 1.131100930854668, 1.160530030246598, 1.064935584242798, 1.188648314948183, 1.102787739350373, 1.040760013236461, 0.758161616384580, 1.072694494043509, 1.043657439917421, 1.103646953294773, 1.067166449562689, 1.108561146082528, 1.073625421456471, 1.166068576642049, 1.107020167988350, 1.108644295692657, 1.074028405660717, 1.159305474516272, 1.100862708921946, 1.299000054907720, 1.251047498283711, 1.473525471256901, 0.962615261280115, 0.947676821810615, 0.928486127243796, 0.951197457006953, 0.920784107560767, 0.959397761955144, 0.926581508095382, 0.965076145659975, 0.938788619680513, 0.960236170192503, 0.926960554438880, 0.963920900310840, 0.937013152861551, 0.972787909944345, 0.946783569868324, 1.003402653777684, 0.973317942522225, 0.959875359868599, 0.926380295381803, 0.961954655493991, 0.934345616424460, 0.970336458363138, 0.942778155192331, 0.998186030430271, 0.968435922663007, 0.988238673384833, 0.959790822109244, 1.008117497565536, 0.983244208646592, 1.092365176637111, 1.049014605393291, 1.029001758398304, 0.994438317203424];
	u = -9.792223812669601;
    u_v  = utility_v(mits, param);

    should(max(abs(u - u_v)) < 1e-8, "utility value");

    mits_h = mits; mits_h(1) = mits_h(1) + 1e-8;
    u_vh   = utility_v(mits_h, param);
    d      = (u_vh - u_v) .* 1e8;

    admat_startup
    optf = setgradopt('forwprod', numel(mits));
    optr = setgradopt('revprod', numel(mits));
    
    [u_vf, g_vf]   = feval(ADfun(@utility_v, 1), mits, param, optf);
    [u_vr, g_vr]   = feval(ADfun(@utility_v, 1), mits, param, optr);

    should(abs(d - g_vf(1)) < 1e-5, "nd vs fd");
    should(max(abs(g_vf - g_vr)) < 1e-8, "fd vs rd");
    should(max(abs(u - u_vf)) < 1e-8, "fd utility value");
    should(max(abs(u - u_vr)) < 1e-8, "rd utility value");

    fprintf('regression test finished\n')

    rmpath ..\data
    rmpath ..\src
    rmpath(genpath('..\dep\admat-2.0'))
end

function should(b, msg)
    if ~b
        fprintf('failed test: %s\n', msg)
    end
end
