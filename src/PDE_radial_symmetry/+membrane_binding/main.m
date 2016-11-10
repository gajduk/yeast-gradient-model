function main()
	params = {'A','b_on','b_off','k_phos','ks_cat','s_phos'};
	for i=1:length(params)
		param = params{i};
		sens_anal(param)
	end
end

function sens_anal(param)
	pp = utils.param_pretty(param);
	kss = linspace(-10,10,7);
	sols = cell(1,length(kss));
	parfor k=1:length(kss)
		model = membrane_binding.CoreBindingModel();
		ks = 2.^kss(k);
		model.(param) = model.(param) * ks;
		plotter = model.run(100,100,50);
		sols{k} = plotter;
	end
	figure('Position',[100,100,800,400]);
	for k=1:length(kss)
		ks = 2.^kss(k);

		subplot(1,2,1)
		sols{k}.plot_gradient(1,sprintf('%s = %.4f',pp,ks),true);
		hold on

		subplot(1,2,2)
		sols{k}.plot_gradient(3,sprintf('%s = %.4f',pp,ks),true);
		hold on
	end
	legend show
	subplot(1,2,1)
	legend show
	suptitle(pp);
	export_fig(sprintf('%s sens_anal.png',param),'-r300')
end