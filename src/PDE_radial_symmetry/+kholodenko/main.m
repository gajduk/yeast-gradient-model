function main()
    all_sens_anal()
end

function all_sens_anal()
  params = {'A','v_mem_kin','D','V1_phos','k2_cat','V2_phos','k3_cat','V3_phos'};
  for i=1:length(params)
    if i > 1
      sens_anal(kholodenko.Cascade(),params{i});
    end
    sens_anal(kholodenko.Mp2Feedback(),params{i});
    sens_anal(kholodenko.Mp1Feedback(),params{i});
  end
end

function sens_anal(model,param,folder)
    base_val = model.(param);
    ks = -10:4:10;
    outputs = cell(1,length(ks));

    parfor i=1:length(ks)
       temp = feval(class(model));
       k = 1.5.^ -ks(i);
       temp.(param) = base_val * k;
       outputs{i} = temp.run(100,10,30);
    end

    figure('Position',[100,100,1500,600]);
    for i=1:length(ks)
      param_s = utils.param_pretty(param);
      k = 1.5.^ -ks(i);
      title_s = sprintf('%s = %3.3f',param_s,k);
      for w=1:3
        subplot(1,3,w)
        outputs{i}.plot_gradient(w,title_s);
        hold on
      end
    end
    
    for w=1:3
      subplot(1,3,w)
      legend show
    end
    
    export_fig(sprintf('figures kholodenko/%s/%s.png',class(model),param),'-transparent','-d300')
end