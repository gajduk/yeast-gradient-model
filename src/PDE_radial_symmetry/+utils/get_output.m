function WTf = get_output(model)
    spatial_steps = 43;
    time_steps = 100;

    end_time = 100;
    output1 = model.run(time_steps,end_time,spatial_steps);
    WTf = fliplr(output1.fraction_active_protein(end,:));
    WTf = WTf-WTf(end);
    WTf = WTf./WTf(2);
    WTf = WTf(:,2:4:43)';
end
