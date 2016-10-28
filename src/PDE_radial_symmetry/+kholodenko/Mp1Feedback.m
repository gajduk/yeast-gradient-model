classdef Mp1Feedback < kholodenko.Cascade
    %same as cascade from box 4 of Kholodenkos paper Signaling dynamics in space and time, with an added negative feedback

    %Membrane          | |---
    %                  V    |
    %Cytoplasm        Mp1   |
    %                  |    |
    %                  V    |
    %                 Mp2   |
    %                  |    |
    %                  V    |
    %                 Mp3 ---

    
    properties
        A = 1;
        kd = .5;
    end
    
    methods
        function [v_mem_kin,v1_phos,v1_kin,v2_phos,v2_kin,v3_phos] = getRates(self,u)
            mp3 = u(3);
            
            [v_mem_kin,v1_phos,v1_kin,v2_phos,v2_kin,v3_phos] = getRates@kholodenko.Cascade(self,u);
            
            v_mem_kin = v_mem_kin .* ( 1 + mp3 ./ self.kd ) ./ ( 1+ mp3 ./ self.kd .* self.A);
       end
    end
    
end

