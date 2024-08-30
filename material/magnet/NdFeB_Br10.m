classdef NdFeB_Br10 < Mat_Magnetic & Mat_Thermal
    properties (Access = public)
        Br;
        Angle;
    end
    methods (Access = public)
        function obj = NdFeB_Br10()
            obj.Mur = 1.05;
            obj.Br = 1.0;
            obj.Sigma = 1/(1.4e-6);
            obj.Rho = 7500;
            obj.PlotColor = "#00715E";
            obj.Angle = [];
         end
    end
end