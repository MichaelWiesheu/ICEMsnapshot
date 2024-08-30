classdef NdFeB_Br14 < Mat_Magnetic & Mat_Thermal & Mat_Mechanical
    properties %(Access = protected)
        Br;
        Angle;
    end
    methods (Access = public)
        function obj = NdFeB_Br14()
            obj.Mur = 1.05;
            obj.Br = 1.4;
            obj.Sigma = 1/(1.4e-6);
            obj.Rho = 7500;
            obj.PlotColor = "#00715E";
            obj.Cost = 50;
            obj.Nu = 0.3;
            obj.E = 180e9;
         end
    end
end