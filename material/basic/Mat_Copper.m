classdef Mat_Copper < Mat_Magnetic & Mat_Thermal

    methods (Access = public)
        function obj = Mat_Copper()
            obj.Mur = 1;
            obj.Sigma = 5.96e7;
            obj.Rho = 8960;
            obj.Cp = 382;
            obj.Kappa = 401;
            obj.PlotColor = "#E9503E";
            obj.Cost = 10;
        end
    end
end