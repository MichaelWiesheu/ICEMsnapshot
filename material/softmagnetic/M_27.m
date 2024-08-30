classdef M_27 < Mat_Magnetic & Mat_Thermal & Mat_Mechanical
    % M27 = 330-50A
    methods (Access = public)
        function obj = M_27()
            obj.Mur = 2000;
            obj.PlotColor = [0.7, 0.7, 0.7];
            obj.Cost = 2;
            obj.Sigma = 1.03e7;
            obj.Rho = 7874;
            obj.Cp = 449;
            obj.type = "nonlinear";
            obj.B = [0.00 ,0.05 ,0.10 ,0.15 ,0.20 ,0.25 ,0.30 ,0.35 ,0.40 ,0.45 ,0.50 ,0.55 ,0.60 ,0.65 ,0.70 ,0.75 ,0.80 ,0.85 ,0.90 ,0.95 ,1.00 ,1.05 ,1.10 ,1.15 ,1.20 ,1.25 ,1.30 ,1.35 ,1.40 ,1.45 ,1.50 ,1.55 ,1.60 ,1.65 ,1.70 ,1.75 ,1.80 ,1.85 ,1.90 ,1.95 ,2.00 ,2.05 ,2.10 ,2.15 ,2.20 ,2.25 ,2.30];
            obj.H = [0.000000, 17.059828, 25.634971, 31.338354, 35.778997, 39.602124, 43.123143, 46.520439, 49.908177, 53.368284, 56.966318, 60.760271, 64.806105, 69.161803, 73.890922, 79.066315, 84.774676, 91.122722, 98.246299, 106.324591, 115.603417, 126.435124, 139.349759, 155.187082, 175.350538, 202.312017, 240.640455, 299.118027, 394.993386, 561.726177, 859.328763, 1375.466888, 2191.246914, 3328.145908, 4760.506172, 6535.339449, 8788.970657, 11670.804347, 15385.186211, 20246.553031, 26995.131141, 38724.496369, 64917.284463, 101489.309338, 137202.828961, 176835.706764, 216374.283609];
            obj.fitHBspline()

            obj.Nu = 0.3;
            obj.E = 180e9;
        end
    end
end