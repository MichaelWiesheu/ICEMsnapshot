classdef Mat_Thermal < handle
    properties (Access = protected)
        Rho
        Cp
        Kappa
    end

    methods (Access = public)
        function rho = getRho(obj)
            rho = obj.Rho;
        end
        function rho = setRho(obj, rho)
            obj.Rho = rho;
        end
        function Cp = getCp(obj)
            Cp = obj.Cp;
        end
        function cp = setCp(obj, cp)
            obj.Cp = cp;
        end
        function kappa = getKappa(obj)
            kappa = obj.Kappa;
        end
        function kappa = setKappa(obj, kappa)
            obj.Kappa = kappa;
        end
    end
end