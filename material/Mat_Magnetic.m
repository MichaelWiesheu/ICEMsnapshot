classdef Mat_Magnetic < handle
    properties (Access = public) %
        PlotColor=[1,1,1];
        Mu0=4*pi*1e-7
        type="linear"
        Cost
        Sigma
        Mur
        B
        H
        HBspline;
        HBsplineDer;
    end

    methods (Access = public)
        function mu = getMuLinear(obj)
            mu = obj.Mur*obj.Mu0;
        end
        function nu = getNuLinear(obj)
            nu = 1./obj.getMuLinear();
        end
        function setMurLinear(obj, mur)
            obj.Mur = mur;
        end
        function setSigma(obj, sigma)
            obj.Sigma = sigma;
        end
        function sig = getSigma(obj)
            sig = obj.Sigma;
        end
        function sig = getCost(obj)
            sig = obj.Cost;
        end
        function setType(obj, type)
            if any(type == ["linear", "nonlinear"])
                obj.type = type;
            else
                warning("Only linear or nonlinear as type allowed.")
            end
        end
        function type = getType(obj)
            type = obj.type;
        end
        function col = getPlotColor(obj)
            col = obj.PlotColor;
        end
        function fitHBspline(obj)
            obj.HBspline = pchip(obj.B, obj.H);
            obj.HBsplineDer = fnder(obj.HBspline);
        end

        function mu = getMuNonlinear(obj, B)
            B(B<=1e-5) = 1e-5; % avoid zero division
            mu = B./ppval(obj.HBspline, B);
            % correct B values higher than measurement (assume behavior like vacuum)
            mu(B>obj.B(end)) = B(B>obj.B(end))./((B(B>obj.B(end))-obj.B(end))./obj.Mu0 + obj.H(end));
        end
        function nu = getNuNonlinear(obj, B)
            nu = 1./obj.getMuNonlinear(B);
        end

        function nuPrime = getNuPrimeNonlinear(obj, B)
            B(B<=0) = 1e-5; % avoid zero division
            nuPrime = (ppval(obj.HBsplineDer, B).*B - ppval(obj.HBspline, B))./(B.^2);
            % correct B values higher than measurement (assume behavior like vacuum)
            B2correct = B(B>obj.B(end));
            nuPrime(B>obj.B(end)) = (B2correct./obj.Mu0 - obj.getNuNonlinear(B2correct).*B2correct)./B2correct.^2; 

            % old version
%             Step = 1e-6;
%             nu1 = obj.getNuNonlinear(B);
%             nu2 = obj.getNuNonlinear(B + Step);
%             nuPrime = (nu2 - nu1) / (Step);
        end

        function nuPrimeB2 = getNuPrimeB2Nonlinear(obj, B)
            % TBD NOT YET FINISHED
%             B2nuSpline = pchip(obj.B.^2, obj.getNuNonlinear(obj.B));
%             B2nuSplineDer = fnder(B2nuSpline);
%             nuPrimeB2 = ppval(B2nuSplineDer, B.^2);
            nuPrimeB2 = obj.getNuPrimeNonlinear(B)./(2*B);
        end

        function C = getCoenergy(obj, B)
            if obj.getType() == "linear"
                C = 1/(2*obj.getMuLinear())*B.^2;
            else
                C = B.^2./obj.getMuNonlinear(B) - obj.getEnergy(B);
            end
        end
        
        function E = getEnergy(obj, B)
            if obj.getType() == "linear"
                E = 1/(2*obj.getMuLinear())*B.^2;
            else
                E = zeros(size(B));
                cfun = @(b) ppval(obj.HBspline, b);
                for iB = 1:numel(B)
                    E(iB) = integral(cfun, 0, B(iB));
                end
            end
        end 

        function plotCharacteristics(obj)
            f = figure("Name", "Plots");
            % B-H curve
            subplot(2,2,1)
            obj.plotBH();
            % B-mu-curve
            subplot(2,2,2)
            obj.plotBMu();
            % B-nu-curve
            subplot(2,2,3)
            obj.plotBNu();
            % B-nuprime-curve
            subplot(2,2,4)
            obj.plotBNuPrime();
        end

        function plotBH(obj)
            b = 0:0.001:max(obj.B)*1.1;
            h = b./obj.getMuNonlinear(b);
            % B-H curve
            plot(h, b, "Color", "black", LineWidth=1.5)
            hold on
            scatter(obj.H, obj.B, "filled", "MarkerEdgeColor","black", "Marker","x","LineWidth",2)
            grid on
            xlabel("H (A/m)")
            ylabel("B (T)")
            title("B-H-curve")
        end

        function plotBMu(obj)
            b = 0:0.001:max(obj.B)*1.1;
            mur = obj.getMuNonlinear(b)/obj.Mu0;
            plot(b, mur, "Color", "black", LineWidth=1.5)
            hold on
            scatter(obj.B, obj.B./obj.H/obj.Mu0, "filled", "MarkerEdgeColor","black", "Marker","x","LineWidth",2)
            grid on
            xlabel("B (T)")
            ylabel("Mu_r")
            title("B-mu-curve")
        end

        function plotBNu(obj)
            b = 0:0.001:max(obj.B)*1.1;
            nu = obj.getNuNonlinear(b);
            plot(b, nu, "Color", "black", LineWidth=1.5)
            hold on
            scatter(obj.B, obj.H./obj.B, "filled", "MarkerEdgeColor","black", "Marker","x","LineWidth",2)
            grid on
            xlabel("B (T)")
            ylabel("Nu ")
            title("B-nu-curve")
        end

        function plotBNuPrime(obj)
            b = 0:0.001:max(obj.B)*1.1;
            nuP = obj.getNuPrimeNonlinear(b);
            plot(b, nuP, "Color", "black", LineWidth=1.5)
            grid on
            xlabel("B (T)")
            ylabel("dNu/dB ")
            title("B-nu'-curve")
        end

        function plotBNuPrimeB2(obj)
            % TBD NOT YET FINISHED
            b = 0:0.001:max(obj.B)*1.1;
            nuP = obj.getNuPrimeB2Nonlinear(b);
            plot(b, nuP, "Color", "black", LineWidth=1.5)
            grid on
            xlabel("B (T)")
            ylabel("dNu/dB^2 ")
            title("B^2-nu'-curve")
        end
    end
end