classdef Mat_BrauerIron < Mat_Magnetic & Mat_Thermal
    properties
        k1 = 0.3774
        k2 = 2.970
        k3 = 388.33
        Bclip = 10;
        muClip;
%         Hclip;
    end

    methods (Access = public)
        function obj = Mat_BrauerIron()
            obj.Mur = 1200;
            obj.type = "nonlinear";
            obj.Sigma = 1.03e7;
            obj.Rho = 7.874;
            obj.correctBclip();
            
        end
        % Calculate the magetic flux density, where a further decreasing of
        % the differential permeability is not physical (Bclip)
        % permeabilities higher than Bclip will increase linearly with mu0
        function correctBclip(obj)
            b = 0:0.001:3;
            h = b.*obj.getNuNonlinear(b);
            mudiff = gradient(b, h);
            obj.Bclip = b(find(mudiff<obj.Mu0, 1));
            obj.muClip = obj.getMuNonlinear(obj.Bclip);
            obj.B = 0:0.1:obj.Bclip;
            obj.H = obj.B.*obj.getNuNonlinear(obj.B);
%             options = optimoptions('fsolve','Display','off');
%             obj.Bclip = fsolve(@(B)obj.getMuNonlinear(B)-obj.Mu0, 0, options);
        end

        function nu = getNuNonlinear(obj, B)
            nu = obj.k1*exp(obj.k2.*B.^2) + obj.k3;
            % correction for too high unphysical values (mu continues linearly)
            nu(B>obj.Bclip) = (B(B>obj.Bclip)-obj.Bclip.*(1-obj.Mu0./obj.muClip))./(B(B>obj.Bclip).*obj.Mu0);
        end
        function mu = getMuNonlinear(obj, B)
            mu = 1./obj.getNuNonlinear(B);
        end

%         function nuPrime2 = getNuPrime2Nonlinear(obj, B)
%             nuPrime2 = obj.k1.*obj.k2.*exp(obj.k2.*B.^2);
%         end
% 
        function nuPrime = getNuPrimeNonlinear(obj, B)
            nuPrime = 2*B.*obj.k1.*obj.k2.*exp(obj.k2.*B.^2);
            B2correct = B(B>obj.Bclip);
            nuPrime(B>obj.Bclip) = (B2correct./obj.Mu0 - obj.getNuNonlinear(B2correct).*B2correct)./B2correct.^2; 
            
%             Step = 1e-8;
% %             B(B==0) = 0.01;
%             nu1 = obj.getNuNonlinear(B);
%             nu2 = obj.getNuNonlinear(B + Step);
%             nuPrime = (nu2 - nu1) / Step;
        end
    end
end