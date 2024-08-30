% hysteresis differential dMu

function muInv = dmu_dB_iron(B, material)
    Step = 1e-8;
    B(B==0) = 0.01;
    mu1 = B./interp1(material.BsIron,material.HsIron, B,"linear", "extrap");
    mu2 = (B+Step)./interp1(material.BsIron,material.HsIron, B+Step,"linear", "extrap");

    muInv = (mu2.^-1-mu1.^-1)/Step;
%     muInv = interp1(material.BsIron, material.MusInvdBIron, B, "linear", "extrap");
%     muInv = min(muInv, 1);
%     muInv = muInv./material.MuVac;
end
