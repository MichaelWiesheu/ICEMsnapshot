% bh-curve 

function mu = mu_hyst_Iron(B, material)
    %mu = interp1(material.BsIron, material.MusIron, B, "linear", "extrap");
    B(B==0) = 0.1;
    mu = B./interp1(material.BsIron,material.HsIron, B,"linear", "extrap");
    
%     mu = max(mu, material.MuVac);
%     mu = mu.*material.MuVac;
end