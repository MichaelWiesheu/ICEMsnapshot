function [x] = PMSM_constraints1 (parameters)
    % Standard Rotor Parameters
    MW          = 19e-3;            % Magnet width
    MH          = 7e-3;             % Magnet height
    MAG         = 7e-3;             % Distance from air gap to magnet

    
    % Standard Stator Parameters
    SLOTS = 36;             % nuber of nuts in the stator
    RD1 = 88e-3;           % outer diameter rotor
    SD1 = 135e-3;           % outer diameter stator
    SD2 = 90e-3;           % inner diameter stator
    
    Sr1 = 1.2e-3;           % radius of outer slot radius
%     Sr2 = 2.5e-3;           % radius of inner slot radius
    SW1 = 4e-3;             % Stator tooth width
    SW2 = 2.3e-3;           % Nut width
    SW3 = 0.64e-3;           % Nut depth
    SW4 = 8.25e-3;             % Stator back iron thickness
    % 
    if nargin == 1
        data_names = fieldnames (parameters);
        for iopts  = 1:numel (data_names)
          eval ([data_names{iopts} '= parameters.(data_names{iopts});']);
        end
    end

    % calculated parameters stator
    center = [0, 0, 0];
    rout_stat = SD1/2;                      % outher radius stator
    rout_rot = RD1/2;                       % outher radius rotor
    rin_stat = SD2/2;                       % inner radius stator
    r_air = (rout_rot+rin_stat)/2;          % Center air gap radius

    beta = 2*pi/SLOTS;
    alpha = beta/2;
    gamma1 = 2*asin(SW2/2 / rin_stat); 

    L2 = SD1/2 - SW4 - Sr1;           % lenght till center of circle 3

    % Equations: 
    % L1 = cos(gamma)*rin_stat + SW3 + Dslit
    % (SW2/2)^2 + Dslit^2 = Sr2^2
    % sin(alpha)*L1 = SW1/2 + Sr2
    const = SW1/2/sin(alpha) - cos(gamma1/2)*rin_stat - SW3;
    a = (1/sin(alpha)^2-1);
    b = (2/sin(alpha)*const);
    c = const^2 + (SW2/2)^2;
    Sr2 = (-b + sqrt(b^2-4*a*c))/2/a;
    Dslot = (Sr2^2 - (SW2/2)^2);
    L1 = (SW1/2 + Sr2)/sin(alpha);
    
    Lslot = sqrt(Sr2^2 - SW2^2/4);
    Rslot = sqrt((L1-Lslot)^2 + (SW2/2)^2);
    gamma2 = 2*asin(SW2/2 / Rslot); 
    epsilon = asin(SW2/2/Sr2);
    
    B = SW1/2 + Sr1;
    angleL2 = atan((L2*sin(alpha)-B)/(L2*cos(alpha)));
    L3 = L2*tan(angleL2);

    % patch definitions
    patches.Iron = [4, 6, 15:19];
    patches.Air = [1:3, 5];
    patches.Magnets = [];
    patches.Windings = [7:14];

    %setting points
    point{1} = [rin_stat*cos(alpha), -rin_stat*sin(alpha)];
    point{2} = [rin_stat*cos(gamma1/2), -rin_stat*sin(gamma1/2)];
    point{3} = [rin_stat*cos(gamma1/2), rin_stat*sin(gamma1/2)];
    point{4} = [rin_stat*cos(alpha), rin_stat*sin(alpha)];
    
    point{5} = [Rslot*cos(alpha), -Rslot*sin(alpha)];
    point{6} = [L1-Lslot, - SW2/2];
    point{7} = [L1-Lslot, SW2/2];
    point{8} = [Rslot*cos(alpha), Rslot*sin(alpha)];

    point{9} = [L1 - sin(alpha)*Sr2, -cos(alpha)*Sr2];
    point{10} = [L1, 0];
    point{11} = [L1 - sin(alpha)*Sr2, cos(alpha)*Sr2];

    point{12} = [L2 - sin(alpha)*Sr1, -L3 - cos(alpha)*Sr1];
    point{13} = [L2, -L3];
    point{14} = [L2, L3];
    point{15} = [L2 - sin(alpha)*Sr1, L3 + cos(alpha)*Sr1];

    point{16} = [L2+Sr1, -L3];
    point{17} = [L2+Sr1, L3];

    point{18} = rout_stat*[cos(alpha), -sin(alpha)];
    point{19} = rout_stat*[cos(alpha), sin(alpha)];

    point{20} = (point{5}+point{18})/2;
    point{21} = (point{9}+point{12})/2;
    point{22} = (point{11}+point{15})/2;
    point{23} = (point{8}+point{19})/2;

    point{20} = point{20} * vecmag(point{21})/vecmag(point{20});
    point{23} = point{23} * vecmag(point{22})/vecmag(point{23});

    pS = 0.6; % point share > 0.5 (centered)
    point{24} = pS*pS*point{9} + pS*(1-pS)*(point{11}) + pS*(1-pS)*point{21} + (1-pS)*(1-pS)*point{22};
    point{25} = pS*pS*point{11} + pS*(1-pS)*(point{9}) + pS*(1-pS)*point{22} + (1-pS)*(1-pS)*point{21};
    
    point{26} = pS*pS*point{12} + pS*(1-pS)*(point{15}) + pS*(1-pS)*point{21} + (1-pS)*(1-pS)*point{22};
    point{27} = pS*pS*point{15} + pS*(1-pS)*(point{22}) + pS*(1-pS)*point{12} + (1-pS)*(1-pS)*point{21};

    point{28} = pS*point{21}+ (1-pS)*point{22};
    point{29} = pS*point{22}+ (1-pS)*point{21};


    x(1) = 3*MW - 2*MAG - 50e-3;
    x(2) = MH + MAG - 15e-3;
    x(3) = - L3;
%     p7 = point{7}; p11 = point{11}; p15 = point{15};
%     x(2) = atan2(p7(2), p7(1)) - atan2(p11(2), p11(1));
    x(4) = 4*a*c - b^2;

    x = x * 100;
end


