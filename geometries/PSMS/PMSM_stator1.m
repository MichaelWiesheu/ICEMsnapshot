function [srf, patches, parameters] = PMSM_stator1(parameters)
    draw_geometry = true;
    % drive standard parameters
    Symmetry = 6;
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

    % overwrite standard parameters by input data
    if nargin == 1
        data_names = fieldnames (parameters);
        for iopts  = 1:numel (data_names)
          eval ([data_names{iopts} '= parameters.(data_names{iopts});']);
        end
    end
    
    % calculated parameters
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

    % creating surfaces
    % Air inside COUPLING INTERFACE IS HARD CODED FOR SAME BOUNDARY REPRESENTATION
    srf(1) = nrbruled(nrbcirc(r_air, center, -beta/2, -beta/6), nrbcirc(rin_stat, center, -beta/2, -gamma1/2));
    srf(2) = nrbruled(nrbcirc(r_air, center, -beta/6, beta/6), nrbcirc(rin_stat, center, -gamma1/2, gamma1/2));
    srf(3) = nrbruled(nrbcirc(r_air, center, beta/6, beta/2), nrbcirc(rin_stat, center, gamma1/2, beta/2));
    % Iron/air layer
    srf(4) = nrbruled(nrbcirc(rin_stat, center, -beta/2, -gamma1/2), nrbcirc(Rslot, center, -beta/2, -gamma2/2));
    srf(5) = nrbruled(nrbcirc(rin_stat, center, -gamma1/2, gamma1/2), nrbline(point{6}, point{7}));
    srf(6) = nrbruled(nrbcirc(rin_stat, center, gamma1/2, beta/2), nrbcirc(Rslot, center, gamma2/2, beta/2));

    % Copper Slot Top Left
    nrbTop1 = nrbdegelev(nrbline(point{22}, point{11}), 1);
    nrbTop2 = nrbcirc(Sr2, point{10}, pi/2+alpha, pi-epsilon);
    nrbTop = nrbglue(nrbTop1, nrbTop2, 2, 1);
    nrbTop.knots = nrbTop.knots/max(nrbTop.knots);
    nrbTop = nrbreverse(nrbTop);
    nrbLeft = nrbline(point{25}, point{7});
    nrbRight = nrbline(point{29}, point{22});
    nrbBot = nrbline(point{25}, point{29});
    srf(7) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);

    % Copper Slot Top Right
    nrbTop1 = nrbdegelev(nrbline(point{15}, point{22}), 1);
    nrbTop2 = nrbcirc(Sr1, point{14}, 0, pi/2+alpha);
    nrbTop = nrbglue(nrbTop1, nrbTop2, 1, 2);
    nrbTop.knots = nrbTop.knots/max(nrbTop.knots);
    nrbTop = nrbreverse(nrbTop);
    nrbLeft = nrbextract(srf(7), 2);
    nrbBot = nrbline(point{29}, point{27});
    nrbRight = nrbline(point{27}, point{17});
    srf(8) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);

    % Copper Slot Bot Left and Right
    mirror = eye(4);
    mirror(2, 2) = -1;
    srf(9) = nrbtform(srf(7), mirror);
    srf(10) = nrbtform(srf(8), mirror);
    
    % Copper Slot Left 
    nrbBot = nrbreverse(nrbextract(srf(9), 1));
    nrbTop = nrbreverse(nrbextract(srf(7), 1));
    nrbLeft = nrbextract(srf(5), 4);
    nrbRight = nrbline(point{24}, point{25});
    srf(11) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);
%     
    % Copper Slot Center Left
    nrbBot = nrbextract(srf(9), 3);
    nrbTop = nrbextract(srf(7), 3);
    nrbLeft = nrbextract(srf(11), 2);
    nrbRight = nrbline(point{28}, point{29});
    srf(12) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);

    % Copper Slot Center Right
    nrbBot = nrbextract(srf(10), 3);
    nrbTop = nrbextract(srf(8), 3);
    nrbLeft = nrbextract(srf(12), 2);
    nrbRight = nrbline(point{26}, point{27});
    srf(13) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);

    % Copper Slot Right
    nrbBot = nrbextract(srf(10), 2);
    nrbTop = nrbextract(srf(8), 2);
    nrbLeft = nrbextract(srf(13), 2);
    nrbRight = nrbline(point{16}, point{17});
    srf(14) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);

    % Iron Top Left
    nrbBot = nrbextract(srf(7), 4);
    nrbTop = nrbline(point{8}, point{23});
    nrbLeft = nrbextract(srf(6), 4);
    nrbRight = nrbline(point{22}, point{23});
    srf(15) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);

    % Iron Top Right
    nrbBot = nrbextract(srf(8), 4);
    nrbTop = nrbline(point{23}, point{19});
    nrbLeft = nrbextract(srf(15), 2);
    nrbRight = nrbline(point{17}, point{19});
    srf(16) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);

    % Bot iron
    srf(17) = nrbtform(srf(15), mirror);
    srf(18) = nrbtform(srf(16), mirror);

    % Right iron 
    nrbBot = nrbextract(srf(18), 2);
    nrbTop = nrbextract(srf(16), 2);
    nrbLeft = nrbextract(srf(14), 2);
    nrbRight = nrbcirc(rout_stat, center, -alpha, alpha);
    srf(19) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);
    
    % Copy and rotate
    nsrf = numel(srf);
    IronInit = patches.Iron;
    AirInit = patches.Air;
    WindingsInit = patches.Windings;
    for I = 1:SLOTS/Symmetry-1
        for i = 1:nsrf
            srf(I*nsrf+i) = nrbtform(srf(i), vecrotz(I*2*pi/SLOTS));
        end
        patches.Iron = [patches.Iron, IronInit+nsrf*I];
        patches.Air = [patches.Air, AirInit+nsrf*I];
        patches.Windings = [patches.Windings, WindingsInit+nsrf*I];
    end
    % Rotate again by half a slot and upward facing
    nsrf = numel(srf);
    for I = 1:nsrf
        srf(I) = nrbtform(srf(I), vecrotz(pi/SLOTS - pi/Symmetry + pi/2));
    end

    if (draw_geometry)
        figure()

        hold on
        for i = 1:numel(srf)
            if ismember(i, patches.Iron)
                nrbplotcol(srf(i), [50, 50], 'color', [0.3, 0.3, 0.3]);
            elseif ismember(i, patches.Air)
                nrbplotcol(srf(i), [50, 50], 'color', [0.1, 0.1, 1]);
            elseif ismember(i, patches.Windings)
                nrbplotcol(srf(i), [50, 50], 'color', [0.9, 0.1, 1]);
            else
                disp("stator patch is not material defined")
                disp(i)
                nrbplotcol(srf(i), [10, 10], 'color', [0.1, 0.1, 0.1]);
            end
        end

%         for i = 1:numel(point)
%             scatter(point{i}(1), point{i}(2), "black", "filled");
%             text(point{i}(1)+1e-5, point{i}(2), string(i))
%         end
        view(2)
        axis equal
    end

end