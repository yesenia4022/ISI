function [Sper] = getRetinaGradient(Vdom)

Mmax = 80;
Ms = 3.2;
Mn = 6.4;
Mmin = 0.8;

diam = 3.32;
Circ = 2*pi*diam/2;
p = -Vdom/360*Circ; %Dorsoventral axis, ventral to dorsal (mm)

Mper = Mmax*((p+2).^Mn) ./ ( (p+2).^Mn + Ms^Mn ) + Mmin;
Sper = 100-Mper;

