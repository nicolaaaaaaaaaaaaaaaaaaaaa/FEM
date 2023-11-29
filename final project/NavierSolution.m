
%Solution for a uniform distributed load over a simply supported
%plate

x1 = 2.5; %Deflection coordinate
x2 = 2.5; %Delfection coordinate
a = 5; %Plate width
b = 5; %Plate height
t = 0.1; %Thickness

EX = 5*10^7;
EY = 1*10^7;
nu = 0.25;
GXY = 6.67*10^6;

p = 1; %Distributed load

%Calculate D entries
alpha = 1-(nu^2*EY/EX);
D1111 = (EX/alpha)*(t^3/12);
D1122 = (EX/alpha)*(nu*EY/EX)*(t^3/12);
D2222 = (EX/alpha)*(EY/EX)*(t^3/12);
D1212 = GXY*t^3/12;

%Summation over w terms
w = 0;
for m = 1:1000
    for n = 1:1000
        eta = m*b/(n*a);
        phi = eta^4*D1111+D2222+eta^2*2*(D1122+2*D1212);
        pmn = 16*p/(pi^2*m*n);
        wmn = pmn*b^4/(pi^4*n^4*phi);

        w = w + (wmn*sin(pi*m*x1/a)*sin(pi*n*x2/b));

    end
end

display(w)