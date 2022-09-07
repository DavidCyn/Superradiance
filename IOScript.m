clear
format long
threshold = 0.000001;

%Radial logarithmic resolution and range
input.dx = 0.0005;
input.xMin = -30;
input.xMax = 10;

%Angular resolution
input.dy = 0.0005;

%Orbital quantum numbers
input.n = 0;
input.l = 1;
input.m = 1;

%Black hole mass and spin
input.M = 1;
input.a = 0.998;

%particle mass
input.mu = 0.49;


Shift = inf;
while abs(Shift) > threshold
    outputR = Radial(input);
    outputA = Angular(input);
    input.domegan = outputR.domegan;
    input.omegan = outputR.omegan;
    input.dLambda = outputA.dLambda;
    input.Lambda = outputA.Lambda;
    Shift = abs((outputR.domegan/outputR.omegan)) + abs(outputA.dLambda);
end
hold on
outputR.omegan
plot(outputR.x,log10(abs(real(outputR.R(1:end/2)))))
plot(outputR.x,log10(abs(imag(outputR.R(1:end/2)))))
plot(outputR.x,log10(abs(outputR.R(1:end/2))))
