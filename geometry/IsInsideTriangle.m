function [t] = IsInsideTriangle(v1, v2, v3, u, thr)

if nargin == 4
    thr = 0.01;
end

v12 = [v2(1) - v1(1); v2(2) - v1(2)];
v13 = [v3(1) - v1(1); v3(2) - v1(2)];
v1u = [u(1) - v1(1); u(2) - v1(2)];

detv12 = det([v12, v13]);
detuv2 = det([v1u, v13]);
detuv1 = det([v1u, v12]);

a = detuv2 / detv12;
b = -detuv1 / detv12;

if a >= -thr && b >= -thr && (a + b) <= 1 + thr
    t = true;
else
    t = false;
end

end