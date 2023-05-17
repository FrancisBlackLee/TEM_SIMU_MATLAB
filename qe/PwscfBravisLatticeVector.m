function [a1, a2, a3, b1, b2, b3] = PwscfBravisLatticeVector(ibrav, varargin)
%PwscfBravisLatticeVector returns the Bravis lattice vectors for the given
%ibrav following PWSCF definition, NOTE that ibrav=0 is not supported.
%   PwscfBravisLatticeVector(ibrav, celldm)
%   PwscfBravisLatticeVector(ibrav, A, B, C, cosAB, cosAC, cosBC)

if nargin == 2
    celldm = varargin{1};
    switch ibrav
        case 1
            a = celldm(1);
            a1 = a * [1,0,0];
            a2 = a * [0,1,0];
            a3 = a * [0,0,1];
        case 2
            a = celldm(1);
            a1 = (a/2) * [-1,0,1];
            a2 = (a/2) * [0,1,1]; 
            a3 = (a/2) * [-1,1,0];
        case 3
            a = celldm(1);
            a1 = (a/2) * [1,1,1];
            a2 = (a/2) * [-1,1,1];
            a3 = (a/2) * [-1,-1,1];
        case -3
            a = celldm(1);
            a1 = (a/2) * [-1,1,1];
            a2 = (a/2) * [1,-1,1];
            a3 = (a/2) * [1,1,-1];
        case 4
            a = celldm(1);
            c = celldm(3) * celldm(1);
            a1 = a * [1,0,0];
            a2 = a * [-1/2,sqrt(3)/2,0];
            a3 = a * [0,0,c/a];
        case 5
            a = celldm(1);
            c = celldm(4);
            tx=sqrt((1-c)/2);
            ty=sqrt((1-c)/6);
            tz=sqrt((1+2 * c)/3);
            a1 = a * [tx,-ty,tz];
            a2 = a * [0,2 * ty,tz];
            a3 = a * [-tx,-ty,tz];
        case -5
            a = celldm(1);
            ap = a / sqrt(3);
            c = celldm(4);
            ty=sqrt((1-c)/6);
            tz=sqrt((1+2 * c)/3);
            u = tz - 2*sqrt(2)*ty;
            v = tz + sqrt(2)*ty;
            a1 = ap * [u,v,v];
            a2 = ap * [v,u,v];
            a3 = ap * [v,v,u];
        case 6
            a = celldm(1);
            a1 = a * [1,0,0];
            a2 = a * [0,1,0];
            a3 = a * [0,0,celldm(3)];
        case 7
            a = celldm(1);
            a1 = (a/2) * [1,-1,celldm(3)];
            a2 = (a/2) * [1,1,celldm(3)];
            a3 = (a/2) * [-1,-1,celldm(3)];
        case 8
            a = celldm(1);
            b = a * celldm(2);
            c = a * celldm(3);
            a1 = [a,0,0];
            a2 = [0,b,0];
            a3 = [0,0,c];
        case 9
            a = celldm(1);
            b = a * celldm(2);
            c = a * celldm(3);
            a1 = [a/2, b/2,0];
            a2 = [-a/2,b/2,0];
            a3 = [0,0,c];
        case -9
            a = celldm(1);
            b = a * celldm(2);
            c = a * celldm(3);
            a1 = [a/2,-b/2,0];
            a2 = [a/2, b/2,0];
            a3 = [0,0,c];
        case 91
            a = celldm(1);
            b = a * celldm(2);
            c = a * celldm(3);
            a1 = [a, 0, 0];
            a2 = [0,b/2,-c/2];
            a3 = [0,b/2,c/2];
        case 10
            a = celldm(1);
            b = a * celldm(2);
            c = a * celldm(3);
            a1 = [a/2,0,c/2];
            a2 = [a/2,b/2,0];
            a3 = [0,b/2,c/2];
        case 11
            a = celldm(1);
            b = a * celldm(2);
            c = a * celldm(3);
            a1 = [a/2,b/2,c/2];
            a2 = [-a/2,b/2,c/2];
            a3 = [-a/2,-b/2,c/2];
        case 12
            a = celldm(1);
            b = a * celldm(2);
            c = a * celldm(3);
            gamma = acos(celldm(4));
            a1 = [a,0,0];
            a2 = [b*cos(gamma),b*sin(gamma),0];
            a3 = [0,0,c];
        case -12
            a = celldm(1);
            b = a * celldm(2);
            c = a * celldm(3);
            beta = acos(celldm(5));
            a1 = [a,0,0];
            a2 = [0,b,0];
            a3 = [c*cos(beta),0,c*sin(beta)];
        case 13
            a = celldm(1);
            b = a * celldm(2);
            c = a * celldm(3);
            gamma = acos(celldm(4));
            a1 = [a/2, 0, -c/2];
            a2 = [b*cos(gamma), b*sin(gamma), 0];
            a3 = [a/2, 0, c/2];
        case -13
            a = celldm(1);
            b = a * celldm(2);
            c = a * celldm(3);
            beta = acos(celldm(5));
            a1 = [a/2, b/2, 0];
            a2 = [-a/2, b/2, 0];
            a3 = [c*cos(beta), 0, c*sin(beta)];
        case 14
            a = celldm(1);
            b = a * celldm(2);
            c = a * celldm(3);
            alpha = acos(celldm(4));
            beta = acos(celldm(5));
            gamma = acos(celldm(6));
            a1 = [a, 0, 0];
            a2 = [b*cos(gamma), b*sin(gamma), 0];
            a3 = [c*cos(beta),  c*(cos(alpha)-cos(beta) * cos(gamma))/sin(gamma),...
                c*sqrt( 1 + 2*cos(alpha) * cos(beta) * cos(gamma) -...
                cos(alpha)^2-cos(beta)^2-cos(gamma)^2 )/sin(gamma) ];
        otherwise
            error('Invalid input');
    end
elseif nargin == 7
    a = varargin{1};
    b = varargin{2};
    c = varargin{3};
    cosab = varargin{4};
    cosac = varargin{5};
    cosbc = varargin{6};
    gamma = acos(cosab);
    beta = acos(cosac);
    alpha = acos(cosbc);
    switch ibrav
        case 1
            a1 = a * [1,0,0];
            a2 = a * [0,1,0];
            a3 = a * [0,0,1];
        case 2
            a1 = (a/2) * [-1,0,1];
            a2 = (a/2) * [0,1,1]; 
            a3 = (a/2) * [-1,1,0];
        case 3
            a1 = (a/2) * [1,1,1];
            a2 = (a/2) * [-1,1,1];
            a3 = (a/2) * [-1,-1,1];
        case -3
            a1 = (a/2) * [-1,1,1];
            a2 = (a/2) * [1,-1,1];
            a3 = (a/2) * [1,1,-1];
        case 4
            a1 = a * [1,0,0];
            a2 = a * [-1/2,sqrt(3)/2,0];
            a3 = a * [0,0,c/a];
        case 5
            tx=sqrt((1-cosab)/2);
            ty=sqrt((1-cosab)/6);
            tz=sqrt((1+2 * cosab)/3);
            a1 = a * [tx,-ty,tz];
            a2 = a * [0,2 * ty,tz];
            a3 = a * [-tx,-ty,tz];
        case -5
            ap = a / sqrt(3);
            ty=sqrt((1-cosab)/6);
            tz=sqrt((1+2 * cosab)/3);
            u = tz - 2*sqrt(2)*ty;
            v = tz + sqrt(2)*ty;
            a1 = ap * [u,v,v];
            a2 = ap * [v,u,v];
            a3 = ap * [v,v,u];
        case 6
            a1 = a * [1,0,0];
            a2 = a * [0,1,0];
            a3 = c * [0,0,1];
        case 7
            a1 = (a/2) * [1,-1,c/a];
            a2 = (a/2) * [1,1,c/a];
            a3 = (a/2) * [-1,-1,c/a];
        case 8
            a1 = [a,0,0];
            a2 = [0,b,0];
            a3 = [0,0,c];
        case 9
            a1 = [a/2, b/2,0];
            a2 = [-a/2,b/2,0];
            a3 = [0,0,c];
        case -9
            a1 = [a/2,-b/2,0];
            a2 = [a/2, b/2,0];
            a3 = [0,0,c];
        case 91
            a1 = [a, 0, 0];
            a2 = [0,b/2,-c/2];
            a3 = [0,b/2,c/2];
        case 10
            a1 = [a/2,0,c/2];
            a2 = [a/2,b/2,0];
            a3 = [0,b/2,c/2];
        case 11
            a1 = [a/2,b/2,c/2];
            a2 = [-a/2,b/2,c/2];
            a3 = [-a/2,-b/2,c/2];
        case 12
            a1 = [a,0,0];
            a2 = [b*cos(gamma),b*sin(gamma),0];
            a3 = [0,0,c];
        case -12
            a1 = [a,0,0];
            a2 = [0,b,0];
            a3 = [c*cos(beta),0,c*sin(beta)];
        case 13
            a1 = [a/2, 0, -c/2];
            a2 = [b*cos(gamma), b*sin(gamma), 0];
            a3 = [a/2, 0, c/2];
        case -13
            a1 = [a/2, b/2, 0];
            a2 = [-a/2, b/2, 0];
            a3 = [c*cos(beta), 0, c*sin(beta)];
        case 14
            a1 = [a, 0, 0];
            a2 = [b*cos(gamma), b*sin(gamma), 0];
            a3 = [c*cos(beta),  c*(cos(alpha)-cos(beta) * cos(gamma))/sin(gamma),...
                c*sqrt( 1 + 2*cos(alpha) * cos(beta) * cos(gamma) -...
                cos(alpha)^2-cos(beta)^2-cos(gamma)^2 )/sin(gamma) ];
        otherwise
            error('Invalid input');
    end
else
    error('Invalid input');
end

[b1, b2, b3] = DirectToReciprocal(a1, a2, a3);

end