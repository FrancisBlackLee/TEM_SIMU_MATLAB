function [regParams,Bfit,ErrorStats]=absorient(A,B,doScale)
% This tool solves the absolute orientation problem, i.e., it finds the 
% rotation, translation, and optionally also the scaling, that best maps one
% collection of 3D point coordinates to another in a least squares sense. 
% Namely,
%  
%            [regParams,Bfit,ErrorStats]=absorientParams(A,B,doScale)
% 
% solves, when doScale=false (the default),
%  
%                min. sum_i ||R*A(:,i) + t - B(:,i)||^2
%  
% where R is a 3D rotation matrix and t is a translation vector.
%  
%
%When doScale=true, it solves the more general problem
% 
%                min. sum_i ||s*R*A(:,i) + t - B(:,i)||^2
% 
%where s is a global scale factor. The registration uses Horn's
%quaternion-based algorithm.
%
%
%in:
%
%  A: a 3xN matrix whos columns are the 3D coords of N source points.
%  B: a 3xN matrix whos columns are the 3D coords of N target points.
%  doScale: Boolean flag. If true (default=false), the registration will 
%          include a scale factor.
%          
%out:
%
%
% regParams: structure output with estimated registration parameters,
%
%     regParams.R:   The estimated rotation
%     regParams.t:   The estimated translation
%     regParams.s:   The estimated scale factor (set to 1 if doScale=false).
%     regParams.M:   4x4 homogenous coordinate transform matrix [s*R,t;[0 0 0 1]].
%     regParams.q:   A unit quaternion [q0 qx qy qz] corresponding to R and
%                    signed to satisfy max(q)=max(abs(q))>0
% 
%
%  Bfit: The rotation, translation, and scaling (as applicable) of A that 
%        best matches B in least squares sense.
%
%
% ErrorStats: structure output with error statistics. In particular,
%             defining err(i)=norm( Bfit(:,i)-B(:,i) ), it contains
%
%      ErrorStats.errlsq = 0.5* norm(err) 
%      ErrorStats.errmax = max(err) 
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Matt Jacobson
% Copyright, Xoran Technologies, Inc.  http://www.xorantech.com



ncols=@(M) size(M,2); %number of columns
matmvec=@(M,v) bsxfun(@minus,M,v); %matrix-minus-vector


nn=ncols(A);

if nargin<3, doScale=false; end

if nn~=ncols(B),
    error 'The number of points to be registered must be the same'
end



lc=mean(A,2);  rc=mean(B,2);  %Centroids

left  = matmvec(A,lc); %Center coordinates at centroids
right = matmvec(B,rc); 

M=left*right.';

[Sxx,Syx,Szx,  Sxy,Syy,Szy,   Sxz,Syz,Szz]=dealr(M(:));

N=[(Sxx+Syy+Szz)  (Syz-Szy)      (Szx-Sxz)      (Sxy-Syx);...
   (Syz-Szy)      (Sxx-Syy-Szz)  (Sxy+Syx)      (Szx+Sxz);...
   (Szx-Sxz)      (Sxy+Syx)     (-Sxx+Syy-Szz)  (Syz+Szy);...
   (Sxy-Syx)      (Szx+Sxz)      (Syz+Szy)      (-Sxx-Syy+Szz)];

[V,D]=eig(N);

[trash,emax]=max(real(  diag(D)  )); emax=emax(1);

q=V(:,emax); %Gets eigenvector corresponding to maximum eigenvalue
q=real(q);   %Get rid of imaginary part caused by numerical error

[trash,ii]=max(abs(q)); sgn=sign(q(ii(1)));
q=q*sgn; %Sign ambiguity

R=quatern2orth(q); %map to orthogonal matrix

if doScale
   
     summ = @(M) sum(M(:));
  
     sss=summ( right.*(R*left))/summ(left.^2);
     t=rc-R*(lc*sss);
     
     
else
    
    sss=1;
    t=rc-R*lc;
   
 
end



regParams.R=R;
regParams.t=t;
regParams.s=sss;
regParams.M=[sss*R,t;[0 0 0 1]];
regParams.q=q/norm(q);

if nargout>1
    
    Bfit=matmvec((sss*R)*A,-t);
    
end

if nargout>2
    
    l2norm = @(M,dim) sqrt(sum(M.^2,dim));

    err=l2norm(Bfit-B,1);
    
    ErrorStats.errlsq=0.5*norm(err);
    ErrorStats.errmax=max(err);
    
   
end
    


function R=quatern2orth(quat)
%Map a quaternion to an orthonormal 3D matrix
%
% R=quatern2orth(quat)
%
%in:
%
% quat: A quaternion [q0 qx qy qz]'
%
%out:
%
% R: The orthonormal 3D matrix induced by the 
%    unit quaternion quat/norm(quat)


quat=quat(:);
nrm=norm(quat);
if ~nrm
 'Quaternion distribution is 0'    
end

quat=quat./norm(quat);

q0=quat(1);
v =quat(2:4);
[qx,qy,qz]=dealr(v);

A=[q0 -qz qy;...
   qz q0 -qx;...
  -qy qx  q0 ];

R=v*v.' + A^2;



function varargout=dealr(v)

  varargout=num2cell(v);

