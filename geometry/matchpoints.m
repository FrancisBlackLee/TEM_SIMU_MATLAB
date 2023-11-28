function [permIndices,Bsorted]=matchpoints(A,B)
%                 
%IN:              
%                 
%    A: an Nx2 matrix of points         
%    B: an Nx2 matrix of points from A transformed and unordered          
%                 
%OUT:             
%                 
%    permIndices: permutation indices of rows of B thought to match A 
%        Bsorted: the Nx2 permuted version of B
    La=landmarks(A);
    Lb=landmarks(B);
    B3=B.'; B3(3,:)=0;
    reg=absor( Lb,La,'doScale',1);
    C3=(reg.s*reg.R)*B3+reg.t;
    C=C3(1:2,:).';
    N=size(A,1);
    e=ones(1,N);
    E=speye(N);
    Aeq=[ kron(E,e) ; kron(e,E) ]; beq=[e,e].';
    Q=kron(E,C.');
    d=reshape(A.',[],1);
    lb=zeros(1,N^2);
    ub=ones(1,N^2);
    intcon=1:N^2;
    
    P=minL1intlin(Q,d,intcon,[],[],Aeq,beq,lb,ub);
    P=round(reshape(P,N,[]));
    permIndices=(1:N)*P;
    
    if nargout>1
     Bsorted=B(permIndices,:);
    end
            function L=landmarks(P)
 
                G=pdist2(P,P); G(~G)=nan;
                
                [i,j]=find( G==min(G(:)) ,1);
                I=P(i,:);
                J=P(j,:);
                K=mean(P,1);
                
                if norm(I-K)<norm(J-K)
                    [I,J]=deal(J,I);
                end
                L=[I;J;K].';
                L(3,:)=0;
                
            end
end