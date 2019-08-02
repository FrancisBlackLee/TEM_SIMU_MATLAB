a=1;% cubic cell size
L=a*[0 1 0 1 0 1 0 1 0.5 0.5 0.5 0 1 0.5 0.25 0.75 0.25 0.75;
     0 0 1 1 0 0 1 1 0.5 0.5 0 0.5 0.5 1 0.25 0.75 0.75 0.25;
     0 0 0 0 1 1 1 1 0 1 0.5 0.5 0.5 0.5 0.25 0.25 0.75 0.75;];
L1=[L;14*ones(1,18)];% the lattice structrue
t=[1 1 1];% projection direction
[Lrp1, SliceInfo1]=CrystalSlicing(L1,t,1);

L2=[L;31*ones(1,14) 33*ones(1,4)];
[Lrp2, SliceInfo2]=CrystalSlicing(L2,t,1);