DH = DHnew;

DH.max(0);
DH.min(L);

lf = pos(mask)*DH/L + pos(1-mask)*lf;
lf.correctBoundaryConditions();
