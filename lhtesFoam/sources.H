DH = DHnew;

DH.max(0);
DH.min(L);

lf = mask*DH/L + (1-mask)*lf;

