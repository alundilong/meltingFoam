DH = DHnew;

DH.max(0);
DH.min(L);

//forAll(mesh.cells(),celli)
//{
//   if(DH[celli] > L.value())
//   {
//       DH[celli] = L.value();
//   }
//
//   if(DH[celli] < 0)
//   {
//       DH[celli] = 0;
//   }
//};

lf = mask*DH/L + (1-mask)*lf;
//lf = DH/L;

