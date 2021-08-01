DH = DHnew;

 forAll(mesh.cells(),celli)
    {

       if (DH[celli] > L.value())
          {
             DH[celli] = L.value();
          }
       if (DH[celli] < 0)
          {
             DH[celli] = 0;
          }

    };

 lf = DH/L;


