#include "mesh.h"

int main(int argc, char* argv[])
{
  int    iTri, iVer;
  double to, ti;

  if (argc < 2) {
    printf(" usage : mesh file \n");
    return 0;
  }

  //--- read a mesh
  to        = clock();
  Mesh* Msh = msh_read(argv[1], 0);
  ti        = clock();

  if (!Msh)
    return 0;

  printf("  Vertices   %10d \n", Msh->NbrVer);
  printf("  Triangles  %10d \n", Msh->NbrTri);
  printf("  time to read the mesh %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  //--- create neigbhors Q2 version
  to = clock();
  msh_neighborsQ2(Msh);
  ti = clock();
  printf("  time q2 neigh.        %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  //--- create neigbhors with hash table
  to = clock();
  msh_neighbors(Msh);
  printf(" nobmre tri %d \n",Msh->NbrTri);
  ti = clock();
  printf("  time hash tab neigh.  %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  //--- TODO: compute mesh quality
 double* Qal = (double*)malloc(sizeof(double) * (Msh->NbrTri + 1));

  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    //Qal[iTri] = (double)iTri / 10.;
    Qal[iTri]=Q1(iTri,Msh);
  }


  // Compter le nombre d'arêtes frontières
  int NbBorderEdge=0;
  for (int iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    int count_neighbours=0;
    for (int iEdg = 0; iEdg < 3; iEdg++) {
      if (Msh->TriVoi[iTri][iEdg]<= Msh->NbrTri &&Msh->TriVoi[iTri][iEdg]>0 )
      {
        count_neighbours+=1;
      }
      
    }
    if (count_neighbours<3)
    {
      NbBorderEdge+=1;
    }
    
  }
   printf("  nombre d'arêtes frontières %d  \n", NbBorderEdge);
  msh_write2dfield_Triangles("quality.solb", Msh->NbrTri, Qal);

  //--- TODO: compute metric field
  double3d* Met = (double3d*)malloc(sizeof(double3d) * (Msh->NbrVer + 1));

  for (iVer = 1; iVer <= Msh->NbrVer; iVer++) {
    Met[iVer][0] = 1.;
    Met[iVer][1] = 0.;
    Met[iVer][2] = 1.;
  }

  msh_write2dmetric("metric.solb", Msh->NbrVer, Met);

  //--- Free memory
  if (Qal != NULL) {
    free(Qal);
    Qal = NULL;
  }
  if (Met != NULL) {
    free(Met);
    Met = NULL;
  }

  return 0;
}
