#include <stdlib.h>
#include <stdio.h>

#include "models.h"
#include "load.h"
#include "errors.h"


/*-----------------------------------------------------------------------
|		      	     Load a xyzi file into a stream
|
-----------------------------------------------------------------------*/
points *loadxyzi(FILE *estream, int i) {
  printf("step 2: Loading point cloud file\n");
  points *cloud;
  cloud = malloc(sizeof(double) * i * 3 + sizeof(int) * 2 * i);
  if (cloud == NULL) {
    error();
  }

  fseek(estream, 0, SEEK_SET);

  char *del = "%lf %lf %lf %d\n";

  int g;
  for (g = 0; g < i; g++) {
    fscanf(estream, del, &cloud[g].Rx, &cloud[g].Ry, &cloud[g].Rz,
           &cloud[g].IN);
  }
  for (g = 0; g < i; g++) {
    cloud[g].Rx = cloud[g].Rx * 1000;
    cloud[g].Ry = cloud[g].Ry * 1000;
    cloud[g].Rz = cloud[g].Rz * 1000;
  }
  fclose(estream);
  return (cloud);
}
/*-----------------------------------------------------------------------
|		      	 Load a orientation file into a stream
|
-----------------------------------------------------------------------*/
info orientation(FILE *ostream) {
  printf("step 3: Loading orientation file\n");
  info orient;

  fseek(ostream, 0, SEEK_SET);

  char *del = "%6s %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf\n";

  fscanf(ostream, del, &orient.name, &orient.c, &orient.xo, &orient.yo,
         &orient.k1, &orient.k2, &orient.imagewidth, &orient.imageheight,
         &orient.Xo, &orient.Yo, &orient.Zo, &orient.omega, &orient.phi,
         &orient.kapa);

  fclose(ostream);
  return (orient);
}
/*-----------------------------------------------------------------------
|				Load a xyzipxpy file into a stream
|
-----------------------------------------------------------------------*/
values *loadt(FILE *istream, int i) {
  printf("step 2: Loading point cloud file\n");
  values *project;
  project = malloc(sizeof(double) * i * 3 + sizeof(int) * 4 * i);
  if (project == NULL) {
    error();
  }

  fseek(istream, 0, SEEK_SET);

  char *del = "%lf %lf %lf %d %d %d\n";

  int g;
  for (g = 0; g < i; g++) {
    fscanf(istream, del, &project[g].Rx, &project[g].Ry, &project[g].Rz,
           &project[g].IN, &project[g].Px, &project[g].Py);
  }

  fclose(istream);
  return (project);
}


