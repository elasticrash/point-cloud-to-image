#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "models.h"
#include "rot.h"
#include "errors.h"
#include "clines.h"
#include "load.h"

#define PI 3.1415926535898

/*-----------------------------------------------------------------------
|                         Function Prototypes |
-----------------------------------------------------------------------*/
int linesoffilex(FILE *stream);
values *loadfile(FILE *stream, int i);
values *sorttod(values real[], int i);
int threshold(values sort[], int i, info orient);
int cmp(const void *a, const void *b);
double collinearity(double xo, double c, double XA, double XO, double YA,
                    double YO, double ZA, double ZO, double omega, double phi,
                    double kapa);

/*-----------------------------------------------------------------------
|                              grad to rad
|
-----------------------------------------------------------------------*/
double gtr(double r) { return ((r / 200) * PI); }
/*-----------------------------------------------------------------------
|                              Collinearity
|
-----------------------------------------------------------------------*/
double collinearity(double xo, double c, double XA, double XO, double YA,
                    double YO, double ZA, double ZO, double omega, double phi,
                    double kapa) {
  double o = gtr(omega);
  double f = gtr(phi);
  double k = gtr(kapa);

  double xa;
  xa =
      xo -
      c * ((XA - XO) * r11(f, k) + (YA - YO) * r21(f, k) + (ZA - ZO) * r31(f)) /
          ((XA - XO) * r13(o, f, k) + (YA - YO) * r23(o, f, k) +
           (ZA - ZO) * r33(o, f));
  return (xa);
}
double collinearityy(double yo, double c, double XA, double XO, double YA,
                     double YO, double ZA, double ZO, double omega, double phi,
                     double kapa) {
  double o = gtr(omega);
  double f = gtr(phi);
  double k = gtr(kapa);

  double ya;
  ya = yo - c *
                ((XA - XO) * r12(o, f, k) + (YA - YO) * r22(o, f, k) +
                 (ZA - ZO) * r32(o, f)) /
                ((XA - XO) * r13(o, f, k) + (YA - YO) * r23(o, f, k) +
                 (ZA - ZO) * r33(o, f));
  return (ya);
}
/*-----------------------------------------------------------------------
|                         Calculate distance
|
-----------------------------------------------------------------------*/
double distance(double xa, double ya, double za, double xb, double yb,
                double zb) {
  double S;
  S = sqrt(pow(xb - xa, 2) + pow(yb - ya, 2) + pow(zb - za, 2));
  return (double)S;
}
/*-----------------------------------------------------------------------
|                   Compare arrays with doubles
|
-----------------------------------------------------------------------*/
int cmp(const void *a, const void *b) {
  const T *x = a, *y = b;
  if (*x > *y) {
    return 1;
  }
  if (*x < *y) {
    return -1;
  }
  return 0;
}
/*-----------------------------------------------------------------------
|                          Locate search window
|
-----------------------------------------------------------------------*/
int locate(values real[], int choise, int pos, int i) {
  int g;
  int *sarray;
  sarray = malloc(sizeof(int) * i);
  if (sarray == NULL) {
    error();
  }

  if (choise == 0) {
    for (g = 0; g < i; g++) {
      sarray[g] = real[g].Px;
    }
  } else if (choise == 1) {
    for (g = 0; g < i; g++) {
      sarray[g] = real[g].Py;
    }
  }

  qsort(sarray, i, sizeof(int), &cmp);

  int result = sarray[pos];

  return sarray[pos];
}
/*-----------------------------------------------------------------------
|					    sort a two d image array
|
-----------------------------------------------------------------------*/
values *sorttod(values real[], int i) {
  int g, k, l, o, b;
  o = 0;

  int rmin = locate(real, 0, 0, i);
  int cmin = locate(real, 1, 0, i);

  int rmax = locate(real, 0, i - 1, i);
  int cmax = locate(real, 1, i - 1, i);

  printf("pixel bounding box %d %d %d %d\n", rmin, rmax, cmin, cmax);
  
  values *temp;
  temp = malloc((sizeof(double) * i) * 3 + (sizeof(int) * i) * 4);

  if (temp == NULL) {
    error();
  }

  printf("step 5: sorting rows\n");
  for (g = rmin; g < rmax + 1; g++) {
    for (l = 0; l < i; l++) {
      if (real[l].Px == g) {
        temp[o].Rx = real[l].Rx;
        temp[o].Ry = real[l].Ry;
        temp[o].Rz = real[l].Rz;
        temp[o].IN = real[l].IN;
        temp[o].Px = real[l].Px;
        temp[o].Py = real[l].Py;

        o++;
      }
	}
  }

  free(real);

  o = 0;

  printf("step 6: sorting columns\n");

  values *sort;
  sort = malloc((sizeof(double) * i) * 3 + (sizeof(int) * i) * 4);

  if (sort == NULL) {
    error();
  }

  for (k = 0; k < i; k = k + g) {
    g = 0;
    while (temp[k + g].Px == temp[k + g + 1].Px) {
      g++;
    }
    g++;
    for (b = cmin; b < cmax + 1; b++) {
      for (l = 0; l < g; l++) {
        if (temp[k + l].Py == b) {
          sort[o].Rx = temp[k + l].Rx;
          sort[o].Ry = temp[k + l].Ry;
          sort[o].Rz = temp[k + l].Rz;
          sort[o].IN = temp[k + l].IN;
          sort[o].Px = temp[k + l].Px;
          sort[o].Py = temp[k + l].Py;
          o++;
        } else {
        }
      }
    }
  }

  free(temp);
  return (sort);
}
/*-----------------------------------------------------------------------
|                   save a xyzi file into a stream
|
-----------------------------------------------------------------------*/
int threshold(values sort[], int i, info orient) {
  int g, k, l, n;
  g = 0;
  n = 0;
  int t = 0;
  FILE *wstream = fopen("exportx", "w");

  printf("step 7: threshold operations \n");

  for (l = 0; l < i - 1; l = l + g) {
    g = 1;
    while (sort[l + g].Px == sort[l + g + 1].Px &&
           sort[l + g].Py == sort[l + g + 1].Py) {
      g++;
    }

    double *tempz;
    tempz = malloc(sizeof(double) * g);
    if (tempz == NULL) {
      error();
    }

    for (k = 0; k < g; k++) {
      tempz[k] = distance(orient.Xo, orient.Yo, orient.Zo, sort[l + k].Rx,
                          sort[l + k].Ry, sort[l + k].Rz);
    }

    qsort(tempz, g, sizeof(tempz[g]), &cmp);

    for (t = 0; t < g; t++) {
      if (tempz[0] + 0.01 > distance(orient.Xo, orient.Yo, orient.Zo,
                                     sort[l + t].Rx, sort[l + t].Ry,
                                     sort[l + t].Rz)) {
        fprintf(wstream, "%lf %lf %lf %d %d %d\n", sort[l + t].Rx,
                sort[l + t].Ry, sort[l + t].Rz, sort[l + t].IN, sort[l + t].Px,
                sort[l + t].Py);
        n++;
      }
    }
    free(tempz);
  }
  free(sort);
  fclose(wstream);
  return n;
}
/*-----------------------------------------------------------------------
|                   save a xyzi file into a stream
|
-----------------------------------------------------------------------*/
values *project(points cloud[], int i, info orient) {
  int g;
  FILE *wstream = fopen("coord", "w");
  printf("step 4: projecting points to plane\n");
  values *real;
  real = malloc((sizeof(double) * i) * 3 + (sizeof(int) * i) * 4);
  if (real == NULL) {
    error();
  }

  for (g = 0; g < i; g++) {

    double xa = collinearity(orient.xo, orient.c, cloud[g].Rx, orient.Xo,
                             cloud[g].Ry, orient.Yo, cloud[g].Rz, orient.Zo,
                             orient.omega, orient.phi, orient.kapa);
    double ya = collinearityy(orient.yo, orient.c, cloud[g].Rx, orient.Xo,
                              cloud[g].Ry, orient.Yo, cloud[g].Rz, orient.Zo,
                              orient.omega, orient.phi, orient.kapa);

    real[g].Rx = cloud[g].Rx;
    real[g].Ry = cloud[g].Ry;
    real[g].Rz = cloud[g].Rz;
    real[g].IN = cloud[g].IN;
    real[g].Px = xa + orient.imagewidth;
    real[g].Py = -ya + orient.imageheight;
    real[g].Pz = 1;

    fprintf(wstream, "%lf %lf %lf %d %d %d\n", cloud[g].Rx, cloud[g].Ry,
            cloud[g].Rz, cloud[g].IN, real[g].Px, real[g].Py);
  }

  fclose(wstream);

  return (real);
}
/*-----------------------------------------------------------------------
|							save image
|
-----------------------------------------------------------------------*/
void saveimage(values project[], info orient) {
  FILE *istream;
  istream = fopen("image.raw", "wb");
  errorf(istream);

  int v;
  int x, y;
  x = orient.imagewidth * 2;
  y = orient.imageheight * 2;

  unsigned char *image;
  image = malloc(sizeof(unsigned char) * x * y);
  if (image == NULL) {
    error();
  }

  int i, k, l;
  l = 0;
  int res = 0;
  int resx = 0;
  int resy = 0;
  int m = 0;
  for (i = 0; i < x * y; i++) {
    m = m + l;
    l = 0;
    v = 0;

    if (res == y) {
      resx++;
      res = 0;
    }
    resy = i - (resx * y);
    res++;

    while (project[m + l].Px == resx && project[m + l].Py == resy) {
      l++;
      v = 1;
    }
    if (v == 0) {
      l = 0;
      image[i] = 0;
    } else if (v == 1) {
      l++;
      double *tempz;
      tempz = malloc(sizeof(double) * l);
      if (tempz == NULL) {
        error();
      }

      for (k = 0; k < l; k++) {
        tempz[k] = distance(orient.Xo, orient.Yo, orient.Zo, project[m + k].Rx,
                            project[m + k].Ry, project[m + k].Rz);
      }
      qsort(tempz, l, sizeof(double), &cmp);

      for (k = 0; k < l; k++) {
        if (distance(orient.Xo, orient.Yo, orient.Zo, project[m + k].Rx,
                     project[m + k].Ry, project[m + k].Rz) == tempz[0]) {
          image[i] = project[m + k].IN;
        }
      }
    }
  }
  printf("\n");
  for (i = 0; i < x * y; i++) {
    fwrite(&image[i], sizeof(unsigned char), 1, istream);
  }
  fclose(istream);
}
/*-----------------------------------------------------------------------
|                                 MAIN
|
-----------------------------------------------------------------------*/
int main() {
  FILE *ostream;
  ostream = fopen("orientation", "r");
  errorf(ostream);

  FILE *estream;
  estream = fopen("export1.xyzi", "r");
  errorf(estream);

  int i = linesoffilex(estream);

  points *cloud;
  cloud = malloc((sizeof(double) * i) * 3 + (sizeof(int) * i));
  if (cloud == NULL) {
    error();
  }
  cloud = loadxyzi(estream, i);

  info orient;
  orient = orientation(ostream);

  values *real;
  real = malloc((sizeof(double) * i) * 3 + (sizeof(int) * i) * 4);
  if (real == NULL) {
    error();
  }

  real = project(cloud, i, orient);


  values *sort;
  sort = malloc((sizeof(double) * i) * 3 + (sizeof(int) * i) * 4);
  if (sort == NULL) {
    error();
  }

  sort = sorttod(real, i);

  int j = threshold(sort, i, orient);
  free(cloud);

  FILE *istream;
  istream = fopen("exportx", "r");
  errorf(istream);

  values *project;
  project = malloc(sizeof(double) * i * 3 + sizeof(int) * 3 * i);
  if (project == NULL) {
    error();
  }

  project = loadt(istream, j);

  saveimage(project, orient);
  
  printf("finished\n");

  free(project);
  fclose(ostream);
  fclose(estream);
  fclose(istream);
  return 0;
}
