#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.1415926535898

typedef struct
{
	char name[6];
	double c, xo, yo;
	double k1, k2;
	int imagewidth, imageheight;
	double Xo, Yo, Zo;
	double omega, phi, kapa;
}info;

typedef struct
{
	double Rx, Ry, Rz;
	int IN;
}points;

typedef struct
{
	double Rx, Ry, Rz;
	int IN;
	int Px, Py, Pz;
}values;

/*-----------------------------------------------------------------------
|                         Function Prototypes           				|
-----------------------------------------------------------------------*/
int linesoffilex(FILE *stream);
values *loadfile(FILE *stream, int i);
values *sorttod(values real[], int i);
int threshold(values sort[], int i, info orient);
int cmp(const void* a, const void* b);
double collinearity(double xo, double c, double XA, double XO, double YA,
					double YO, double ZA, double ZO, double omega,
					double phi, double kapa);
//load.h
points *loadxyzi(FILE *estream, int i);
values *loadxyzipxpypz(FILE *stream, int i);
info orientation(FILE *ostream);
//errors.h
void error();