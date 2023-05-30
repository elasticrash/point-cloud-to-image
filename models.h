typedef double T;

typedef struct {
  char name[6];
  double c, xo, yo;
  double k1, k2;
  int imagewidth, imageheight;
  double Xo, Yo, Zo;
  double omega, phi, kapa;
} info;

typedef struct {
  double Rx, Ry, Rz;
  int IN;
} points;

typedef struct {
  double Rx, Ry, Rz;
  int IN;
  int Px, Py, Pz;
} values;
