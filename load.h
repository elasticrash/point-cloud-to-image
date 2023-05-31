#ifndef LOAD_H
#define LOAD_H
points *loadxyzi(FILE *estream, int i);
info orientation(FILE *ostream);
values *loadt(FILE *istream, int i);
#endif
