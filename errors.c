#include <stdio.h>
#include <stdlib.h>

#include "errors.h"
/*-----------------------------------------------------------------------
|					       allocate memory error
|
-----------------------------------------------------------------------*/
void error() { printf("not enough memory\n"); }
/*-----------------------------------------------------------------------
|					           file open error
|
-----------------------------------------------------------------------*/
void errorf(FILE *stream) {
  if (!stream) {
    printf("file failed to open");
  }
}
