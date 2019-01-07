/*
The MIT License (MIT)

Copyright (c) 2016 Tomasz Koziara

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "macros.h"
#include "parmec.h"
#include "output.h"
#include "parmec_ispc.h"
#include "version.h"

using namespace parmec;
using namespace ispc; /* ispc_num_cores */

int main (int argc, char *argv[])
{
  if (argc == 1)
  {
    printf ("VERSION: %s %s\n", VERSION_DATE, VERSION_HASH);
    printf ("SYNOPSIS: parmec [-ntasks n] path/to/file.py\n");
    printf ("         -ntasks n: number of tasks (default: hardware supported maximum)\n");
    return 1;
  }
  else
  {
    init();

    if (strcmp (argv[1], "-ntasks") == 0 && argc > 2)
    {
      ntasks = atoi (argv[2]);
    }

    for (int i = 1; i < argc; i ++)
    {
      int len = strlen(argv[i]);
      if (len > 3 && strcmp(argv[i]+len-3, ".py") == 0)
      {
        input (argv[i], argv, argc);
	break;
      }
    }

    output_reset ();
  }

  return 0;
}
