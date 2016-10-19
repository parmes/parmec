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
#include "parmec_ispc.h"

using namespace parmec;
using namespace ispc; /* ispc_num_cores */

int main (int argc, char *argv[])
{
  if (argc == 1)
  {
    printf ("SYNOPSIS: parmec [-threads n] path/to/file.py\n");
    printf ("         -threads n: number of threads (default: hardware supported maximum)\n");
    return 1;
  }
  else
  {
    init();

    if (strcmp (argv[1], "-threads") == 0 && argc > 2)
    {
      threads = atoi (argv[2]);
    }
    else threads = ispc_num_cores();

    int len = strlen (argv[argc-1]);
    ERRMEM (outpath = new char [len+1]);
    strcpy (outpath, argv[argc-1]);
   
    if (outpath[len-3] != '.' ||
        outpath[len-2] != 'p' ||
	outpath[len-1] != 'y')
    {
      fprintf (stderr, "ERROR: input file does not have '.py' extension!\n");
      return 2;
    }
    else outpath[len-3] = '\0';

    input (argv[argc-1]);
  }

  return 0;
}
