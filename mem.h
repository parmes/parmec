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

#ifndef __mem__
#define __mem__

typedef struct memory_pool MEM;

struct memory_pool 
{
  void *blocks; /* list of allocated memory blocks */
  char *freechunk; /* next free chunk of memory */
  char *lastchunk; /* last chunk in current block */
  void *deadchunks; /* list of dealocated chunks of memory */
  size_t chunksize; /* size of a chunk */
  size_t chunksinblock; /* number of memory chunks in a block */
};

/* allocate global zero'd memory */
void* MEM_CALLOC (size_t size);

/* initialize memory pool */
void MEM_Init (MEM *pool, size_t chunksize, size_t chunksinblock);

/* allocate a chunk of memory from the pool */
void* MEM_Alloc (MEM *pool);

/* free a chunk of memory to the pool */
void MEM_Free (MEM *pool, void *chunk);

/* return amount of memory in the pool */
size_t MEM_Size (MEM *pool);

/* release memory pool memory back to system */
void MEM_Release (MEM *pool);

#endif
