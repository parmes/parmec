/*
The MIT License (MIT)

Copyright (c) 2015 Tomasz Koziara

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

#include "mem.h"

#ifndef __map__
#define __map__

typedef struct map MAP; /* map type */
typedef int (*MAP_Compare) (void*, void*); /* comparison callback */

struct map
{
  MAP *p, *l, *r; /* parent, left, right */
  char colour; /* red, blue */
  void *key, *data; /* user data */
};

/* insert an element into the map (return NULL if the key is already mapped) */
MAP* MAP_Insert (MEM *pool, MAP **root,
  void *key, void *data, MAP_Compare compare); /* 'compare' == NULL 
						  => direct pointer 'key' comparison 
						  (this applies below as well) */
/* find data for a specific key */
void* MAP_Find (MAP *root, void *key, MAP_Compare compare);

/* find map node for a specific key */
MAP* MAP_Find_Node (MAP *root, void *key, MAP_Compare compare);

/* delete an element from the map */
void* MAP_Delete (MEM *pool, MAP **root, void *key, MAP_Compare compare);

/* delete a specific map node => return the next node by data order */
MAP* MAP_Delete_Node (MEM *pool, MAP **root, MAP *node);

/* postorder traverse and free map memory */
void MAP_Free (MEM *pool, MAP **root);

/* return number of items */
int MAP_Size (MAP *root);

/* first element */
MAP* MAP_First (MAP *root);

/* last element */
MAP* MAP_Last (MAP *root);

/* previous element */
MAP* MAP_Prev (MAP *node);

/* next element */
MAP* MAP_Next (MAP *node);

#endif
