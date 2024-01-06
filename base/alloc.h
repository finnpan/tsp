/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

#ifndef __ALLOC_H
#define __ALLOC_H

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CC_SAFE_MALLOC(nnum,type)                                          \
    (type *) CCutil_allocrus (((size_t) (nnum)) * sizeof (type))

#define CC_FREE(object,type) {                                             \
    CCutil_freerus ((void *) (object));                                    \
    object = (type *) NULL;                                                \
}

#define CC_IFFREE(object,type) {                                           \
    if ((object)) CC_FREE ((object),type);                                 \
}


void
   *CCutil_allocrus (size_t size),
    CCutil_freerus (void *p);


#ifdef __cplusplus
}
#endif

#endif /* __ALLOC_H */
