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

/****************************************************************************/
/****************************************************************************/
/*                                                                          */
/*                      PROTOTYPES FOR FILES IN UTIL                        */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  CC_SAFE_MALLOC(nnum,type)                                               */
/*    int nnum (the number of objects to be malloced)                       */
/*    data type (the sort of objects to be malloced)                        */
/*    RETURNS a pointer to the allocated space. If out of memory,           */
/*            it prints an error message and returns NULL.                  */
/*                                                                          */
/*  CC_FREE(object,type)                                                    */
/*    type *object (pointer to previously allocated space)                  */
/*    data type (the sort of object)                                        */
/*    ACTION: frees the memory and sets the object to NULL.                 */
/*                                                                          */
/*  CC_IFFREE(object,type)                                                  */
/*    type *object (pointer to previously allocated space)                  */
/*    data type (the sort of object)                                        */
/*    ACTION: if *object is not NULL, frees the memory and sets             */
/*            the object to NULL.                                           */
/*                                                                          */
/*  CC_PTR_ALLOC_ROUTINE (type, functionname, chunklist, freelist)          */
/*    data type (the sort of objects)                                       */
/*    string functionname (the generated function)                          */
/*    CCbigchunkptr *chunklist (used to accumulate bigchunks)               */
/*    type *freelist (used for the linked list of objects)                  */
/*    ACTION: Generates a function ("functionname") that returns            */
/*            (type *) objects, keeping the free ones on freelist           */
/*            and getting its space from calls to                           */
/*            CCutil_bigchunkalloc.                                         */
/*                                                                          */
/*  CC_PTR_FREE_ROUTINE (type, functionname, freelist)                      */
/*    Parameters as above.                                                  */
/*    ACTION: Generates a function that adds an object to the               */
/*            freelist.                                                     */
/*                                                                          */
/*  CC_PTR_FREE_LIST_ROUTINE (type, functionname, freefunction)             */
/*    Parameters defined as above, with freefunction the function           */
/*    generated by CC_PTR_FREE_ROUTINE.                                     */
/*    ACTION: Generates a function to free a linked list of                 */
/*            objects using calls to freefunction.                          */
/*                                                                          */
/*  CC_PTR_FREE_WORLD_ROUTINE (type, functionname, chunklist, freelist)     */
/*    Parameters defined as above.                                          */
/*    ACTION: Generates a function that returns all of the                  */
/*            memory used in the CC_PTR_ALLOC_ROUTINE allocations           */
/*            back to the global supply of CCbigchunkptrs.                  */
/*                                                                          */
/*  CC_PTR_LEAKS_ROUTINE (type, name, chunklist, freelist, field,           */
/*      fieldtype)                                                          */
/*    As above, with "field" the name of a "fieldtype" field in the         */
/*    object type that can be set to 0 or to 1.                             */
/*    ACTION: Generates a function that checks to see that we have          */
/*            not leaked any of the objects.                                */
/*                                                                          */
/*  CC_PTR_STATUS_ROUTINE (type, name, chunklist, freelist)                 */
/*       ACTION: Like LEAKS, but does not check for duplicates (and so      */
/*               does not corrupt the objects).                             */
/*                                                                          */
/*    NOTES:                                                                */
/*       These routines use the functions in allocrus.c.  The PTR macros    */
/*    generate the functions for allocating objects for linked lists. They  */
/*    get their raw memory from the bigchunk supply, so foo_free_world      */
/*    (generated by CC_PTR_FREE_WORLD_ROUTINE) should be called for each    */
/*    type of linked object "foo" when closing down the local memory.       */
/*       To use these functions, put the macros near the top of the file    */
/*    before any calls to the functions (since the macros also write the    */
/*    function prototypes). If you use CC_PTR_FREE_LIST_ROUTINE for foo,    */
/*    you must also use CC_PTR_FREE_ROUTINE, and                            */
/*    CC_PTR_FREE_LIST_ROUTINE must be listed after CC_PTR_FREE_ROUTINE     */
/*    (to get the prototype).                                               */
/*                                                                          */
/****************************************************************************/

#ifndef __UTIL_H
#define __UTIL_H

#include "config.h"

#define CC_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

#define CCutil_MAXDOUBLE (1e30)
#define CCutil_MAXINT    (2147483647)

#define CCcheck_rval(rval,msg) {                                          \
    if ((rval)) {                                                          \
        fprintf (stderr, "%s\n", (msg));                                   \
        goto CLEANUP;                                                      \
    }                                                                      \
}

#define CCcheck_NULL(item,msg) {                                           \
    if ((!item)) {                                                         \
        fprintf (stderr, "%s\n", (msg));                                   \
        rval = 1;                                                          \
        goto CLEANUP;                                                      \
    }                                                                      \
}


#define CC_SBUFFER_SIZE (4000)
#define CC_SFNAME_SIZE (32)

typedef struct CC_SFILE {
    int           status;
    int           desc;
    int           type;
    int           chars_in_buffer;
    int           current_buffer_char;     /* only used for reading */
    int           bits_in_last_char;       /* writing: number of empty bits in
                                            * buffer[chars_in_buffer];
                                            * reading: number of full bits in
                                            * buffer[?] */
    int           pos;
    char          fname[CC_SFNAME_SIZE];
    char          hname[CC_SFNAME_SIZE];
    unsigned char buffer[CC_SBUFFER_SIZE];
} CC_SFILE;

typedef struct CCrandstate {
    int a;
    int b;
    int arr[55];
} CCrandstate;



/****************************************************************************/
/*                                                                          */
/*                             allocrus.c                                   */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                   MEMORY ALLOCATION MACROS                               */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 24, 1995 (cofeb24)                                       */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#define CC_SAFE_MALLOC(nnum,type)                                          \
    (type *) CCutil_allocrus (((size_t) (nnum)) * sizeof (type))

#define CC_FREE(object,type) {                                             \
    CCutil_freerus ((void *) (object));                                    \
    object = (type *) NULL;                                                \
}

#define CC_IFFREE(object,type) {                                           \
    if ((object)) CC_FREE ((object),type);                                 \
}

#define CC_PTRWORLD_ALLOC_ROUTINE(type, ptr_alloc_r, ptr_bulkalloc_r)        \
                                                                             \
static int ptr_bulkalloc_r (CCptrworld *world, int nalloc)                   \
{                                                                            \
    CCbigchunkptr *bp;                                                       \
    int i;                                                                   \
    int count = CC_BIGCHUNK / sizeof ( type );                               \
    type *p;                                                                 \
                                                                             \
    while (nalloc > 0) {                                                     \
        bp = CCutil_bigchunkalloc ();                                        \
        if (bp == (CCbigchunkptr *) NULL) {                                  \
            fprintf (stderr, "ptr alloc failed\n");                          \
            return 1;                                                        \
        }                                                                    \
        bp->next = world->chunklist ;                                        \
        world->chunklist = bp;                                               \
                                                                             \
        p = ( type * ) bp->this_one;                                         \
        for (i=count-2; i>=0; i--) {                                         \
            p[i].next = &p[i+1];                                             \
        }                                                                    \
        p[count - 1].next = (type *) world->freelist;                        \
        world->freelist = (void *) p;                                        \
        nalloc -= count;                                                     \
    }                                                                        \
    return 0;                                                                \
}                                                                            \
                                                                             \
static type *ptr_alloc_r (CCptrworld *world)                                 \
{                                                                            \
    type *p;                                                                 \
                                                                             \
    if (world->freelist == (void *) NULL) {                                  \
        if (ptr_bulkalloc_r (world, 1)) {                                    \
            fprintf (stderr, "ptr alloc failed\n");                          \
            return ( type * ) NULL;                                          \
        }                                                                    \
    }                                                                        \
    p = (type *) world->freelist ;                                           \
    world->freelist = (void *) p->next;                                      \
                                                                             \
    return p;                                                                \
}

#define CC_PTRWORLD_FREE_ROUTINE(type, ptr_free_r)                           \
                                                                             \
static void ptr_free_r (CCptrworld *world, type *p)                          \
{                                                                            \
    p->next = (type *) world->freelist ;                                     \
    world->freelist = (void *) p;                                            \
}

#define CC_PTRWORLD_LISTADD_ROUTINE(type, entrytype, ptr_listadd_r, ptr_alloc_r) \
                                                                             \
static int ptr_listadd_r (type **list, entrytype x, CCptrworld *world)       \
{                                                                            \
    if (list != (type **) NULL) {                                            \
        type *p = ptr_alloc_r (world);                                       \
                                                                             \
        if (p == (type *) NULL) {                                            \
            fprintf (stderr, "ptr list add failed\n");                       \
            return 1;                                                        \
        }                                                                    \
        p->this = x;                                                         \
        p->next = *list;                                                     \
        *list = p;                                                           \
    }                                                                        \
    return 0;                                                                \
}

#define CC_PTRWORLD_LISTFREE_ROUTINE(type, ptr_listfree_r, ptr_free_r)       \
                                                                             \
static void ptr_listfree_r (CCptrworld *world, type *p)                      \
{                                                                            \
    type *next;                                                              \
                                                                             \
    while (p != (type *) NULL) {                                             \
        next = p->next;                                                      \
        ptr_free_r (world, p);                                               \
        p = next;                                                            \
    }                                                                        \
}

#define CC_PTRWORLD_LEAKS_ROUTINE(type, ptr_leaks_r, field, fieldtype)       \
                                                                             \
static int ptr_leaks_r (CCptrworld *world, int *total, int *onlist)          \
{                                                                            \
    int count = CC_BIGCHUNK / sizeof ( type );                               \
    int duplicates = 0;                                                      \
    type * p;                                                                \
    CCbigchunkptr *bp;                                                       \
                                                                             \
    *total = 0;                                                              \
    *onlist = 0;                                                             \
                                                                             \
    for (bp = world->chunklist ; bp; bp = bp->next)                          \
        (*total) += count;                                                   \
                                                                             \
    for (p = (type *) world->freelist ; p; p = p->next) {                    \
        (*onlist)++;                                                         \
        p-> field = ( fieldtype ) 0;                                         \
    }                                                                        \
    for (p = (type *) world->freelist ; p; p = p->next) {                    \
        if ((unsigned long long) p-> field == (unsigned long long) (size_t) 1)                           \
            duplicates++;                                                    \
        else                                                                 \
            p-> field = ( fieldtype ) (size_t) 1;                            \
    }                                                                        \
    if (duplicates) {                                                        \
        fprintf (stderr, "WARNING: %d duplicates on ptr free list \n",       \
                 duplicates);                                                \
    }                                                                        \
    return *total - *onlist;                                                 \
}

#define CC_PTRWORLD_ROUTINES(type, ptr_alloc_r, ptr_bulkalloc_r, ptr_free_r) \
CC_PTRWORLD_ALLOC_ROUTINE (type, ptr_alloc_r, ptr_bulkalloc_r)               \
CC_PTRWORLD_FREE_ROUTINE (type, ptr_free_r)

#define CC_PTRWORLD_LIST_ROUTINES(type, entrytype, ptr_alloc_r, ptr_bulkalloc_r, ptr_free_r, ptr_listadd_r, ptr_listfree_r) \
CC_PTRWORLD_ROUTINES (type, ptr_alloc_r, ptr_bulkalloc_r, ptr_free_r)        \
CC_PTRWORLD_LISTADD_ROUTINE (type, entrytype, ptr_listadd_r, ptr_alloc_r)    \
CC_PTRWORLD_LISTFREE_ROUTINE (type, ptr_listfree_r, ptr_free_r)

#define CC_BIGCHUNK ((int) ((1<<16) - sizeof (CCbigchunkptr) - 16))

struct CCbigchunk;

typedef struct CCbigchunkptr {
    void                 *this_one;
    struct CCbigchunk    *this_chunk;
    struct CCbigchunkptr *next;
} CCbigchunkptr;


typedef struct CCptrworld {
    int refcount;
    void *freelist;
    CCbigchunkptr *chunklist;
} CCptrworld;



void
   *CCutil_allocrus (size_t size),
   *CCutil_reallocrus (void *ptr, size_t size),
    CCutil_freerus (void *p),
    CCutil_bigchunkfree (CCbigchunkptr *bp),
    CCptrworld_init (CCptrworld *world),
    CCptrworld_add (CCptrworld *world),
    CCptrworld_delete (CCptrworld *world);

int
    CCutil_reallocrus_scale (void **pptr, int *pnnum, int count, double scale,
        size_t size),
    CCutil_reallocrus_count (void **pptr, int count, size_t size);

CCbigchunkptr
    *CCutil_bigchunkalloc (void);



/****************************************************************************/
/*                                                                          */
/*                             dheaps_i.c                                   */
/*                                                                          */
/****************************************************************************/

typedef struct CCdheap {
    double  *key;
    int     *entry;
    int     *loc;
    int     total_space;
    int     size;
} CCdheap;


void
    CCutil_dheap_free (CCdheap *h),
    CCutil_dheap_delete (CCdheap *h, int i),
    CCutil_dheap_changekey (CCdheap *h, int i, double newkey);

int
    CCutil_dheap_init (CCdheap *h, int k),
    CCutil_dheap_resize (CCdheap *h, int newsize),
    CCutil_dheap_findmin (CCdheap *h),
    CCutil_dheap_deletemin (CCdheap *h),
    CCutil_dheap_insert (CCdheap *h, int i);




/****************************************************************************/
/*                                                                          */
/*                             edgelen.c                                    */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*  Before defining CCUTIL_EDGELEN_FUNCTIONPTR, read the notes at the top   */
/*  of edgelen.c, and carefully consider the consequences.  You probably    */
/*  do not want CCUTIL_EDGELEN_FUNCTIONPTR defined.                         */
/*                                                                          */
/****************************************************************************/

#undef  CCUTIL_EDGELEN_FUNCTIONPTR

typedef struct CCdata_user {
    double  *x;
    double  *y;
} CCdata_user;

typedef struct CCdata_rhvector {
    int dist_00;
    int dist_01;
    int dist_02;
    int dist_12;
    int dist_22;
    double p;
    int rhlength;
    char *space;
    char **vectors;
} CCdata_rhvector;

typedef struct CCdatagroup {
    int    (*edgelen) (int i, int j, struct CCdatagroup *dat);
    double  *x;
    double  *y;
    double  *z;
    int    **adj;
    int     *adjspace;
    int    **len;
    int     *lenspace;
    int     *degree;
    int      norm;
    int      dsjrand_param;
    int      default_len;     /* for edges not in sparse graph   */
    int      sparse_ecount;   /* number of edges in sparse graph */
    double   gridsize;        /* for toroidal norm */
    double   dsjrand_factor;
    CCdata_rhvector rhdat;
    CCdata_user     userdat;
    int      ndepot;          /* used with the subdivision code   */
    int      orig_ncount;     /* just ncount-ndepot               */
    int     *depotcost;       /* cost from each node to the depot */
    int     *orig_names;      /* the nodes names from full problem */
} CCdatagroup;



#ifdef CCUTIL_EDGELEN_FUNCTIONPTR
extern int
  (*CCutil_dat_edgelen) (int i, int j, CCdatagroup *dat);
#else  /* CCUTIL_EDGELEN_FUNCTIONPTR */
int
    CCutil_dat_edgelen (int i, int j, CCdatagroup *dat);
#endif /* CCUTIL_EDGELEN_FUNCTIONPTR */

int
    CCutil_dat_setnorm (CCdatagroup *dat, int norm);

void
    CCutil_dat_getnorm (CCdatagroup *dat, int *norm),
    CCutil_dsjrand_init (CCdatagroup *dat, int maxdist, int seed),
    CCutil_init_datagroup (CCdatagroup *dat),
    CCutil_freedatagroup (CCdatagroup *dat);


#define CC_KD_NORM_TYPE    128            /* Kdtrees work      */
#define CC_X_NORM_TYPE     256            /* Old nearest works */
#define CC_JUNK_NORM_TYPE  512            /* Nothing works     */

#define CC_D2_NORM_SIZE      1024         /* x,y coordinates   */
#define CC_D3_NORM_SIZE      2048         /* x,y,z coordinates */
#define CC_MATRIX_NORM_SIZE  4096         /* adj matrix        */

#define CC_NORM_BITS      (CC_KD_NORM_TYPE | CC_X_NORM_TYPE | \
                           CC_JUNK_NORM_TYPE)
#define CC_NORM_SIZE_BITS (CC_D2_NORM_SIZE | CC_D3_NORM_SIZE | \
                           CC_MATRIX_NORM_SIZE)

#define CC_MAXNORM        (0 |   CC_KD_NORM_TYPE |     CC_D2_NORM_SIZE)
#define CC_EUCLIDEAN_CEIL (1 |   CC_KD_NORM_TYPE |     CC_D2_NORM_SIZE)
#define CC_EUCLIDEAN      (2 |   CC_KD_NORM_TYPE |     CC_D2_NORM_SIZE)
#define CC_EUCLIDEAN_3D   (3 |    CC_X_NORM_TYPE |     CC_D3_NORM_SIZE)
#define CC_USER           (4 | CC_JUNK_NORM_TYPE |                   0)
#define CC_ATT            (5 |    CC_X_NORM_TYPE |     CC_D2_NORM_SIZE)
#define CC_GEOGRAPHIC     (6 |    CC_X_NORM_TYPE |     CC_D2_NORM_SIZE)
#define CC_MATRIXNORM     (7 | CC_JUNK_NORM_TYPE | CC_MATRIX_NORM_SIZE)
#define CC_DSJRANDNORM    (8 | CC_JUNK_NORM_TYPE |                   0)
#define CC_CRYSTAL        (9 |    CC_X_NORM_TYPE |     CC_D3_NORM_SIZE)
#define CC_SPARSE        (10 | CC_JUNK_NORM_TYPE |                   0)
#define CC_RHMAP1        (11 | CC_JUNK_NORM_TYPE |                   0)
#define CC_RHMAP2        (12 | CC_JUNK_NORM_TYPE |                   0)
#define CC_RHMAP3        (13 | CC_JUNK_NORM_TYPE |                   0)
#define CC_RHMAP4        (14 | CC_JUNK_NORM_TYPE |                   0)
#define CC_RHMAP5        (15 | CC_JUNK_NORM_TYPE |                   0)
#define CC_EUCTOROIDAL   (16 | CC_JUNK_NORM_TYPE |     CC_D2_NORM_SIZE)
#define CC_GEOM          (17 |    CC_X_NORM_TYPE |     CC_D2_NORM_SIZE)
#define CC_MANNORM       (18 |   CC_KD_NORM_TYPE |     CC_D2_NORM_SIZE)
#define CC_SUBDIVISION   (99 | CC_JUNK_NORM_TYPE |                   0)

/* Distances CC_RHMAP1 through CC_RHMAP5 are for an application to          */
/* radiation hybrid mapping in genetics, explained in: Agarwala R,          */
/* Applegate DL,  Maglott D, Schuler GD, Schaffer AA: A Fast and Scalable   */
/* Radiation Hybrid Map Construction and Integration Strategy. Genome       */
/* Research, 10:350-364, 2000.  The correspondence to the distance function */
/* terms used in that paper is: CC_RMAP1 (weighted_ocb), CC_RHMAP2          */
/* (normalized_mle), CC_RHMAP3 (base_mle), CC_RHMAP4 (extended_mle),        */
/* CC_RHMAP5 (normalized_ocb)                                               */

/* For X-NORMS, scales are such that |x[i] - x[j]| * scale <= edgelen(i,j). */
/* Geographic is slightly off, since the fractional part of x[i] is really  */
/* really minutes, not fractional degrees.                                  */




/****************************************************************************/
/*                                                                          */
/*                             util.c                                    */
/*                                                                          */
/****************************************************************************/


void
    CCutil_int_perm_quicksort (int *perm, int *len, int n),
    CCutil_double_perm_quicksort (int *perm, double *len, int n),
    CCutil_rselect (int *arr, int l, int r, int m, double *coord,
        CCrandstate *rstate);

/* since urandom's generator does everything modulo CC_PRANDMAX, if two
 * seeds are congruent mod x and x|CC_PRANDMAX, then the resulting numbers
 * will be congruent mod x.  One example was if CC_PRANDMAX = 1000000000 and
 * urandom is used to generate a point set from a 1000x1000 grid, seeds
 * congruent mod 1000 generate the same point set.
 *
 * For this reason, we use 1000000007 (a prime)
 */
#define CC_PRANDMAX 1000000007

void
   CCutil_sprand (int seed, CCrandstate *r);

int
   CCutil_lprand (CCrandstate *r);

double CCutil_zeit (void);

int
    CCutil_bix_getopt (int argc, char **argv, const char *def, int *p_optind,
        char **p_optarg);


#define CC_BIX_GETOPT_UNKNOWN -3038


int
    CCutil_edge_to_cycle (int ncount, int *elist, int *yesno, int *cyc);

int
    CCutil_gettsplib (char *datname, int *ncount, CCdatagroup *dat);


#endif /* __UTIL_H */
