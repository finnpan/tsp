set(
    source_files
    allocrus.c bgetopt.c  dheaps_i.c edgelen.c  edgemap.c  edgeutil.c
    eunion.c   fastread.c getdata.c  safe_io.c  sortrus.c  subdiv.c
    urandom.c  util.c     zeit.c     genhash.c  priority.c signal.c
)

set(library_name util)

add_library(${library_name} STATIC ${source_files})
