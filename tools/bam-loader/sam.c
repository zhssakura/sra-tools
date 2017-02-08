/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <regex.h>
#include <stdint.h>
#include <unistd.h>
#include <zlib.h>

char *tags=NULL; // Space delimited tags seen in current line
char *seqnames=NULL;
char *ids=NULL;

char * mmapbuf=NULL;
off_t  mmapbuf_sz=0;
char * mmapbuf_cur=NULL;

extern int SAMparse();

int moredata(char * buf,int * numbytes, int maxbytes)
{
    fprintf(stderr,"moredata %d %d\n", *numbytes, maxbytes);
    *numbytes=maxbytes;
    size_t bytesremain=mmapbuf_sz-(mmapbuf_cur-mmapbuf);
    if (maxbytes > bytesremain) *numbytes=bytesremain;

    memcpy(buf,mmapbuf_cur,*numbytes);
    mmapbuf_cur+=*numbytes;

    return 0;
}

void SAMerror(const char * s) 
{
    fprintf(stderr,"ERR: %s\n", s);
}

int main(int argc, char * argv[])
{
    tags=strdup("");
    seqnames=strdup("");
    ids=strdup("");

    if (argc!=2)
    {
        fprintf(stderr,"Usage: %s filename.{sb}am\n", argv[0]);
        return 1;
    }

    const char *fname=argv[1];
    int fd=open(fname, O_RDONLY);
    if (fd==-1)
    {
        fprintf(stderr,"Cannot open %s:%s\n", fname, strerror(errno));
        return 1;
    }

    struct stat st;
    if (fstat(fd,&st))
    {
        fprintf(stderr,"Cannot stat %s:%s\n", fname, strerror(errno));
        return 1;
    }
    mmapbuf_sz=st.st_size;

    mmapbuf=mmap(NULL, (size_t)mmapbuf_sz, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd); // TODO: Still need to munmap
//    madvise(mmapbuf, mmapbuf_sz, MADV_SEQUENTIAL | MADV_WILLNEED);
    madvise(mmapbuf, mmapbuf_sz, MADV_SEQUENTIAL);

    if (mmapbuf==MAP_FAILED)
    {
        fprintf(stderr,"Cannot mmap %s:%s\n", fname, strerror(errno));
        return 1;
    }


    if (!memcmp(mmapbuf,"\x1f\x8b\x08",3))
    {
        fprintf(stderr,"gzip file\n");
        z_stream strm;

    } else if (mmapbuf[0]=='@')
    {
        fprintf(stderr,"SAM file\n");
        mmapbuf_cur=mmapbuf;
        return SAMparse();
    }

}

