/* ===========================================================================
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

#include <kapp/args.h>
#include <kapp/main.h>
#include <kfs/file.h>
#include <klib/log.h>
#include <klib/rc.h>
#include <klib/vector.h>
#include <kproc/queue.h>
#include <kproc/thread.h>
#include <kproc/timeout.h>
//#include <memview.hpp>
//#include <samextract.hpp>
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


extern "C" {
extern int SAMparse();
int moredata(char * buf,int * numbytes, int maxbytes)
{
    fprintf(stderr,"moredata %d %d\n", *numbytes, maxbytes);
    *numbytes=maxbytes;
    int bytesremain=mmapbuf_sz-(mmapbuf_cur-mmapbuf);
    if (maxbytes > bytesremain) *numbytes=bytesremain;

    memmove(buf,mmapbuf_cur,*numbytes);
    mmapbuf_cur+=*numbytes;

    return 0;
}

void SAMerror(const char * s)
{
    fprintf(stderr,"ERR: %s\n", s);
}
} // extern C

rc_t inflater(const KThread * self, void * in)
{
    fprintf(stderr,"\tThread started.\n");
    KQueue *que=(KQueue *)in;
    struct timeout_t tm;
    TimeoutInit(&tm, 100000); // 1 second

    while (1)
    {
        void * where=0;
        fprintf(stderr,"\tchecking queue...");
        rc_t rc=KQueuePop(que, &where, &tm);
        if (rc==0)
        {
            fprintf(stderr,"\t\tGot work: %p\n", where);
            size_t f=(long)where;
            for (int i=0; i!=10000000; ++i) f*=32131;
            fprintf(stderr,"result=%lld\n",f);
        } else
        if ((int)GetRCObject(rc)== rcTimeout) 
        {
            fprintf(stderr,"queue empty\n");
            return 0;
        } else
        if ((int)GetRCObject(rc)== rcData) 
        {
            fprintf(stderr,"queue data\n");
            return 0;
        } else
        {
            fprintf(stderr,"rc=%d\n",rc);
        }
    }

    return 0;
}

int threadinflate(void)
{
    z_stream strm;

    KQueue * que;

    const int numthreads=16;

    rc_t rc;

    rc=KQueueMake(&que, numthreads);
    if (rc) return 1;

    KThread * threads[16];
    for (auto i=0; i!=numthreads; ++i)
        rc=KThreadMake(&threads[i], inflater, que);


    struct timeout_t tm;
    TimeoutInit(&tm, 1000); // 1 second

    mmapbuf_cur=mmapbuf;

    while(1)
    {
        memset(&strm,0,sizeof strm);
        strm.next_in=(Bytef*)mmapbuf_cur;
        if (mmapbuf_sz > INT32_MAX)
            strm.avail_in=INT32_MAX;
        else
            strm.avail_in=(uInt)mmapbuf_sz;
        int zrc=inflateInit2(&strm, MAX_WBITS + 16); // Only gzip format
        switch (zrc)
        {
            case Z_OK:
                break;
            case Z_MEM_ERROR:
                fprintf(stderr,"error: Out of memory in zlib\n");
                break;
            case Z_VERSION_ERROR:
                fprintf(stderr,"zlib version is not compatible; need version %s but have %s\n", ZLIB_VERSION,zlibVersion());
                break;
            case Z_STREAM_ERROR:
                fprintf(stderr,"zlib stream error\n");
                break;
            default:
                fprintf(stderr,"zlib error\n");
                break;
        }

        gz_header head;
        uint8_t extra[256];
        memset(&head,0,sizeof head);
        head.extra=extra;
        head.extra_max=sizeof(extra);
        char outbuf[64];
        strm.next_out=(Bytef*)outbuf;
        strm.avail_out=64;
        zrc=inflateGetHeader(&strm,&head);
        while (head.done==0)
        {
            zrc=inflate(&strm,Z_BLOCK);
            if (zrc!=Z_OK)
            {
                fprintf(stderr,"inflate error\n");
            }
        }
        fprintf(stderr,"got header\n");
        if (head.extra && head.extra_len==6 && head.extra[0]=='B' && head.extra[1]=='C' && head.extra[2]==2 && head.extra[3]==0)
        {
            uint16_t bsize=head.extra[4]+head.extra[5]*256;
            inflateEnd(&strm);
            fprintf(stderr,"extra:   %d\n",bsize);
            fprintf(stderr,"total_in:%d\n",strm.avail_in);
            fprintf(stderr,"buf   in:%p\n",mmapbuf_cur);
            fprintf(stderr,"next_in: %p\n",strm.next_in);
            fprintf(stderr,"offset:  %ld\n",mmapbuf_cur-mmapbuf);
            while(1)
            {
                rc=KQueuePush(que, strm.next_in, &tm);
                fprintf(stderr,"queued:%d %d\n",rc, rcTimeout);
                if ((int)GetRCObject(rc)== rcTimeout) fprintf(stderr,"queue full\n");
                if (rc == 0) break;
            }
            if (bsize<=28) return 0;
            mmapbuf_cur=mmapbuf_cur+bsize+1;
        } else
        {
            fprintf(stderr,"error: BAM required extra extension not found\n");
        }
    }

    KQueueSeal(que);
    // TODO: Wait and Release threads
    KQueueRelease(que);

    return 0;
}


static
void samload(char const path[])
{
    tags=strdup("");
    seqnames=strdup("");
    ids=strdup("");

    int fd=open(path, O_RDONLY);
    if (fd==-1)
    {
        fprintf(stderr,"Cannot open %s:%s\n", path, strerror(errno));
    }

    struct stat st;
    if (fstat(fd,&st))
    {
        fprintf(stderr,"Cannot stat %s:%s\n", path, strerror(errno));
    }
    mmapbuf_sz=st.st_size;

    mmapbuf=(char *)mmap(NULL, (size_t)mmapbuf_sz, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd); // TODO: Still need to munmap
//    madvise(mmapbuf, mmapbuf_sz, MADV_SEQUENTIAL | MADV_WILLNEED);
    madvise(mmapbuf, mmapbuf_sz, MADV_SEQUENTIAL);

    if (mmapbuf==MAP_FAILED)
    {
        fprintf(stderr,"Cannot mmap %s:%s\n", path, strerror(errno));
        return;
    }

    if (!memcmp(mmapbuf,"\x1f\x8b\x08",3))
    {
        fprintf(stderr,"gzip file\n");
        threadinflate();
        return;
    } else if (mmapbuf[0]=='@')
    {
        fprintf(stderr,"SAM file\n");
        mmapbuf_cur=mmapbuf;
        SAMparse();
        return;
    }
}

rc_t CC UsageSummary(char const *name)
{
    fprintf(stderr,"Usage: %s file.{sb}am [file...]\n", name);
    return 0;
}

rc_t CC Usage(Args const *args)
{
    return 0;
}

rc_t CC KMain(int argc, char *argv[])
{
    if (argc == 1) {
        UsageSummary(argv[0]);
        return 0;
    }
    while (--argc) {
        samload(*++argv);
    }
    return 0;
}

