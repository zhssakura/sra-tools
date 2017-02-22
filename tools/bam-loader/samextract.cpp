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
#include <klib/defs.h>
#include <klib/vector.h>
#include <kproc/queue.h>
#include <kproc/thread.hpp>
#include <kproc/timeout.h>
//#include <memview.hpp>
//#include <samextract.hpp>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <regex.h>
#include <stdint.h>
#include <unistd.h>
#include <zlib.h>


// TODO: Remove globals, make into state struct
char *tags=NULL; // Space delimited tags seen in current line
char *seqnames=NULL;
char *ids=NULL;

char * mmapbuf=NULL;
off_t  mmapbuf_sz=0;
char * mmapbuf_cur=NULL;


enum chunk_state
{
    empty, compressed, uncompressed
};

typedef struct chunk_s
{
    Bytef * in=NULL;
    uInt insize=0;
    uInt outsize=0;
    chunk_state state=empty;
    Bytef out[66000];
} chunk;

typedef struct bamalign
{
    int32_t block_size;
    int32_t refID;
    int32_t pos;
    uint32_t bin_mq_nl;
    uint32_t flag_nc;
    int32_t l_seq;
    int32_t next_refID;
    int32_t next_pos;
    int32_t tlen;
} bamalign;

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


rc_t inflater(const KThread * kt, void * in)
{
    KQueue *inflatequeue=(KQueue *)in;
    struct timeout_t tm;

    z_stream strm;
    pthread_t threadid=pthread_self();
    fprintf(stderr,"\tThread %lu started.\n",threadid);

    while (1)
    {
        void * where=NULL;
        fprintf(stderr,"\t\tthread %lu checking queue\n",threadid);
        TimeoutInit(&tm, 5000); // 5 seconds
        rc_t rc=KQueuePop(inflatequeue, &where, &tm);
        if (rc==0)
        {
            chunk * c=(chunk *)where;
            fprintf(stderr,"\t\tthread %lu chunk %p size %u\n", threadid, c->in, c->insize);
            if (c->state!=compressed)
            {
                fprintf(stderr,"Inflater bad state\n");
                return 1;
            }

            memset(&strm,0,sizeof strm);
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
                    fprintf(stderr,"zlib error %s\n",strm.msg);
                    break;
            }
            strm.next_in=c->in;
            strm.avail_in=c->insize;
            strm.next_out=c->out;
            strm.avail_out=c->outsize;

            zrc=inflate(&strm,Z_NO_FLUSH);
            switch (zrc)
            {
                case Z_OK:
                    fprintf(stderr,"\t\tthread %lu OK %d %d %lu\n", threadid, strm.avail_in, strm.avail_out,strm.total_out);
                    c->outsize=strm.total_out;
                    c->state=uncompressed;
                    break;
                case Z_MEM_ERROR:
                    fprintf(stderr,"error: Out of memory in zlib\n");
                    break;
                case Z_VERSION_ERROR:
                    fprintf(stderr,"zlib version is not compatible; need version %s but have %s\n", ZLIB_VERSION,zlibVersion());
                    break;
                case Z_STREAM_ERROR:
                    fprintf(stderr,"zlib stream error %s\n",strm.msg);
                    break;
                case Z_STREAM_END:
                    fprintf(stderr,"zlib stream end%s\n",strm.msg);
                    break;
                default:
                    fprintf(stderr,"inflate error %d %s\n",zrc, strm.msg);
                    break;
            }
            inflateEnd(&strm);
        } else if ((int)GetRCObject(rc)== rcTimeout) 
        {
            fprintf(stderr,"\t\tthread %lu queue empty\n",threadid);
            if (KQueueSealed(inflatequeue))
            {
                fprintf(stderr,"\t\tQueue sealed, thread %lu complete\n", threadid);
                return 0;
            }
//            usleep(1);
        } else if ((int)GetRCObject(rc)== rcData) 
        {
            fprintf(stderr,"\t\tthread %lu queue data (empty?)\n",threadid);
            break;
        } else
        {
            fprintf(stderr,"rc=%d\n",rc);
        }
    }

    fprintf(stderr,"\t\tthread %lu terminating.\n",threadid);
    return 0;
}

class chunkview
{
  public:
    chunkview(KQueue * queue) : que(queue), c(NULL), cur(NULL) {};

    ~chunkview()
    {
        que=NULL;
        c=NULL;
        cur=NULL;
    }

    chunkview(const chunkview &) = delete;
    chunkview(const chunkview &&) = delete;
    chunkview & operator=(const chunkview &) = delete;
    chunkview & operator=(const chunkview &&) = delete;

  private:
    bool getnextchunk()
    {
        void * where=NULL;
        struct timeout_t tm;
        fprintf(stderr,"\t\tBlocker thread checking blocker queue\n");
        if (c && c->state!=uncompressed)
        {
            fprintf(stderr,"\t\tblocker bad state\n");
            return false;
        }
        if (c) c->state=empty;
        c=NULL;

        while (1)
        {
            TimeoutInit(&tm, 5000); // 5 seconds
            rc_t rc=KQueuePop(que, &where, &tm);
            if (rc==0)
            {
                c=(chunk *)where;
                if (c->state==empty)
                {
                    fprintf(stderr,"\t\tBlocker bad state\n");
                    return false;
                }
                while (c->state!=uncompressed)
                {
                    fprintf(stderr,"\t\tBlocker busy\n");
                    usleep(1);
                }
                fprintf(stderr,"\t\tBlocker thread chunk %p size %u\n", c->in, c->insize);
                break;
            } else if ((int)GetRCObject(rc)== rcTimeout) 
            {
                fprintf(stderr,"\t\tBlocker thread queue empty\n");
                if (KQueueSealed(que))
                {
                    fprintf(stderr,"\t\tQueue sealed, Blocker thread complete\n");
                    return false;
                }
    //            usleep(1);
            } else if ((int)GetRCObject(rc)== rcData) 
            {
                fprintf(stderr,"\t\tBlocker thread queue data (empty?)\n");
                return false;
            } else
            {
                fprintf(stderr,"blocker rc=%d\n",rc);
                return false;
            }
        }
        return true;
    }
  public:
    bool getbytes(void * dest, uInt len)
    {
        char * where=(char *)dest;
        while (len)
        {
            if (c==NULL || c->outsize==0)
            {
                if (!getnextchunk())
                    return false;
                cur=c->out;
            }

            size_t howmany=c->outsize;
            if (len < c->outsize) howmany=len;

            memmove(where,cur,howmany);
            len-=howmany;
            c->outsize-=howmany;
            where+=howmany;
            cur+=howmany;
        }
        return true;
    }

  private:
    KQueue * que;
    chunk * c;
    Bytef * cur;
};


rc_t blocker(const KThread * kt, void * in)
{
    KQueue *blockqueue=(KQueue *)in;

    pthread_t threadid=pthread_self();
    fprintf(stderr,"\tBlocker thread %lu started.\n",threadid);

    chunkview cv(blockqueue);

    char magic[4];
    if (!cv.getbytes(magic,4)) return 1;
    if (memcmp(magic,"BAM\x01",4))
    {
        fprintf(stderr,"BAM magic not found\n");
        return 2;
    }
    int32_t l_text;
    if (!cv.getbytes(&l_text,4)) return 1;
    if (l_text<0)
    {
        fprintf(stderr,"error: invalid l_text\n");
        return 1;
    }

    char *text=(char *)calloc(1,l_text+1);
    if (!cv.getbytes(text,l_text)) return 1;
    text[l_text]='\0';

    fprintf(stderr,"SAM header:%s\n",text);

    int32_t n_ref;
    if (!cv.getbytes(&n_ref,4)) return 1;
    if (n_ref<0)
    {
        fprintf(stderr,"error: invalid n_ref\n");
        return 1;
    }
    fprintf(stderr,"# references %d\n", n_ref);

    for (int i=0; i!=n_ref; ++i)
    {
        int32_t l_name;
        if (!cv.getbytes(&l_name,4)) return 1;
        fprintf(stderr,"%d: l_name=%d\n",i, l_name);
        if (l_name < 0)
        {
            fprintf(stderr,"error: invalid reference name length\n");
            return 1;
        }
        if (l_name > 256)
        {
            fprintf(stderr,"warning: Long reference name\n");
            return 1;
        }
        char *name=(char *)calloc(1,l_name+1);
        if (!cv.getbytes(name,l_name)) return 1;
        fprintf(stderr,"%d: reference name %s\n",i, name);
        free(name); // TODO, persist?
        int32_t l_ref;
        if (!cv.getbytes(&l_ref,4)) return 1;
        fprintf(stderr,"length of reference sequence %d=%d\n",i,l_ref);
    }

    while (1)
    {
        bamalign align;
        if (!cv.getbytes(&align,sizeof(align))) return 1;
        fprintf(stderr,"\n\n\nalignment block_size=%d refID=%d pos=%d\n", 
                align.block_size, 
                align.refID, 
                align.pos);

        if (align.block_size < 0)
        {
            fprintf(stderr,"error: invalid block_size\n");
            return 1;
        }
        if (align.pos < 0)
        {
            fprintf(stderr,"error: invalid pos\n");
            return 1;
        }
        if (align.refID < -1 || align.refID > n_ref)
        {
            fprintf(stderr,"error: bad refID\n");
            return 1;
        }
        if (align.next_refID < -1 || align.next_refID > n_ref)
        {
            fprintf(stderr,"error: bad next_refID\n");
            return 1;
        }
        fprintf(stderr,"align.bin_mq_nl=%d\n",align.bin_mq_nl);
        uint16_t bin=align.bin_mq_nl >> 16;
        uint8_t mapq=(align.bin_mq_nl >> 8) & 0xff;
        uint8_t l_read_name=align.bin_mq_nl & 0xff;
        fprintf(stderr,"bin=%d mapq=%d l_read_name=%d\n", bin, mapq, l_read_name);

        uint16_t flag=align.flag_nc >> 16;
        uint16_t n_cigar_op=align.flag_nc & 0xffff;
        fprintf(stderr,"flag=%x n_cigar_op=%d\n", flag, n_cigar_op);

        char * read_name=(char *)calloc(1,l_read_name);
        if (!cv.getbytes(read_name,l_read_name)) return 1;
        fprintf(stderr,"read_name=%s\n",read_name);

        static const char opmap[]="MIDNSHP=X???????";
        uint32_t * cigar=(uint32_t *)calloc(n_cigar_op,sizeof(uint32_t));
        if (!cv.getbytes(cigar,n_cigar_op*4)) return 1;
        for (int i=0; i!=n_cigar_op; ++i)
        {
            int32_t oplen=cigar[i] >> 4;
            int32_t op=cigar[i] & 0xf;
            fprintf(stderr,"\tcigar %d=%x len=%d %d(%c)\n", i, cigar[i], oplen, op, opmap[op]);
        }

        static const char seqmap[]="=ACMGRSVTWYHKDBN";
        char * seq=(char *)calloc(1,align.l_seq+1);
        int bytesofseq=(align.l_seq+1)/2;
        uint8_t * seqbytes=(uint8_t *)calloc(1,bytesofseq);
        if (!cv.getbytes(seqbytes,bytesofseq)) return 1;
        int i=0;
        int j=0;
        while (i < align.l_seq)
        {
            seq[i]=seqmap[seqbytes[j] >> 4];
            i+=2;
            j+=1;
        }
        i=1;
        j=0;
        while (i < align.l_seq)
        {
            seq[i]=seqmap[seqbytes[j] & 0xf];
            i+=2;
            j+=1;
        }
        free(seqbytes);

        char *qual=(char*)calloc(1,align.l_seq);
        if (!cv.getbytes(qual,align.l_seq)) return 1;

        fprintf(stderr,"%d sequences\n",align.l_seq);
        for (int i=0; i!=align.l_seq; ++i)
            fprintf(stderr," seq#%d %c %.2x ",i, seq[i], qual[i]);
        fprintf(stderr,"\nseq=%s\n",seq);

        int remain=align.block_size-(sizeof(align)+l_read_name+n_cigar_op*4+bytesofseq+align.l_seq)+4; // TODO
        fprintf(stderr,"%d bytes remaining for ttvs\n",remain);
        char * ttvs=(char*)calloc(1,remain); // TODO: alloca? check <64K guarantee
        if (!cv.getbytes(ttvs,remain)) return 1;
        char * cur=ttvs;
        while (cur<ttvs+remain)
        {
            char tag[2];
            char c;
            int8_t i8;
            uint8_t u8;
            int16_t i16;
            uint16_t u16;
            int32_t i32;
            uint32_t u32;
            char * z;
            tag[0]=*cur++;
            tag[1]=*cur++;
            char val_type=*cur++;
            fprintf(stderr,"ttv: %c%c:%c\n", tag[0], tag[1], val_type);
            switch (val_type)
            {
              case 'A':
                  fprintf(stderr,"val='%c'\n",*cur++);
                  break;
              case 'c':
                  i8=*cur++;
                  break;
              case 'C':
                  u8=*cur++;
                  fprintf(stderr,"val=%d\n",u8);
                  break;
              case 's':
                  memmove(&i16,cur,2);
                  cur+=2;
                  break;
              case 'S':
                  memmove(&u16,cur,2);
                  cur+=2;
                  break;
              case 'i':
                  memmove(&i32,cur,4);
                  cur+=4;
                  break;
              case 'I':
                  memmove(&u32,cur,4);
                  cur+=4;
                  break;
              case 'f':
                  float f;
                  break;
              case 'Z':
                  z=cur;
                  while (isprint(*cur))
                      ++cur;
                  fprintf(stderr,"val='%s'\n",z);
                  ++cur;
                  break;
              case 'H':
                  z=cur;
                  while (isalnum(*cur))
                      ++cur;
                  fprintf(stderr,"val='%s'\n",z);
                  // TODO: Convert to ?
                  ++cur;
                  break;
              case 'B':
                  val_type=*cur++;
                  memmove(&u32,cur,4);
                  cur+=4;
                  cur+=u32*1; // TODO, based on size of val_type
                  break;
              default:
                  fprintf(stderr,"Bad val_type:%c\n", val_type);
                  return 1;
            }
        }
        fprintf(stderr,"no more ttvs\n");
        free(ttvs);
    }


    return 0;
}


void waitforthreads(Vector * threads)
{
    for (uint32_t i=0; i!=VectorLength(threads); ++i)
    {
        rc_t rc;
        fprintf(stderr,"waiting for thread %d to complete\n", i);
        KThreadWait((KThread*)VectorGet(threads,i),&rc);
    }
    fprintf(stderr,"all threads completed\n");
}

void releasethreads(Vector * threads)
{
    fprintf(stderr,"Release threads\n");
    for (uint32_t i=0; i!=VectorLength(threads); ++i)
    {
        KThreadRelease((KThread*)VectorGet(threads,i));
    }
}

bool threadinflate(void)
{
    rc_t rc;
    struct timeout_t tm;
    z_stream strm;

    const int numthreads=16;
    KQueue * inflatequeue;
    KQueue * blockqueue;
    Vector chunks;
    Vector threads;

    KQueueMake(&inflatequeue, numthreads);
    KQueueMake(&blockqueue, numthreads);
    VectorInit(&chunks,0,numthreads);
    VectorInit(&threads,0,numthreads);

    for (uint32_t i=0; i!=numthreads; ++i)
    {
        chunk * c=(chunk *)calloc(1,sizeof(chunk));
        VectorAppend(&chunks,NULL,c);

        KThread * thread=(KThread *)calloc(1,sizeof(KThread));
        KThreadMake(&thread, inflater, (void *)inflatequeue);
        VectorAppend(&threads, NULL, thread);
    }
    KThread * blockthread=(KThread *)calloc(1,sizeof(KThread));
    KThreadMake(&blockthread, blocker, (void *)blockqueue);
    VectorAppend(&threads, NULL, blockthread);

    mmapbuf_cur=mmapbuf;

    while(1)
    {
        memset(&strm,0,sizeof strm);
        if (!memcmp(mmapbuf_cur,"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00",28))
        {
            fprintf(stderr,"EOF marker found\n");
            break;
        }
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
        fprintf(stderr,"found gzip header\n");
        if (head.extra && head.extra_len==6 && head.extra[0]=='B' && head.extra[1]=='C' && head.extra[2]==2 && head.extra[3]==0)
        {
            uint16_t bsize=head.extra[4]+head.extra[5]*256;
            inflateEnd(&strm);
            fprintf(stderr,"extra:   %d\n",bsize);
            fprintf(stderr,"total_in:%d\n",strm.avail_in);
            fprintf(stderr,"buf   in:%p\n",mmapbuf_cur);
            fprintf(stderr,"next_in: %p\n",strm.next_in);
            fprintf(stderr,"offset:  %ld\n",mmapbuf_cur-mmapbuf);
            if (bsize<=28) 
            {
                fprintf(stderr,"small block found\n");
                break;
            }

            chunk * c;
            // Find unused chunk, busy polls
            uint32_t i=0;
            while (1)
            {
                c=(chunk *)VectorGet(&chunks,i);
                if (c->state==empty) break;
                ++i;
                if (i==VectorLength(&chunks))
                {
                    fprintf(stderr,"looping\n");
                    usleep(1);
                }
                i %= VectorLength(&chunks);
            }

            c->state=compressed;
            c->in=(Bytef*)mmapbuf_cur;
            c->insize=bsize;
            c->outsize=sizeof(c->out);

            while(1)
            {
                TimeoutInit(&tm, 1000); // 1 second
                rc=KQueuePush(inflatequeue, (void *)c, &tm);
                if ((int)GetRCObject(rc)== rcTimeout) 
                {
                    fprintf(stderr,"queue full\n");
//                    usleep(1);
                } else if (rc == 0)
                {
                    fprintf(stderr,"queued: %p %d:%d %d\n", c->in, c->insize, rc, rcTimeout);
                    break;
                } else
                {
                    fprintf(stderr,"something queue %d\n", rc);
                }
            }

            while(1)
            {
                TimeoutInit(&tm, 1000); // 1 second
                rc=KQueuePush(blockqueue, (void *)c, &tm);
                if ((int)GetRCObject(rc)== rcTimeout) 
                {
                    fprintf(stderr,"block queue full\n");
//                    usleep(1);
                } else if (rc == 0)
                {
                    fprintf(stderr,"block queued: %p %d:%d %d\n", c->in, c->insize, rc, rcTimeout);
                    break;
                } else
                {
                    fprintf(stderr,"block something queue %d\n", rc);
                }
            }

            mmapbuf_cur=mmapbuf_cur+bsize+1;
        } else
        {
            fprintf(stderr,"error: BAM required extra extension not found\n");
        }
    }

    KQueueSeal(inflatequeue);
    KQueueSeal(blockqueue);
    waitforthreads(&threads);
    KQueueRelease(inflatequeue);
    KQueueRelease(blockqueue);
    releasethreads(&threads);

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

