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
#include <klib/rc.h>
#include <klib/defs.h>
#include <klib/vector.h>
#include <kproc/queue.h>
#include <kproc/thread.hpp>
#include <kproc/timeout.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <regex.h>
#include <stdint.h>
#include <unistd.h>
#include "samextract.h"
#include "samextract-lib.h"

extern "C" {

void logmsg (const char * fname, int line, const char * func, const char * severity, const char * fmt, ...)
{
    char * buf;
    size_t bufsize=0;
    FILE * buffd;

    const char * basename=strrchr(fname,'/');
    if (!basename) basename=strrchr(fname,'\\');
    if (basename) ++basename;
    if (!basename) basename=fname;
    va_list args;
    va_start(args, fmt);

    buffd=open_memstream(&buf,&bufsize); 
    if (buffd==NULL)
    {
        fprintf(stderr,"can't open memstream\n");
        return;
    }
    fprintf(buffd, "%s:", severity);
    vfprintf(buffd, fmt, args);
    va_end(args);
    fprintf(buffd, "\t[%s:%s():%d]\n", basename, func, line);
    fclose(buffd);
    size_t r=fwrite(buf, bufsize, 1, stderr);
    if (r!=1) fprintf(stderr,"previous %zd log message truncated\n", bufsize);
    free(buf);
}


rc_t inflater(const KThread * kt, void * in)
{
    KQueue *inflatequeue=(KQueue *)in;
    struct timeout_t tm;

    z_stream strm;
    pthread_t threadid=pthread_self();
    INFO("\tThread %lu started.",threadid);

    while (1)
    {
        void * where=NULL;
        DBG("\t\tthread %lu checking queue",threadid);
        TimeoutInit(&tm, 5000); // 5 seconds
        rc_t rc=KQueuePop(inflatequeue, &where, &tm);
        if (rc==0)
        {
            chunk * c=(chunk *)where;
            DBG("\t\tthread %lu chunk %p size %u", threadid, c->in, c->insize);
            if (c->state!=compressed)
            {
                ERR("Inflater bad state");
                return 1;
            }

            memset(&strm,0,sizeof strm);
            int zrc=inflateInit2(&strm, MAX_WBITS + 16); // Only gzip format
            switch (zrc)
            {
                case Z_OK:
                    break;
                case Z_MEM_ERROR:
                    ERR("Out of memory in zlib");
                    break;
                case Z_VERSION_ERROR:
                    ERR("zlib version is not compatible; need version %s but have %s", ZLIB_VERSION,zlibVersion());
                    break;
                case Z_STREAM_ERROR:
                    ERR("zlib stream error");
                    break;
                default:
                    ERR("zlib error %s",strm.msg);
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
                    DBG("\t\tthread %lu OK %d %d %lu", threadid, strm.avail_in, strm.avail_out,strm.total_out);
                    c->outsize=strm.total_out;
                    c->state=uncompressed;
                    break;
                case Z_MEM_ERROR:
                    ERR("error: Out of memory in zlib");
                    break;
                case Z_VERSION_ERROR:
                    ERR("zlib version is not compatible; need version %s but have %s", ZLIB_VERSION,zlibVersion());
                    break;
                case Z_STREAM_ERROR:
                    ERR("zlib stream error %s",strm.msg);
                    break;
                case Z_STREAM_END:
                    ERR("zlib stream end%s",strm.msg);
                    break;
                default:
                    ERR("inflate error %d %s",zrc, strm.msg);
                    break;
            }
            inflateEnd(&strm);
        } else if ((int)GetRCObject(rc)== rcTimeout) 
        {
            INFO("\t\tthread %lu queue empty",threadid);
            if (KQueueSealed(inflatequeue))
            {
                INFO("\t\tQueue sealed, thread %lu complete", threadid);
                return 0;
            }
//            usleep(1);
        } else if ((int)GetRCObject(rc)== rcData) 
        {
            INFO("\t\tthread %lu queue data (empty?)",threadid);
            break;
        } else
        {
            ERR("rc=%d",rc);
        }
    }

    INFO("\t\tthread %lu terminating.",threadid);
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
        DBG("\t\tBlocker thread checking blocker queue");
        if (c && c->state!=uncompressed)
        {
            ERR("\t\tblocker bad state");
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
                    ERR("\t\tBlocker bad state");
                    return false;
                }
                while (c->state!=uncompressed)
                {
                    DBG("\t\tBlocker busy");
                    usleep(1);
                }
                DBG("\t\tBlocker thread chunk %p size %u", c->in, c->insize);
                break;
            } else if ((int)GetRCObject(rc)== rcTimeout) 
            {
                INFO("\t\tBlocker thread queue empty");
                if (KQueueSealed(que))
                {
                    INFO("\t\tQueue sealed, Blocker thread complete");
                    return false;
                }
    //            usleep(1);
            } else if ((int)GetRCObject(rc)== rcData) 
            {
                INFO("\t\tBlocker thread queue data (empty?)");
                return false;
            } else
            {
                ERR("blocker rc=%d",rc);
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
    INFO("\tBlocker thread %lu started.",threadid);

    chunkview cv(blockqueue);

    char magic[4];
    if (!cv.getbytes(magic,4)) return 1;
    if (memcmp(magic,"BAM\x01",4))
    {
        ERR("BAM magic not found");
        return 2;
    }
    i32 l_text;
    if (!cv.getbytes(&l_text,4)) return 1;
    if (l_text<0)
    {
        ERR("error: invalid l_text");
        return 1;
    }

    char *text=(char *)calloc(1,l_text+1);
    if (!cv.getbytes(text,l_text)) return 1;
    text[l_text]='\0';

    DBG("SAM header:%s",text);
    // TODO: Parse SAM header, make sure no alignments

    i32 n_ref;
    if (!cv.getbytes(&n_ref,4)) return 1;
    if (n_ref<0)
    {
        ERR("error: invalid n_ref");
        return 1;
    }
    DBG("# references %d", n_ref);

    for (int i=0; i!=n_ref; ++i)
    {
        i32 l_name;
        if (!cv.getbytes(&l_name,4)) return 1;
        DBG("%d: l_name=%d",i, l_name);
        if (l_name < 0)
        {
            ERR("error: invalid reference name length");
            return 1;
        }
        if (l_name > 256)
        {
            WARN("warning: Long reference name");
            return 1;
        }
        char *name=(char *)calloc(1,l_name+1);
        if (!cv.getbytes(name,l_name)) return 1;
        DBG("%d: reference name %s",i, name);
        free(name); // TODO, persist?
        i32 l_ref;
        if (!cv.getbytes(&l_ref,4)) return 1;
        DBG("length of reference sequence %d=%d",i,l_ref);
    }

    while (1)
    {
        bamalign align;
        if (!cv.getbytes(&align,sizeof(align))) return 1;
        DBG("\n\n\nalignment block_size=%d refID=%d pos=%d", 
                align.block_size, 
                align.refID, 
                align.pos);

        if (align.block_size < 0)
        {
            ERR("error: invalid block_size");
            return 1;
        }
        if (align.pos < 0)
        {
            ERR("error: invalid pos");
            return 1;
        }
        if (align.refID < -1 || align.refID > n_ref)
        {
            ERR("error: bad refID");
            return 1;
        }
        if (align.next_refID < -1 || align.next_refID > n_ref)
        {
            ERR("error: bad next_refID");
            return 1;
        }
        DBG("align.bin_mq_nl=%d",align.bin_mq_nl);
        u16 bin=align.bin_mq_nl >> 16;
        uint8_t mapq=(align.bin_mq_nl >> 8) & 0xff;
        uint8_t l_read_name=align.bin_mq_nl & 0xff;
        DBG("bin=%d mapq=%d l_read_name=%d", bin, mapq, l_read_name);

        u16 flag=align.flag_nc >> 16;
        u16 n_cigar_op=align.flag_nc & 0xffff;
        DBG("flag=%x n_cigar_op=%d", flag, n_cigar_op);

        char * read_name=(char *)calloc(1,l_read_name);
        if (!cv.getbytes(read_name,l_read_name)) return 1;
        DBG("read_name=%s",read_name);

        static const char opmap[]="MIDNSHP=X???????";
        u32 * cigar=(u32 *)calloc(n_cigar_op,sizeof(u32));
        if (!cv.getbytes(cigar,n_cigar_op*4)) return 1;
        for (int i=0; i!=n_cigar_op; ++i)
        {
            i32 oplen=cigar[i] >> 4;
            i32 op=cigar[i] & 0xf;
            DBG("\tcigar %d=%x len=%d %d(%c)", i, cigar[i], oplen, op, opmap[op]);
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

        INFO("%d pairs in sequence",align.l_seq);
        for (int i=0; i!=align.l_seq; ++i)
            DBG(" seq#%d %c %.2x ",i, seq[i], qual[i]);
        DBG("\nseq=%s",seq);

        int remain=align.block_size-(sizeof(align)+l_read_name+n_cigar_op*4+bytesofseq+align.l_seq)+4; // TODO
        DBG("%d bytes remaining for ttvs",remain);
        char * ttvs=(char*)calloc(1,remain); // TODO: alloca? check <64K guarantee
        if (!cv.getbytes(ttvs,remain)) return 1;
        char * cur=ttvs;
        while (cur<ttvs+remain)
        {
            char tag[2];
            char c;
            i8 i8;
            u8 u8;
            i16 i16;
            u16 u16;
            i32 i32;
            u32 u32;
            char * z;
            tag[0]=*cur++;
            tag[1]=*cur++;
            char val_type=*cur++;
            DBG("ttv: %c%c:%c", tag[0], tag[1], val_type);
            switch (val_type)
            {
              case 'A':
                  c=*cur++;
                  DBG("val='%c'",c);
                  break;
              case 'c':
                  i8=*cur++;
                  DBG("val=%d",i8);
                  break;
              case 'C':
                  u8=*cur++;
                  DBG("val=%d",u8);
                  break;
              case 's':
                  memmove(&i16,cur,2);
                  DBG("val=%d",i16);
                  cur+=2;
                  break;
              case 'S':
                  memmove(&u16,cur,2);
                  DBG("val=%d",u16);
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
                  //float f;
                  break;
              case 'Z':
                  z=cur;
                  while (isprint(*cur))
                      ++cur;
                  DBG("val='%s'",z);
                  ++cur;
                  break;
              case 'H':
                  z=cur;
                  while (isalnum(*cur))
                      ++cur;
                  DBG("val='%s'",z);
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
                  ERR("Bad val_type:%c", val_type);
                  return 1;
            }
        }
        DBG("no more ttvs");
        free(ttvs);
    }


    return 0;
}


void waitforthreads(Vector * threads)
{
    for (u32 i=0; i!=VectorLength(threads); ++i)
    {
        rc_t rc;
        INFO("waiting for thread %d to complete", i);
        KThreadWait((KThread*)VectorGet(threads,i),&rc);
    }
    INFO("all threads completed");
}

void releasethreads(Vector * threads)
{
    INFO("Release threads");
    for (u32 i=0; i!=VectorLength(threads); ++i)
    {
        KThreadRelease((KThread*)VectorGet(threads,i));
    }
}

bool threadinflate(Extractor * state, int numthreads)
{
    rc_t rc;
    struct timeout_t tm;
    z_stream strm;

    if (numthreads==-1)
        numthreads=sysconf(_SC_NPROCESSORS_ONLN) - 1;
    KQueue * inflatequeue;
    KQueue * blockqueue;
    Vector chunks;
    Vector threads;

    KQueueMake(&inflatequeue, numthreads);
    KQueueMake(&blockqueue, numthreads);
    VectorInit(&chunks,0,numthreads);
    VectorInit(&threads,0,numthreads);

    for (int i=0; i!=numthreads; ++i)
    {
        chunk * c=(chunk *)calloc(1,sizeof(chunk));
        c->state=empty;
        VectorAppend(&chunks,NULL,c);

        KThread * thread=(KThread *)calloc(1,sizeof(KThread));
        KThreadMake(&thread, inflater, (void *)inflatequeue);
        VectorAppend(&threads, NULL, thread);
    }
    KThread * blockthread=(KThread *)calloc(1,sizeof(KThread));
    KThreadMake(&blockthread, blocker, (void *)blockqueue);
    VectorAppend(&threads, NULL, blockthread);

    while(1)
    {
        memset(&strm,0,sizeof strm);
        if (!memcmp(state->mmapbuf_cur,"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00",28))
        {
            INFO("EOF marker found");
            break;
        }
        strm.next_in=(Bytef*)state->mmapbuf_cur;
        if (state->mmapbuf_sz > INT32_MAX)
            strm.avail_in=INT32_MAX;
        else
            strm.avail_in=(uInt)state->mmapbuf_sz;
        int zrc=inflateInit2(&strm, MAX_WBITS + 16); // Only gzip format
        switch (zrc)
        {
            case Z_OK:
                break;
            case Z_MEM_ERROR:
                ERR("error: Out of memory in zlib");
                break;
            case Z_VERSION_ERROR:
                ERR("zlib version is not compatible; need version %s but have %s", ZLIB_VERSION,zlibVersion());
                break;
            case Z_STREAM_ERROR:
                ERR("zlib stream error");
                break;
            default:
                ERR("zlib error");
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
                ERR("inflate error");
            }
        }
        INFO("found gzip header");
        if (head.extra && head.extra_len==6 && head.extra[0]=='B' && head.extra[1]=='C' && head.extra[2]==2 && head.extra[3]==0)
        {
            u16 bsize=head.extra[4]+head.extra[5]*256;
            inflateEnd(&strm);
            DBG("extra:   %d",bsize);
            DBG("total_in:%d",strm.avail_in);
            DBG("buf   in:%p",state->mmapbuf_cur);
            DBG("next_in: %p",strm.next_in);
            DBG("offset:  %ld",state->mmapbuf_cur-state->mmapbuf);
            if (bsize<=28) 
            {
                INFO("small block found");
                break;
            }

            chunk * c;
            // Find unused chunk, busy polls
            u32 i=0;
            while (1)
            {
                c=(chunk *)VectorGet(&chunks,i);
                if (c->state==empty) break;
                ++i;
                if (i==VectorLength(&chunks))
                {
                    DBG("looping");
                    usleep(1);
                }
                i %= VectorLength(&chunks);
            }

            c->state=compressed;
            c->in=(Bytef*)state->mmapbuf_cur;
            c->insize=bsize;
            c->outsize=sizeof(c->out);

            while(1)
            {
                TimeoutInit(&tm, 1000); // 1 second
                rc=KQueuePush(inflatequeue, (void *)c, &tm);
                if ((int)GetRCObject(rc)== rcTimeout) 
                {
                    DBG("queue full");
//                    usleep(1);
                } else if (rc == 0)
                {
                    DBG("queued: %p %d:%d %d", c->in, c->insize, rc, rcTimeout);
                    break;
                } else
                {
                    DBG("something queue %d", rc);
                }
            }

            while(1)
            {
                TimeoutInit(&tm, 1000); // 1 second
                rc=KQueuePush(blockqueue, (void *)c, &tm);
                if ((int)GetRCObject(rc)== rcTimeout) 
                {
                    DBG("block queue full");
//                    usleep(1);
                } else if (rc == 0)
                {
                    DBG("block queued: %p %d:%d %d", c->in, c->insize, rc, rcTimeout);
                    break;
                } else
                {
                    DBG("block something queue %d", rc);
                }
            }

            state->mmapbuf_cur=state->mmapbuf_cur+bsize+1;
        } else
        {
            ERR("error: BAM required extra extension not found");
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


rc_t MakeExtractor(Extractor **state, const char * fname, uint32_t num_threads=-1)
{
    Extractor * s=(Extractor *)calloc(1,sizeof(Extractor));
    *state=s;
    s->mmapbuf=NULL;
    s->mmapbuf_sz=0;
    s->mmapbuf_cur=NULL;

    VectorInit(&s->headers,0,0);
    VectorInit(&s->alignments,0,0);

    s->tags=strdup("");
    s->seqnames=strdup("");
    s->ids=strdup("");

    int fd=open(fname, O_RDONLY);
    if (fd==-1)
    {
        ERR("Cannot open %s:%s", fname, strerror(errno));
    }

    struct stat st;
    if (fstat(fd,&st))
    {
        ERR("Cannot stat %s:%s", fname, strerror(errno));
    }
    s->mmapbuf_sz=st.st_size;
    DBG("stat_sz = %ld", s->mmapbuf_sz);

    s->mmapbuf=(char *)mmap(NULL, (size_t)s->mmapbuf_sz, PROT_READ, MAP_PRIVATE, fd, 0);
    s->mmapbuf_cur=s->mmapbuf;
    close(fd); // Still need to munmap

//    madvise(mmapbuf, mmapbuf_sz, MADV_SEQUENTIAL | MADV_WILLNEED);
    madvise(s->mmapbuf, s->mmapbuf_sz, MADV_SEQUENTIAL);

    if (s->mmapbuf==MAP_FAILED)
    {
        ERR("Cannot mmap %s:%s", fname, strerror(errno));
        return 1;
    }

    if (!memcmp(s->mmapbuf,"\x1f\x8b\x08",3))
    {
        INFO("gzip file");
        threadinflate(s,num_threads);
    } else if (s->mmapbuf[0]=='@')
    {
        INFO("SAM file");

        SAM_parsebegin(s);
        while (1)
        {
            ssize_t remain=s->mmapbuf_sz-(s->mmapbuf_cur-s->mmapbuf);
            if (remain<=0) 
            {
                INFO("Buffer complete");
                break;
            }
            if (s->mmapbuf_cur[0]!='@') 
            {
                INFO("out of headers");
                break;
            }
            char * nl=(char *)memchr(s->mmapbuf_cur, '\x0a', remain);
            if (nl)
            {
                ++nl;
                size_t linesize=nl-s->mmapbuf_cur;
                DBG("linesize=%d",linesize);
                SAM_parsebuffer(s,s->mmapbuf_cur,linesize);
                s->mmapbuf_cur+=linesize;
            } else
            {
                INFO("No newline");
                // No newline, TODO: final line terminated?
                break;
            }
        }
        DBG("Done parsing headers");
    }

    return 0;
}

rc_t ReleaseExtractor(Extractor **state)
{
    DBG("release_Extractor");
    // TODO: invalidate
    Extractor *s=*state;
    SAM_parseend(s);

    munmap(s->mmapbuf,s->mmapbuf_sz);
    s->mmapbuf=NULL;
    s->mmapbuf_sz=0;
    s->mmapbuf_cur=NULL;

    VectorWhack(&s->headers,NULL,NULL);
    VectorWhack(&s->alignments,NULL,NULL);

    free(s->tags);
    s->tags=NULL;
    free(s->seqnames);
    s->seqnames=NULL;
    free(s->ids);
    s->ids=NULL;

    free(s);

    return 0;
}

rc_t ExtractorGetHeaders(Extractor **state, Vector *headers)
{
    Extractor * s=*state;
    DBG("get_headers");
    VectorInit(headers,0,0);
    for (u32 i=0; i!=VectorLength(&s->headers); ++i)
    {
        Header * hdr=(Header *)VectorGet(&s->headers,i);
        VectorAppend(headers,NULL,hdr);
    }
    return 0;
}

rc_t ExtractorInvalidateHeaders(Extractor **state)
{
    Extractor * s=*state;
    DBG("invalidate_headers");
    for (u32 i=0; i!=VectorLength(&s->headers); ++i)
    {
        Header * hdr=(Header *)VectorGet(&s->headers,i);
        free((void*)hdr->headercode);
        hdr->headercode=NULL;
        free((void*)hdr->tag);
        hdr->tag=NULL;
        free((void*)hdr->value);
        hdr->value=NULL;
        free(hdr);
    }
    VectorInit(&s->headers,0,0);
    return 0;
}

rc_t ExtractorGetAlignments(Extractor **state, Vector *alignments)
{
    Extractor * s=*state;
    DBG("get_alignments");
    VectorInit(alignments,0,0);
    VectorInit(&s->alignments,0,0);
    int numaligns=20;
    while (numaligns--)
    {
        ssize_t remain=s->mmapbuf_sz-(s->mmapbuf_cur-s->mmapbuf);
        if (remain<=0) 
        {
            INFO("Buffer complete");
            break;
        }
        if (s->mmapbuf_cur[0]=='@') 
        {
            INFO("headers restarted");
            break;
        }
        char * nl=(char *)memchr(s->mmapbuf_cur, '\x0a', remain);
        if (nl)
        {
            ++nl;
            size_t linesize=nl-s->mmapbuf_cur;
            DBG("alignment #%d linesize=%d",numaligns,linesize);
            SAM_parsebuffer(s,s->mmapbuf_cur,linesize);
            s->mmapbuf_cur+=linesize;
        } else
        {
            INFO("No newline");
            // No newline, TODO: final line terminated?
            break;
        }
    }

    VectorWhack(alignments,NULL,NULL);
    VectorCopy(&s->alignments,alignments);
    VectorWhack(&s->alignments,NULL,NULL);
    DBG("got_alignments");

    return 0;
}

rc_t ExtractorInvalidateAlignments(Extractor **state)
{
    DBG("invalidate_alignments");
    Extractor * s=*state;
    for (uint32_t i=0; i!=VectorLength(&s->alignments); ++i)
    {
        Alignment * align=(Alignment *)VectorGet(&s->alignments,i);
        free((void*)align->read);
        align->read=NULL;
        free((void*)align->cigar);
        align->cigar=NULL;
        free((void*)align->rname);
        align->rname=NULL;
        align->pos=0;
        free(align);
        align=NULL;
        VectorSet(&s->alignments,i,align);
    }
    VectorWhack(&s->alignments,NULL,NULL);
    VectorInit(&s->alignments,0,0);

    return 0;
}

} // extern C

