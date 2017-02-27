/* ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnologmsgy Information
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

#ifndef _h_samextract_
#define _h_samextract_
#include <klib/rc.h>
#include <klib/defs.h>
#include <klib/vector.h>
#include <zlib.h>
#include <sys/types.h>
#include "samextract-lib.h"

typedef enum chunk_state { empty, compressed, uncompressed } chunk_state;

typedef  int8_t i8;
typedef uint8_t u8;
typedef  int16_t i16;
typedef uint16_t u16;
typedef  int32_t i32;
typedef uint32_t u32;

typedef struct chunk_s
{
    Bytef * in;
    uInt insize;
    uInt outsize;
    chunk_state state;
    Bytef out[66000];
} chunk;

typedef struct bamalign
{
    i32 block_size;
    i32 refID;
    i32 pos;
    u32 bin_mq_nl;
    u32 flag_nc;
    i32 l_seq;
    i32 next_refID;
    i32 next_pos;
    i32 tlen;
} bamalign;

#ifdef __cplusplus
extern "C" {
#endif
    int moredata(char * buf, int * numbytes, size_t maxbytes);
    int SAMerror(const char *);
    void logmsg (const char * fname, int line, const char * func, const char * severity, const char * fmt, ...);
    void samload(char const path[]);
    extern bool SAM_parsebegin(Extractor * state);
    extern int SAM_parsestring(Extractor * state, char * str);
    extern int SAM_parsebuffer(Extractor * state, char * str, size_t size);
    extern void SAM_parseend(Extractor * state);

    extern Extractor * globstate;
#ifdef __cplusplus
}
#endif

//#ifndef DEBUG
#define DEBUG 0
//#endif

#define ERR(...) logmsg(__FILE__, __LINE__, __func__, "Error",  __VA_ARGS__)

#define WARN(...) logmsg(__FILE__, __LINE__, __func__, "Warning", __VA_ARGS__)
#define INFO(...) logmsg(__FILE__, __LINE__, __func__, "Info", __VA_ARGS__)
#define DBG(...) \
    do { if (DEBUG) logmsg(__FILE__, __LINE__, __func__, "Debug",  __VA_ARGS__); } while (0)

#define MIN(a,b)    (((a) < (b)) ? (a) : (b))

#endif
