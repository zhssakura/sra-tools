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

#ifndef _h_samextract_lib_
#define _h_samextract_lib_
#include <klib/rc.h>
#include <klib/defs.h>
#include <klib/vector.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif
typedef struct Extractor
{
    char *mmapbuf;
    off_t mmapbuf_sz;
    char *mmapbuf_cur;

    Vector headers;
    Vector alignments;

    char * read;
    char * cigar;
    char * rname;
    uint32_t pos;

    char * tags; // Space delimited tags seen in current line
    char * seqnames;
    char * ids;
} Extractor;

typedef struct Header
{
    const char * headercode; // HD, SQ, RG, PG, CO
    const char * tag; // VN, SN, LN, ID, ...
    const char * value;
} Header;

typedef struct Alignment
{
    const char * read;
    const char * cigar;
    const char * rname;
    uint32_t pos;
} Alignment;

rc_t MakeExtractor(Extractor **state, const char * fname, uint32_t num_threads);
rc_t ReleaseExtractor(Extractor **state); // dtor

rc_t ExtractorGetHeaders(Extractor **state, Vector *headers);
rc_t ExtractorInvalidateHeaders(Extractor **state);

rc_t ExtractorGetAlignments(Extractor **state, Vector *alignments);
rc_t ExtractorInvalidateAlignments(Extractor **state);

#ifdef __cplusplus
}
#endif
#endif // __h_sam_extract_lib_

