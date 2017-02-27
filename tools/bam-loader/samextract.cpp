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
#include "samextract-lib.h"


rc_t CC UsageSummary(char const *name)
{
    fprintf(stderr, "Usage: %s file.{sb}am [file...]", name);
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
        Extractor * extractor;
        rc_t rc=MakeExtractor(&extractor, *(++argv), -1);
        printf("Made extractor\n");

        Vector headers;
        rc=ExtractorGetHeaders(&extractor, &headers);
        for (uint32_t i=0; i!=VectorLength(&headers); ++i)
        {
            Header * hdr=(Header *)VectorGet(&headers,i);
            printf("\tHeader%d: %s %s %s\n", i, hdr->headercode, hdr->tag, hdr->value);
        // Do stuff with headers
        }
        ExtractorInvalidateHeaders(&extractor);


        uint32_t vlen;
        do
        {
            Vector alignments;
            rc=ExtractorGetAlignments(&extractor, &alignments);
            vlen=VectorLength(&alignments);
            printf("Returned %d alignments\n",vlen);
            for (uint32_t i=0; i!=vlen; ++i)
            {
                Alignment * align=(Alignment *)VectorGet(&alignments,i);
                printf("\tAlignment%d: %s\n", i, align->read);
            // Do stuff with headers
            }
            printf("\n");
            ExtractorInvalidateAlignments(&extractor);
        } while (vlen);

        ReleaseExtractor(&extractor);
        printf("Done with file\n");
    }
    printf("KMain done\n");
    return 0;
}

