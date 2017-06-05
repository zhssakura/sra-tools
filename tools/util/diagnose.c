/*==============================================================================
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
* ==============================================================================
*
*/
#include "test-sra-priv.h" /* endpoint_to_string */
#include <kfs/file.h> /* KFile */
#include <klib/out.h> /* KOutMsg */
#include <klib/printf.h> /* string_vprintf */
#include <klib/rc.h>
#include <klib/text.h> /* String */
#include <kns/endpoint.h> /* KNSManagerInitDNSEndpoint */
#include <kns/http.h> /* KHttpRequest */
#include <kns/manager.h> /* KNSManager */
#include <kns/kns-mgr-priv.h> /* KNSManagerMakeReliableHttpFile */
#include <kns/stream.h> /* KStream */
#include <ctype.h> /* isprint */
#include <limits.h> /* PATH_MAX */
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

typedef struct {
    int n [ 5 ];
    int ended;
} STest;
static void STestInit ( STest * self ) {
    assert ( self );
    memset ( self, 0, sizeof * self );
    self -> ended = -1;
}
static rc_t STestVStart ( STest * self, bool checking,
                          const char * fmt, va_list args  )
{
    char b [ 512 ] = "";
    rc_t rc = string_vprintf ( b, sizeof b, NULL, fmt, args );
    if ( rc != 0 )
        OUTMSG ( ( "CANNOT PRINT: %R", rc ) );
    else {
        assert ( self );
        int i = 0;
        if ( self -> ended == -1 ) {
            self -> ended = sizeof self -> n / sizeof self -> n [ 0 ] - 1;
            for ( i = 0; i < sizeof self -> n / sizeof self -> n [ 0 ]; ++ i )
                if ( self -> n [ i ] == 0 ) {
                    self -> ended = i;
                    break;
                }
        }
        assert ( self -> ended >= 0 );
        ++ self -> n [ self -> ended ];
        self -> ended = -1;
        OUTMSG ( ( "> %d", self -> n [ 0 ] ) );
        for ( i = 1; i < sizeof self -> n / sizeof self -> n [ 0 ]; ++ i )
            if ( self -> n [ i ] == 0 )
                break;
            else
                OUTMSG ( ( ".%d", self -> n [ i ] ) );
    }
    if ( rc == 0 )
        rc = OUTMSG ( ( " %s%s%s", checking ? "Checking " : "", b,
                        checking ? "...\n" : " " ) );
    return rc;
}
typedef enum {
    eFAIL,
    eOK,
    eMGS,
    eEND,
    eDONE,
} EOK;
static rc_t STestVEnd ( STest * self, EOK ok,
                        const char * fmt, va_list args )
{
    int i = 0;
    assert ( self );
    if ( ok != eMGS ) {
        if ( self -> ended != -1 )
            -- self -> ended;
        else {
            self -> ended = sizeof self -> n / sizeof self -> n [ 0 ] - 1;
            for ( i = 0; i < sizeof self -> n / sizeof self -> n [ 0 ]; ++ i )
                if ( self -> n [ i ] == 0 ) {
                    self -> ended = i - 1;
                    break;
                }
        }
    }
    char b [ 512 ] = "";
    rc_t rc = string_vprintf ( b, sizeof b, NULL, fmt, args );
    if ( rc != 0 )
        OUTMSG ( ( "CANNOT PRINT: %R", rc ) );
    else {
        if ( ok == eFAIL || ok == eOK || ok == eDONE ) {
            rc = OUTMSG ( ( "< %d", self -> n [ 0 ] ) );
            for ( i = 1; i <= self -> ended; ++ i )
                OUTMSG ( ( ".%d", self -> n [ i ] ) );
            for ( ; i < sizeof self -> n / sizeof self -> n [ 0 ]; ++ i )
                if ( self -> n [ i ] != 0 )
                    self -> n [ i ] = 0;
                else
                    break;
            OUTMSG ( ( " " ) );
        }
        OUTMSG ( ( b ) );
        switch ( ok ) {
            case eFAIL: OUTMSG ( ( ": FAILURE\n" ) ); break;
            case eOK  : OUTMSG ( ( ": OK\n"      ) ); break;
            case eEND :
            case eDONE: OUTMSG ( (     "\n"      ) ); break;
            default   :                               break;
        }
    }
    return rc;
}
static rc_t STestEnd ( STest * self, EOK ok, const char * fmt, ...  )  {
    va_list args;
    va_start ( args, fmt );
    rc_t rc = STestVEnd ( self, ok, fmt, args );
    va_end ( args );
    return rc;
}
static rc_t STestStart ( STest * self, bool checking,
                         const char * fmt, ...  )
{
    va_list args;
    va_start ( args, fmt );
    rc_t rc = STestVStart ( self, checking, fmt, args );
    va_end ( args );
    return rc;
}
static rc_t STestCheckUrl ( STest * self, const KNSManager * mgr,
    const String * host, const char * path, bool print,
    const char * exp, size_t esz )
{
    char buffer [ 1024 ] = "";
    char full [ PATH_MAX ] = "";
    ver_t http_vers = 0x01010000;
    KHttpRequest * req = NULL;
    KHttpResult * rslt = NULL;
    size_t num_read = 0;
    assert ( path );
    bool have_size = path [ 0 ] != '\0';
    bool not_exist = false;
    String gap;
    CONST_STRING ( & gap, "gap-download.ncbi.nlm.nih.gov" );
    if ( StringEqual ( host, & gap ) )
        not_exist = true;
    if ( have_size ) {
        const char cgi [] = "names.cgi";
        size_t csz = sizeof cgi - 1;
        uint32_t m =  string_measure ( path, NULL );
        if ( m > csz && string_cmp ( path + m - csz, csz, cgi, csz, csz ) == 0 )
        {
            have_size = false;
        }
    }
    rc_t rc = string_printf ( full, sizeof full, NULL,
                              "https://%S%s", host, path );
    if ( rc != 0 )
        OUTMSG ( ( "CANNOT PRINT PATH: %R\n", rc ) );
    if ( rc == 0 )
        rc = STestStart ( self, true, "Access to %s", full );
    uint64_t sz = 0;
    if ( rc == 0 ) {
        const KFile * file = NULL;
        rc = STestStart ( self, false,
                          "KFile = KNSManagerMakeReliableHttpFile(%s):", full );
        if ( rc == 0 ) {
            rc = KNSManagerMakeReliableHttpFile ( mgr, & file, NULL,
                                                  http_vers, full);
            if ( rc != 0 ) {
                if ( ( not_exist || ! have_size ) &&
                 (rc == SILENT_RC ( rcNS, rcFile, rcOpening, rcFile, rcNotFound)
                  ||
                  rc == SILENT_RC ( rcNS, rcFile, rcOpening, rcSize, rcUnknown))
                )
                {
                    rc = 0;
                    STestEnd ( self, eEND, "Skipped (does not exist)" );
                }
                else
                    STestEnd ( self, eEND, "FAILURE: %R", rc );
            }
            else {
                if ( rc == 0 ) {
                    STestEnd ( self, eEND, "OK" );
                    rc = STestStart ( self, false,
                          "KFileSize(KFile(%s)) =", full );
                    rc = KFileSize ( file, & sz );
                    if ( rc == 0 )
                        STestEnd ( self, eEND, "%lu: OK", sz );
                    else
                        STestEnd ( self, eEND, "FAILURE: %R", rc );
                }
                else if ( have_size )
                    STestEnd ( self, eEND, "FAILURE: %R", rc );
                else if ( rc == SILENT_RC ( rcNS,
                                    rcFile, rcOpening, rcSize, rcUnknown ) )
                {
                    rc = 0;
                    STestEnd ( self, eEND, "Skipped (Size Unknown)" );
                }
                else
                    STestEnd ( self, eEND, "FAILURE: Size Unknown but %R", rc );
            }
        }
        KFileRelease ( file );
        file = NULL;
    }
    if ( rc == 0 ) {
        bool skipped = false;
        uint64_t pos = 0;
        size_t bytes = 4096;
        size_t ebytes = bytes;
        if ( sz < ebytes )
            ebytes = sz;
        if ( sz > bytes * 2 )
            pos = sz / 2;
        KClientHttp * http = NULL;
        rc = STestStart ( self, true, "Support of Range requests" );
        if ( rc == 0 ) {
            rc = STestStart ( self, false,
                "KClientHttp = KNSManagerMakeClientHttps(%S):", host );
            if ( rc == 0 ) {
                rc = KNSManagerMakeClientHttps ( mgr, & http, NULL,
                                                 http_vers, host, 0 );
                if ( rc == 0 )
                    STestEnd ( self, eEND, "OK" );
                else
                    STestEnd ( self, eEND, "FAILURE: %R", rc );
            }
        }
        if ( rc == 0 ) {
            if ( path [ 0 ] == '\0' )
                path = "/";
            rc = KHttpMakeRequest ( http, & req, path );
        }
        if ( rc == 0 ) {
            rc = STestStart ( self, false, "KHttpResult = "
                "KHttpRequestHEAD(KHttpMakeRequest(KClientHttp)):" );
            if ( rc == 0 ) {
                rc = KHttpRequestHEAD ( req, & rslt );
                if ( rc == 0 )
                    STestEnd ( self, eEND, "OK" );
                else
                    STestEnd ( self, eEND, "FAILURE: %R", rc );
            }
        }
        if ( rc == 0 ) {
            rc = STestStart ( self, false,
                "KHttpResultGetHeader(KHttpResult, Accept-Ranges) =" );
            if ( rc == 0 ) {
                rc = KHttpResultGetHeader ( rslt, "Accept-Ranges",
                                            buffer, sizeof buffer, & num_read );
                if ( rc == 0 ) {
                    const char bytes [] = "bytes";
                    if ( string_cmp ( buffer, num_read, bytes, sizeof bytes - 1,
                            sizeof bytes - 1 ) == 0 )
                    {
                        STestEnd ( self, eEND, "'%.*s': OK",
                                   ( int ) num_read, buffer );
                    }
                    else {
                        STestEnd ( self, eEND, "'%.*s': FAILURE",
                                   ( int ) num_read, buffer );
                        rc = RC ( rcExe, rcFile,
                                  rcOpening, rcFunction, rcUnsupported );
                    }
                }
                else if ( have_size )
                    STestEnd ( self, eEND, "FAILURE: %R", rc );
                else if ( sz == 0 && rc == SILENT_RC
                             ( rcNS, rcTree, rcSearching, rcName, rcNotFound ) )
                {
                    rc = 0;
                    skipped = true;
                    STestEnd ( self, eEND, "Skipped (Not Required for %S)",
                                           host );
                }
                else
                    STestEnd ( self, eEND, "FAILURE: Not Required but %R", rc );
            }
        }
        KHttpResultRelease ( rslt );
        rslt = NULL;
        if ( rc == 0 ) {
            rc = STestStart ( self, false, "KHttpResult = KHttpRequestByteRange"
                "(KHttpMakeRequest, %lu, %zu):", pos, bytes );
            if ( rc == 0 ) {
                rc = KHttpRequestByteRange ( req, pos, bytes );
                if ( rc == 0 )
                    STestEnd ( self, eEND, "OK" );
                else
                    STestEnd ( self, eEND, "FAILURE: %R", rc );
            }
        }
        if ( rc == 0 ) {
            rc = STestStart ( self, false, "KHttpResult = "
                "KHttpRequestGET(KHttpMakeRequest(KClientHttp)):" );
            if ( rc == 0 ) {
                rc = KHttpRequestGET ( req, & rslt );
                if ( rc == 0 )
                    STestEnd ( self, eEND, "OK" );
                else
                    STestEnd ( self, eEND, "FAILURE: %R", rc );
            }
        }
        if ( rc == 0 ) {
            uint64_t po = 0;
            size_t byte = 0;
            rc = KClientHttpResultRange ( rslt, & po, & byte );
            if ( rc == 0 ) {
                if ( po != pos || ( ebytes > 0 && byte != ebytes ) ) {
                    STestStart ( self, false,
                                 "KClientHttpResultRange(KHttpResult,&p,&b):" );
                    STestEnd ( self, eEND, "FAILURE: expected:{%lu,%zu}, "
                        "got:{%lu,%zu}", pos, ebytes, po, byte );
                    rc = RC ( rcExe, rcFile, rcReading, rcRange, rcOutofrange );
                }
            }
            else if ( ! have_size && sz == 0 && rc == SILENT_RC
                    ( rcNS, rcNoTarg, rcValidating, rcError, rcUnsupported ) )
            {
                rc = 0;
                skipped = true;
            }
            else {
                STestStart ( self, false,
                             "KClientHttpResultRange(KHttpResult):" );
                STestEnd ( self, eEND, "FAILURE: %R", rc );
            }
        }
        if ( rc == 0 ) {
            rc = STestStart ( self, false,
                "KHttpResultGetHeader(KHttpResult, Content-Range) =" );
            if ( rc == 0 ) {
                rc = KHttpResultGetHeader ( rslt, "Content-Range",
                                            buffer, sizeof buffer, & num_read );
                if ( rc == 0 )
                    STestEnd ( self, eEND, "'%.*s': OK",
                               ( int ) num_read, buffer );
                else if ( have_size )
                    STestEnd ( self, eEND, "FAILURE: %R", rc );
                else if ( skipped && sz == 0 && rc == SILENT_RC
                             ( rcNS, rcTree, rcSearching, rcName, rcNotFound ) )
                {
                    rc = 0;
                    STestEnd ( self, eEND, "Skipped (Not Required for %S)",
                                           host );
                }
                else
                    STestEnd ( self, eEND, "FAILURE: Not Required but %R", rc );
            }
        }
        KHttpResultRelease ( rslt );
        rslt = NULL;
        KHttpRequestRelease ( req );
        req = NULL;
        KHttpRelease ( http );
        http = NULL;
        if ( rc == 0 && skipped )
            STestEnd ( self, eDONE,
                       "Support of Range requests: Skipped (Not Required)" );
        else
            STestEnd ( self, rc == 0 ? eOK : eFAIL,
                       "Support of Range requests" );
    }
    if ( rc == 0 ) {
        rc = STestStart ( self, false,
                          "KHttpRequest = KNSManagerMakeRequest(%s):", full );
        if ( rc == 0 ) {
            rc = KNSManagerMakeRequest ( mgr, & req, http_vers, NULL, full );
            if ( rc == 0 )
                STestEnd ( self, eEND, "OK"  );
            else
                STestEnd ( self, eEND, "FAILURE: %R", rc );
        }
    }
    if ( rc == 0 ) {
        rc = STestStart ( self, false,
            "KHttpResult = KHttpRequestGET(KHttpRequest):" );
        if ( rc == 0 ) {
            rc = KHttpRequestGET ( req, & rslt );
            if ( rc == 0 )
                STestEnd ( self, eEND, "OK" );
            else
                STestEnd ( self, eEND, "FAILURE: %R", rc );
        }
    }
    if ( rc == 0 ) {
        uint32_t code = 0;
        rc = STestStart ( self, false, "KHttpResultStatus(KHttpResult) =" );
        if ( rc == 0 ) {
            rc = KHttpResultStatus ( rslt, & code, NULL, 0, NULL );
            if ( rc != 0 )
                STestEnd ( self, eEND, "FAILURE: %R", rc );
            else {
                OUTMSG ( ( "%u: ", code ) );
                if ( code == 200 )
                    STestEnd ( self, eEND, "OK" );
                else if ( code == 404 && not_exist )
                    STestEnd ( self, eEND, "Skipped (does not exist)" );
                else {
                    STestEnd ( self, eEND, "FAILURE" );
                    rc = RC ( rcExe, rcFile, rcReading, rcFile, rcInvalid );
                }
            }
        }
    }
    if ( rc == 0 ) {
        KStream * response = NULL;
        rc = KHttpResultGetInputStream ( rslt, & response );
        if ( rc != 0 )
            OUTMSG ( (
                "KHttpResultGetInputStream(KHttpResult) = %R\n", rc ) );
        else {
            size_t total = 0;
            rc = STestStart ( self, false, "KStreamRead(KHttpResult):" );
            while ( rc == 0 ) {
                rc = KStreamRead ( response, buffer, sizeof buffer,
                                   & num_read );
                if ( rc != 0 )
                    if ( not_exist && rc ==
                        SILENT_RC ( rcNS, rcFile, rcReading, rcSelf, rcNull ) )
                    {
                        STestEnd ( self, eEND, "Skipped (does not exist)" );
                        rc = 0;
                        break;
                    }
                    else
                        STestEnd ( self, eEND, "FAILURE: %R", rc );
                else if ( num_read != 0 ) {
                    if ( total == 0 && esz > 0 ) {
                        int i = 0;
                        int s = esz;
                        if ( num_read < esz )
                            s = num_read;
                        STestEnd ( self, eMGS, "'" );
                        for ( i = 0; i < s; ++ i ) {
                            if ( isprint ( buffer [ i ] ) )
                                STestEnd ( self, eMGS, "%c", buffer [ i ] );
                            else if ( buffer [ i ] == 0 )
                                STestEnd ( self, eMGS, "\\0" );
                            else
                                STestEnd ( self, eMGS, "\\%03o",
                                           ( unsigned char ) buffer [ i ] );
                        }
                        STestEnd ( self, eMGS, "': " );
                        if ( string_cmp ( buffer, num_read, exp, esz, esz )
                                != 0 )
                        {
                            STestEnd ( self, eEND, " FAILURE: bad content" );
                            rc = RC ( rcExe,
                                      rcFile, rcReading, rcString, rcUnequal );
                        }
                    }
                    total += num_read;
                }
                else {
                    if ( ! have_size && sz == 0 )
                        sz = total;
                    if ( total == sz ) {
                        if ( print ) {
                            if ( total >= sizeof buffer )
                                buffer [ sizeof buffer - 1 ] = '\0';
                            else {
                                buffer [ total ] = '\0';
                                while ( total > 0 ) {
                                    -- total;
                                    if ( buffer [ total ] == '\n' )
                                        buffer [ total ] = '\0';
                                    else
                                        break;
                                }
                            }
                            STestEnd ( self, eMGS, "%s: ", buffer );
                        }
                        STestEnd ( self, eEND, "OK" );
                    }
                    else
                        STestEnd ( self, eEND,
                                   "%s: SIZE DO NOT MATCH (%zu)\n", total );
                    break;
                }
            }
        }
        KStreamRelease ( response );
        response = NULL;
    }
    STestEnd ( self, rc == 0 ? eOK : eFAIL, "Access to %s", full );
    return rc;
}
/******************************************************************************/
static rc_t STestCheckNetwork ( STest * self, const KNSManager * mgr,
    const String * domain, const char * path, const char * exp, size_t esz,
    const char * path2, const char * fmt, ... )
{
    uint16_t port = 443;
    KEndPoint ep;
    va_list args;
    va_start ( args, fmt );
    char b [ 512 ] = "";
    rc_t rc = string_vprintf ( b, sizeof b, NULL, fmt, args );
    va_end ( args );
    rc = STestStart ( self, true, b );
    if ( rc == 0 )
        rc = STestStart ( self, false, "KNSManagerInitDNSEndpoint(%S:%hu) =",
                          domain, port );
    if ( rc == 0 ) {
        rc = KNSManagerInitDNSEndpoint ( mgr, & ep, domain, port );
        if ( rc != 0 )
            STestEnd ( self, eEND, "FAILURE: %R", rc );
        else {
            char endpoint [ 1024 ] = "";
            rc_t rx = endpoint_to_string ( endpoint, sizeof endpoint, & ep );
            if ( rx != 0 )
                STestEnd ( self, eEND, "CANNOT CONVERT TO STRING" );
            else
                STestEnd ( self, eEND, "'%s': OK", endpoint );
        }
    }
    if ( rc == 0 ) {
        assert ( path );
        rc = STestCheckUrl ( self, mgr, domain, path, false, exp, esz );
        if ( path2 != NULL ) {
            rc_t r2 = STestCheckUrl ( self, mgr, domain, path2, true, 0, 0 );
            if ( rc == 0 )
                rc = r2;
        }
    }
    STestEnd ( self, rc == 0 ? eOK : eFAIL, b );
    return rc;
}
rc_t MainQuickCheck ( const KNSManager * mgr ) {
    const char exp [] = "NCBI.sra\210\031\003\005\001\0\0\0";
    rc_t rc = 0;
    rc_t r1 = 0;
    STest t;
    STestInit ( & t );
    rc = STestStart ( & t, true, "Network" );
    {
        String d;
        CONST_STRING ( & d, "gap-download.ncbi.nlm.nih.gov" );
        rc_t r2 = STestCheckNetwork ( & t, mgr, & d, "", NULL, 0, 
                                      NULL, "Access to %S", & d );
        if ( r1 == 0 )
            r1 = r2;
    }
    {
        String d;
        CONST_STRING ( & d, "ftp-trace.ncbi.nlm.nih.gov" );
        rc_t r2 = STestCheckNetwork ( & t, mgr, & d,
            "/sra/refseq/KC702174.1", exp, sizeof exp - 1,
            "/sra/sdk/current/sratoolkit.current.version",
            "Access to %S", & d );
        if ( r1 == 0 )
            r1 = r2;
    }
    {
        String d;
        CONST_STRING ( & d, "sra-download.ncbi.nlm.nih.gov" );
        rc_t r2 = STestCheckNetwork ( & t, mgr, & d, "/srapub/SRR042846",
            exp, sizeof exp - 1, NULL, "Access to %S", & d );
        if ( r1 == 0 )
            r1 = r2;
    }
    {
        String d;
        CONST_STRING ( & d, "www.ncbi.nlm.nih.gov" );
        rc_t r2 = STestCheckNetwork ( & t, mgr, & d, "", 0, 0,
            "/Traces/names/names.cgi", "Access to %S", & d );
        if ( r1 == 0 )
            r1 = r2;
    }
    STestEnd ( & t, r1 == 0 ? eOK : eFAIL, "Network" );
  if(0){
    OUTMSG ( ( "> 1.1 Checking Access to ftp-trace...\n" ) );
    KHttpRequest * req = NULL;
    OUTMSG ( ( "KHttpRequest = KNSManagerMakeRequest("
"https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current.version"
        ")... " ) );
    rc_t r2 = KNSManagerMakeRequest ( mgr, & req, 0x01010000, NULL,
"https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current.version"
        );
    if ( r2 == 0 )
        OUTMSG ( ( "OK\n" ) );
    else
        OUTMSG ( ( "FAILURE: %R\n", r2 ) );
    KHttpResult * rslt = NULL;
    if ( r2 == 0 ) {
        OUTMSG ( ( "KHttpResult = KHttpRequestGET(KHttpRequest)... " ) );
        r2 = KHttpRequestGET ( req, & rslt );
        if ( r2 == 0 )
            OUTMSG ( ( "OK\n" ) );
        else
            OUTMSG ( ( "FAILURE: %R\n", r2 ) );
    }
    if ( r2 == 0 ) {
        uint32_t code = 0;
        OUTMSG ( ( "KHttpResultStatus(KHttpResult) = " ) );
        r2 = KHttpResultStatus ( rslt, & code, NULL, 0, NULL );
        if ( r2 != 0 )
            OUTMSG ( ( "FAILURE: %R\n", r2 ) );
        else {
            OUTMSG ( ( "%u: ", code ) );
            if ( code == 200 )
                OUTMSG ( ( "OK\n" ) );
            else {
                OUTMSG ( ( "FAILURE\n" ) );
                r2 = RC ( rcExe, rcFile, rcReading, rcFile, rcInvalid );
            }
         }
    }
    if ( r2 == 0 ) {
        KStream * response = NULL;
        r2 = KHttpResultGetInputStream ( rslt, & response );
        char base [ 512 ] = "";
        size_t num_read = 0;
        OUTMSG ( ( "KStreamRead(KHttpResult) = " ) );
        r2 = KStreamRead ( response, base, sizeof base, & num_read );
        if ( r2 != 0 )
            OUTMSG ( ( "FAILURE: %R\n", r2 ) );
        else {
            if ( num_read >= sizeof base )
                base [ sizeof base - 1 ] = '\0';
            else {
                base [ num_read ] = '\0';
                while ( num_read > 0 ) {
                    -- num_read;
                    if ( base [ num_read ] == '\n' )
                        base [ num_read ] = '\0';
                    else
                        break;
                }
            }
            OUTMSG ( ( "%s: OK\n", base ) );
        }
        KStreamRelease ( response );
        response = NULL;
    }
    KHttpResultRelease ( rslt );
    rslt = NULL;
    KHttpRequestRelease ( req );
    req = NULL;
    if ( r2 == 0 )
        OUTMSG ( ( "< 1.1 Access to ftp-trace: OK\n" ) );
    else
        OUTMSG ( ( "< 1.1 Access to ftp-trace: FAILURE\n" ) );
    if ( r1 == 0 )
        r1 = r2;
  }/* 100000000
SRR042846 !!!
0 13769 => SRR053325 Twice weekly longitudinal vaginal sampling
1 13997 => SRR045450 2010 Human Microbiome Project filtered data is public 
srapath SRR015685 200000000
      0 30461 => SRR013401 2011 proteobacteria Brucella genome shotg. sequencing
srapath SRR010945
      0 31609 => SRR002749 2008 454 sequencing of Human immunodeficiency virus 1
srapath SRR002682
      0 57513 => SRR000221 Nitrosopumilus maritimus SCM1 FGBT genomic fragment*/
  if(0){
    String domain;
    uint16_t port = 443;
    CONST_STRING ( & domain, "www.ncbi.nlm.nih.gov" );
    OUTMSG ( ( "> 1.2 Checking Access to %S...\n", & domain ) );
    KEndPoint ep;
    OUTMSG ( ( "KNSManagerInitDNSEndpoint(%S:%hu) = ", & domain, port ) );
    rc_t r2 = KNSManagerInitDNSEndpoint ( mgr, & ep,
                                          & domain, port );
    if ( r2 == 0 ) {
        char endpoint [ 1024 ] = "";
        rc_t rx = endpoint_to_string ( endpoint, sizeof endpoint, & ep );
        if ( rx != 0 )
            OUTMSG ( ( "CANNOT CONVERT TO STRING\n" ) );
        else
            OUTMSG ( ( "'%s': OK\n", endpoint ) );
    }
    else
        OUTMSG ( ( "FAILURE: %R\n", r2 ) );
    if ( r2 == 0 )
        OUTMSG ( ( "< 1.2 Access to %S: OK\n", & domain ) );
    else
        OUTMSG ( ( "< 1.2 Access to %S: FAILURE\n", & domain ) );
    if ( r1 == 0 )
        r1 = r2;
  }
  if(0){
    //const char domain [] = "sra-download.ncbi.nlm.nih.gov";
    const char domain [] = "ftp-trace.ncbi.nlm.nih.gov";
    //const char domain [] = "ftp32-1.st-va.ncbi.nlm.nih.gov";
    //const char domain [] = "ftp12.be-md.ncbi.nlm.nih.gov";
    const char dir [] = "/sra/refseq";
    OUTMSG ( ( "> 1.3 Checking Access to %s...\n", domain ) );
    ver_t http_vers = 0x01010000;
    const char proto [] = "https://";
    const char file [] = "/KC702174.1";
    //const char path [] = "/srapub/SRR002749";
    char path [ PATH_MAX ] = "";
    KClientHttpRequest * kns_req = NULL;
    rc_t r2 = string_printf ( path, sizeof path, NULL, "%s%s%s%s",
                         proto, domain, dir, file );
    uint64_t sz = 0;
    if ( r2 != 0 )
        OUTMSG ( ( "CANNOT PRINT PATH: %R\n", r2 ) );
    else {
        const KFile * file = NULL;
        OUTMSG ( ( "KFile = KNSManagerMakeReliableHttpFile(%s)... ", path ) );
        r2 = KNSManagerMakeReliableHttpFile ( mgr, & file, NULL,
                                              http_vers, path);
        if ( r2 == 0 )
            r2 = KFileSize ( file, & sz );
        if ( r2 == 0 )
            OUTMSG ( ( " Size(KFile)=%lu: OK\n", sz ) );
        else
            OUTMSG ( ( "FAILURE: %R\n", r2 ) );
        KFileRelease ( file );
        file = NULL;
    }
    if ( r2 == 0 ) {
        OUTMSG ( ( "KHttpRequest = KNSManagerMakeRequest(%s)... ", path ) );
        r2 = KNSManagerMakeRequest ( mgr,
            & kns_req, http_vers, NULL, "%s", path );
        if ( r2 == 0 )
            OUTMSG ( ( "OK\n" ) );
        else
            OUTMSG ( ( "FAILURE: %R\n", r2 ) );
    }
    KHttpResult * rslt = NULL;
    if ( r2 == 0 ) {
        OUTMSG ( ( "KHttpResult = KHttpRequestGET(KHttpRequest)... " ) );
        r2 = KHttpRequestGET ( kns_req, & rslt );
        if ( r2 == 0 )
            OUTMSG ( ( "OK\n" ) );
        else
            OUTMSG ( ( "FAILURE: %R\n", r2 ) );
    }
    KStream * s = NULL;
    if ( r2 == 0 )
        r2 = KClientHttpResultGetInputStream ( rslt, & s );
    size_t total = 0;
    while ( r2 == 0 ) {
        char buffer [ 1024 ];
        size_t num_read = 0;
        r2 = KStreamRead ( s, buffer, sizeof buffer, & num_read );
        if ( r2 != 0 ) {
            OUTMSG ( ( "KStreamRead(KClientHttpResult)... FAILURE: %R\n",
                        r2 ) );
        }
        else if ( num_read != 0 )
            total += num_read;
        else {
            if ( total == sz )
                OUTMSG ( ( "KStreamRead(KClientHttpResult)... OK\n" ) );
            else
                OUTMSG ( (
                    "KStreamRead(KClientHttpResult) = %zu: SIZE DO NOT MATCH\n",
                    total ) );
            break;
        }
    }
    KStreamRelease ( s );
    s = NULL;
    KClientHttpResultRelease ( rslt );
    rslt = NULL;
    KClientHttpRequestRelease ( kns_req );
    kns_req = NULL;
    if ( r2 == 0 )
        OUTMSG ( ( "< 1.3 Access to %s: OK\n", domain ) );
    else
        OUTMSG ( ( "< 1.3 Access to %s: FAILURE\n", domain ) );
    if ( r1 == 0 )
        r1 = r2;
  }
    //if ( r1 == 0 )        OUTMSG ( ( "< 1 Network: OK\n" ) );    else        OUTMSG ( ( "< 1 Network: FAILURE\n" ) );
    /*OUTMSG ( ( ">>> Checking configuration:\n" ) );
//  OUTMSG ( ( "  Site repository: " ) );    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( ">> Checking remote repository: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "> Checking main/resolver-cgi node: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "> Checking protected/resolver-cgi node: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "> Checking /repository/remote/disabled node: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "> Checking that caching is enabled: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "/repository/user/cache-disabled" ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( ">> Checking main user repository: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "> Checking /user/main/root node: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "> Checking user/main/apps/sra node: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "> Checking user/main/apps/refseq node: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "> Checking user/main/apps/wgs node: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "> Checking user/main/apps/file node: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "> Checking user/main/apps/nannot node: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "> Checking user/main/apps/nakmer node: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "  Protected repositories: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "Network:\n" ) );
    OUTMSG ( ( "  Proxy: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "  Access to NCBI toolkit version file: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "    InitDNSEndpoint: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "    Toolkit version: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "  Toolkit version CGI: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "  Resolver CGI: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "  Access to SRA download site: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "  HTTP download test: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "  Ascp download test: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "  Site resolving: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "  Remote resolving: " ) );
    OUTMSG ( ( "\n" ) );
    OUTMSG ( ( "  Cache resolving: " ) );
    OUTMSG ( ( "\n" ) );*/
/*    rc_t rc2 = 0;
    const char acc[] = "SRR000001";
    const char path[] = "/repository/remote/protected/CGI/resolver-cgi";
    const KConfigNode *node = NULL;
    assert(self);
    rc = KConfigOpenNodeRead(self->cfg, &node, "%s", path);
    if (rc == 0) {
        OUTMSG(("configuration: found\n"));
    }
    else {
        OUTMSG(("ERROR: configuration not found or incomplete\n"));
    }
    if (rc == 0) {
        rc_t rc3 = MainCallCgi(self, node, acc);
        if (rc3 != 0 && rc2 == 0) {
            rc2 = rc3;
        }
        rc3 = MainQuickResolveQuery(self, acc);
        if (rc3 != 0 && rc2 == 0) {
            rc2 = rc3;
        }
    }
    RELEASE(KConfigNode, node);
    if (rc2 != 0 && rc == 0) {
        rc = rc2;
    }
        */
    return rc;
}
