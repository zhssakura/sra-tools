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

/* %{
   Prologue
   %}
   Declarations
   %%
   Grammar rules
   %%
   Epilogue
   */

%{
    #include <stdio.h>
    #include <ctype.h>
    #include <stdlib.h>
    #include <string.h>
    #include <sys/stat.h>
    #include <sys/types.h>
    #include <regex.h>
    #include <stdint.h>

    #include "samextract.h"
    #include "samextract-lib.h"
    #include "samextract-tokens.h"

    #define YYDEBUG 1
/*    #define SAMdebug 1 */

    size_t alignfields=2; // 1 based, QNAME is #1

    int SAMerror(const char * s)
    {
        ERR("%s",s);
        return 0;
    }

    // Returns 1 if match found
    int regexcheck(const char *regex, const char * value)
    {
        regex_t preg;

        int result=regcomp(&preg, regex, REG_EXTENDED);
        if (result)
        {
            size_t s=regerror(result, &preg, NULL, 0);
            char *errmsg=malloc(s);
            regerror(result, &preg, errmsg, s);
            ERR("regcomp error on '%s': %s", regex, errmsg);
            free(errmsg);
            regfree(&preg);
            return 0;
        }

        regmatch_t matches[1];
        if (regexec(&preg, value, 1, matches, 0))
        {
            ERR("Value: '%s' doesn't match regex '%s'", value, regex);
            regfree(&preg);
            return 0;
        }
        regfree(&preg);
        return 1;
    }

    // Returns 1 if OK
    int validate(const char * tag, const char * value)
    {
        /* Pair of TAG, regexp: "/..." TODO: or integer range "1-12345" */
        const char * validations[] =
        {
            "VN", "/.*", // @PG also has "/[0-9]+\\.[0-9]+",
            "SO", "/unknown|unsorted|queryname|coordinate",
            "GO", "/none|query|reference",
            "SN", "/[!-)+-<>-~][!-~]*",
            "LN", "/[0]*[1-9][0-9]{0,10}", // TODO: range check 1..2**31-1
            "AS", "/.*",
            "MD", "/[0-9A-Z\\*]{32}", // bam.c treats same as M5
            "M5", "/[0-9A-Za-z\\*]{32}", // TODO: lowercase acceptable?
            "SP", "/.*",
            "UR", "/.*",
            "ID", "/.*",
            "CN", "/.*",
            "DS", "/.*",
            "DT", "/.*",
            "FO", "/\\*|[ACMGRSVTWYHKDBN]+",
            "KS", "/.*",
            "LB", "/.*",
            "PG", "/.*",
            "PI", "/.*",
            "PL", "/.*",
            "PM", "/.*",
            "PU", "/.*",
            "SM", "/.*",
            "PN", "/.*",
            "CL", "/.*",
            "PP", "/.*",
            "DS", "/.*",
            "\0", "\0"
        };

        int ok=0;

        for (size_t i=0;;++i)
        {
            const char *valtag=validations[i*2];
            const char *valval=validations[i*2+1];
            if (*valtag=='\0')
            {
                WARN("No validation for tag %s", tag);
                ok=1;
                break;
            }
            if (!strcmp(tag, valtag))
            {
                if (valval[0]=='/')
                {
                    ok=regexcheck(valval+1, value);
                    break;
                } else
                {
                // Parse integer range
                    WARN("range not implemented");
                    ok=1;
                }
            }
        }

        return ok;
    }

    void check_required_tag(const char * tags, const char * tag)
    {
        if (!strstr(tags,tag))
        {
            ERR("%s tag not seen in header", tag);
        }
    }

    // Returns 1 if OK
    int checkopttagtype(const char * optfield)
    {
        const char *opttypes="AMi ASi BCZ BQZ CCZ CMi COZ CPi CQZ CSZ CTZ E2Z FIi FSZ FZZ H0i H1i H2i HIi IHi LBZ MCZ MDZ MQi NHi NMi OCZ OPi OQZ PGZ PQi PTZ PUZ QTZ Q2Z R2Z RGZ RTZ SAZ SMi TCi U2Z UQi";
        const char type=optfield[3];
        char tag[3];

        tag[0]=optfield[0];
        tag[1]=optfield[1];
        tag[2]='\0';

        if (tag[0]=='X' ||
            tag[0]=='Y' ||
            tag[0]=='Z') return 1;

        const char *p=strstr(opttypes,tag);
        if (p==NULL) return 1;

        if (p[2]!=type)
        {
            ERR("tag %s should have type %c, not %c", tag, p[2], type); return 0;
        }

        return 1;
    }

    void mark_headers(const char * type)
    {
        for (u32 i=0; i!=VectorLength(&globstate->headers); ++i)
        {
            Header * hdr;
            hdr=VectorGet(&globstate->headers,i);
            if (!strcmp(hdr->headercode,"TBD"))
                hdr->headercode=strdup(type);
//            DBG("Vector header[%d]: %s %s %s", i, hdr->headercode,hdr->tag,hdr->value);
        }
    }

    void process_tagvalue(const char * tag, const char * value)
    {
        if (strlen(tag)!=2)
        {
            ERR("tag '%s' must be 2 characters", tag);
        }

        if (islower(tag[0] &&
            islower(tag[1])))
        {
            DBG("optional tag");
        } else
        {
            if (!validate(tag, value))
            {
                ERR("Tag validataion %s failed",tag);
            }
            globstate->tags=realloc(globstate->tags, strlen(globstate->tags) + strlen(tag) + 1 + 1);
            strcat(globstate->tags,tag); strcat(globstate->tags," ");

            if (!strcmp(tag,"SN"))
            {
                char * s=malloc(strlen(value)+2);
                strcpy(s,value);
                strcat(s," ");
                if (strstr(globstate->seqnames,s))
                {
                    ERR("duplicate sequence %s", value);
                }
                globstate->seqnames=realloc(globstate->seqnames,strlen(globstate->seqnames) + strlen(value) + 1 + 1);
                strcat(globstate->seqnames,s);
                free(s);
            }
            if (!strcmp(tag,"ID"))
            {
                char * s=malloc(strlen(value)+2);
                strcpy(s,value);
                strcat(s," ");
                if (strstr(globstate->ids,s))
                {
                    ERR("duplicate id %s", value);
                }
                globstate->ids=realloc(globstate->ids,strlen(globstate->ids) + strlen(value) + 1 + 1);
                strcat(globstate->ids,s);
                free(s);
            }
        }
        Header * hdr=calloc(1,sizeof(Header));
        hdr->headercode="TBD";
        hdr->tag=strdup(tag);
        hdr->value=strdup(value);
        VectorAppend(&globstate->headers,NULL,hdr);
    }

    void process_align(const char *field)
    {
        const char * opt="(required)";
        if (alignfields>=12) opt="(optional)";
        DBG("alignvalue #%zu%s: %s", alignfields, opt, field);
        switch (alignfields)
        {
            case 2: // FLAG
            {
                int flag;
                if (sscanf(field, "%d", &flag)!=1 ||
                    flag < 0 ||
                    flag > 4095)
                {
                    ERR("error parsing FLAG: %s", field);
                }
                DBG("flag is %d",flag);
                break;
            }
            case 3: // RNAME
            {
                const char * rname=field;
                if (!regexcheck("\\*|[!-)+-<>-~][!-~]*",rname))
                {
                    ERR("error parsing RNAME");
                }
                DBG("rname is %s",rname);
                globstate->rname=strdup(rname);
                break;
            }
            case 4: // POS
            {
                int pos;
                if (sscanf(field, "%d", &pos)!=1 ||
                    pos < 0 ||
                    pos > INT32_MAX)
                {
                    ERR("error parsing POS: %s", field);
                }
                DBG("pos is %d",pos);
                globstate->pos=pos;
                break;
            }
            case 5: // MAPQ
            {
                int mapq;
                if (sscanf(field, "%d", &mapq)!=1 ||
                    mapq < 0 ||
                    mapq > UINT8_MAX)
                {
                    ERR("error parsing MAPQ: %s", field);
                }
                DBG("mapq is %d", mapq);
                break;
            }
            case 6: // CIGAR
            {
                const char * cigar=field;
                if (!regexcheck("\\*|([0-9]+[MIDNSHPX=])+",cigar))
                {
                    ERR("error parsing cigar");
                }
                DBG("cigar is %s",cigar);
                globstate->cigar=strdup(cigar);
                break;
            }
            case 7: // RNEXT
            {
                const char * rnext=field;
                if (!regexcheck("\\*|=|[!-)+-<>-~][!-~]*",rnext))
                {
                    ERR("error parsing rnext");
                }
                DBG("rnext is %s",rnext);
                break;
            }
            case 8: // PNEXT
            {
                int pnext;
                if (sscanf(field, "%d", &pnext)!=1 ||
                    pnext < 0 ||
                    pnext > INT32_MAX)
                {
                    ERR("error parsing PNEXT: %s", field);
                }
                DBG("pnext is %d",pnext);
                break;
            }
            case 9: // TLEN
            {
                int tlen;
                if (sscanf(field, "%d", &tlen)!=1 ||
                    tlen < INT32_MIN ||
                    tlen > INT32_MAX)
                {
                    ERR("error parsing TLEN: %s", field);
                }
                DBG("tlen is %d", tlen);
                break;
            }
            case 10: // SEQ
            {
                const char * seq=field;
                if (!regexcheck("\\*|[A-Za-z=.]+",seq))
                {
                    ERR("error parsing seq");
                }
                DBG("seq is %s",seq);
                globstate->read=strdup(seq);
                break;
            }
            case 11: // QUAL
            {
                const char * qual=field;
                if (!regexcheck("[!-~]+",qual))
                {
                    ERR("error parsing qual");
                }
                DBG("qual is %s", qual);
                break;
            }
            default: // Optional
            {
 //               /TT:t:
                if ((strlen(field)<5) ||
                  field[2]!=':' ||
                  field[4]!=':')
                  {
                    ERR("invald tagtypevalue:%s", field);
                  }
                const char type=field[3];
                if (!checkopttagtype(field))
                {
                    WARN("Optional field tag %s doesn't match type", field);

                }
                const char * value=&field[5];
                switch (type)
                {
                    case 'A':
                        if (!regexcheck("[!-~]", value))
                            ERR("value doesn't match A type:%s",value);
                        break;
                    case 'i':
                        if (!regexcheck("[-+]?[0-9]+", value))
                            ERR("value doesn't match i type:%s",value);
                        break;
                    case 'f':
                        if (!regexcheck("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", value))
                            ERR("value doesn't match f type:%s",value);
                        break;
                    case 'Z':
                        if (!regexcheck("[ !-~]*", value))
                            ERR("value doesn't match Z type:%s",value);
                        break;
                    case 'H':
                        if (!regexcheck("([0-9A-F][0-9A-F])*", value))
                            ERR("value doesn't match H type:%s",value);
                        break;
                    case 'B':
                        if
                        (!regexcheck("[cCsSiIf](,[-+]?[0-9]*\\.?[0-9]+(eE][-+]?[0-9]+)?)+", value))
                            ERR("value doesn't match B type:%s",value);
                        break;
                    default:
                        break;
                }
                DBG("optional field:%s", field);
                break;
            }
        }
        ++alignfields;

    }

%}

/* Declarations */
%union {
 int intval;
 char * strval;
 double floatval;
}

%token <strval> HEADER
%token <strval> SEQUENCE
%token <strval> READGROUP
%token <strval> PROGRAM
%token <strval> COMMENT
%token <strval> TAG
 /* %token <strval> BADTAG */
%token <strval> VALUE
%token <strval> ALIGNVALUE
/* %token DIGITS */
%token <strval> QNAME
%token COLON
%token TAB
%token CONTROLCHAR
%token EOL
%token END 0 "end of file"
%error-verbose
%name-prefix="SAM"
/* %define api.pure // was pure-parser */
/* %lex-param   { extractor * state } */
/* %parse-param { extractor * state } */
%require "2.5" //TODO: 3.0.4+
%error-verbose


%%
 /* Grammar rules */
sam: /* beginning of input */
   /* empty %empty */
   | sam line
   ;

line:
   EOL /* Spec is unclear about empty lines, accept for now */
   | CONTROLCHAR { ERR("CONTROLCHAR"); }
   | comment { DBG("comment"); }
   | header EOL { DBG("header"); }
   | sequence { DBG("sequence"); }
   | program { DBG("program"); }
   | readgroup { DBG("readgroup"); }
   | alignment { DBG("alignment"); }
   ;

comment:
    COMMENT { 
        mark_headers("CO");
    }
    ;

header:
    HEADER tagvaluelist
    {
        DBG("header tagvaluelist");
        check_required_tag(globstate->tags,"VN");
        if (!strcmp(globstate->tags,"SO ") &&
            !strcmp(globstate->tags,"GO "))
           WARN("Both SO and GO tags present");
        if (!(strcmp(globstate->tags,"SO ") ||
              strcmp(globstate->tags,"GO ")))
           WARN("neither SO or GO tags present");
        free(globstate->tags);
        globstate->tags=strdup("");

        mark_headers("HD");
    }
    ;

sequence:
    SEQUENCE tagvaluelist
    {
        DBG("sequence");
        DBG(" sequences were: %s", globstate->seqnames);
        check_required_tag(globstate->tags,"SN");
        check_required_tag(globstate->tags,"LN");
        free(globstate->tags);
        globstate->tags=strdup("");
        mark_headers("SQ");
    }
    ;

program:
     PROGRAM tagvaluelist
     {
        DBG("ids were: %s", globstate->ids);
        DBG("program");
        check_required_tag(globstate->tags,"ID");
        free(globstate->tags);
        globstate->tags=strdup("");
        mark_headers("PG");
     }
     ;


readgroup:
     READGROUP tagvaluelist
     {
        DBG("readgroup");
        DBG("ids were: %s", globstate->ids);
        check_required_tag(globstate->tags,"ID");
        free(globstate->tags);
        globstate->tags=strdup("");
        mark_headers("RG");
     }
     ;

tagvaluelist: tagvalue { DBG(" one tagvaluelist"); }
  | tagvaluelist tagvalue { DBG(" many tagvaluelist"); }
  ;

tagvalue: TAB TAG COLON VALUE {
        DBG("tagvalue:%s=%s", $2, $4);
        const char * tag=$2;
        const char * value=$4;
        process_tagvalue(tag,value);
        free($2);
        free($4);
        };
  | TAB TAB TAG COLON VALUE { ERR("two tabs"); }
  | TAB TAB EOL { ERR("empty tags"); }
  | TAB TAG TAG {
        const char * tag=$2;
        WARN("malformed TAG:VALUE 'TAB %s(NOT COLON)...'", tag);
        }
  | TAB EOL { WARN("empty tags"); }
  ;

alignment:
    QNAME avlist
    {
        DBG(" avlist qname:%s fields=%zu", $1, alignfields);
        alignfields=2;
        Alignment * align=calloc(1,sizeof(Alignment));
        align->read=globstate->read;
        align->cigar=globstate->cigar;
        align->rname=globstate->rname;
        align->pos=globstate->pos;
        VectorAppend(&globstate->alignments,NULL,align);
        free($1);
    }
    ;

avlist:
      av { DBG(" one av"); }
 |    avlist av {
           // TODO"bison: many avlist");
            }
    ;

av:
    TAB ALIGNVALUE
    {
        const char * field=$2;
        process_align(field);
        free($2);
    }
    ;

%%


 /* Epilogue */

