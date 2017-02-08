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

    #include "sam.tab.h"

    #define YYDEBUG 1

    extern int SAMlex(void);
    extern int SAMerror(const char * s);
    extern int moredata(char * buf,int * numbytes, int maxbytes);

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
        "M5", "/[0-9A-Z\\*]{32}",
        "SP", "/.*",
        "UR", "/.*:/.*",
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

    // Returns 1 if match found
    int regexcheck(const char *regex, const char * value)
    {
        regex_t preg;

//        fprintf(stderr,"Compiling %s\n", regex);
        int result=regcomp(&preg, regex, REG_EXTENDED);
        if (result)
        {
            size_t s=regerror(result, &preg, NULL, 0);
            char *errmsg=malloc(s);
            regerror(result, &preg, errmsg, s);
            fprintf(stderr,"regcomp error on '%s': %s\n", regex, errmsg);
            free(errmsg);
            regfree(&preg);
            return 0;
        }

        regmatch_t matches[1];
        if (regexec(&preg, value, 1, matches, 0))
        {
            fprintf(stderr,"Value: '%s' doesn't match regex '%s'\n", value, regex);
            regfree(&preg);
            return 0;
        }
        //fprintf(stderr,"match\n");
        regfree(&preg);
        return 1;
    }

    // Returns 1 if OK
    int validate(const char * tag, const char * value)
    {
        int ok=0;

        for (size_t i=0;;++i)
        {
            const char *valtag=validations[i*2];
            const char *valval=validations[i*2+1];
            if (*valtag=='\0')
            {
                fprintf(stderr,"No validation for tag %s\n", tag);
                ok=1;
                break;
            }
            if (!strcmp(tag, valtag))
            {
//                fprintf(stderr,"Checking %s\n", valtag);
                if (valval[0]=='/')
                {
                    ok=regexcheck(valval+1, value);
                    break;
                } else
                {
                // Parse integer range
                    fprintf(stderr,"range not implemented\n");
                    ok=1;
                }
            }
        }

        return ok;
    }

    size_t alignfields=2; // 1 based, QNAME is #1
    extern char *tags;
    extern char *seqnames;
    extern char *ids;

    void check_required_tag(const char * tag)
    {
        if (!strstr(tags,tag))
        {
            fprintf(stderr,"error: %s tag not seen in header\n", tag);
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
            fprintf(stderr,"error: tag %s should have type %c, not %c\n", tag, p[2], type); return 0;
        }

        return 1;
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
%require "2.5"


%%
 /* Grammar rules */
sam: /* beginning of input */
   /* empty %empty */
   | sam line
   ;

line:
   EOL /* Spec is unclear about empty lines, accept for now */
   | CONTROLCHAR { fprintf(stderr,"error: CONTROLCHAR\n"); }
   | comment { fprintf(stderr,"comment\n"); }
   | header EOL { fprintf(stderr,"header\n\n"); }
   | sequence {fprintf(stderr,"sequence\n\n"); }
   | program {fprintf(stderr,"program\n\n"); }
   | readgroup {fprintf(stderr,"readgroup\n\n"); }
   | alignment { fprintf(stderr,"alignment\n\n"); }

comment:
    COMMENT { }

header:
    HEADER tagvaluelist
    {
        fprintf(stderr,"header tagvaluelist\n");
        check_required_tag("VN");
        if (!strcmp(tags,"SO ") &&
            !strcmp(tags,"GO "))
           fprintf(stderr,"warn: Both SO and GO tags present\n");
        if (!(strcmp(tags,"SO ") ||
              strcmp(tags,"GO ")))
           fprintf(stderr,"warn: neither SO or GO tags present\n");
        free(tags);
        tags=strdup("");
    }

sequence:
    SEQUENCE tagvaluelist
    {
        fprintf(stderr, "sequence\n");
        fprintf(stderr," sequences were: %s\n", seqnames);
        check_required_tag("SN");
        check_required_tag("LN");
        free(tags);
        tags=strdup("");
        }

program:
     PROGRAM tagvaluelist
     {
        fprintf(stderr,"ids were: %s\n", ids);
        fprintf(stderr, "program\n");
        check_required_tag("ID");
        free(tags);
        tags=strdup("");
     }


readgroup:
     READGROUP tagvaluelist
     {
        fprintf(stderr, "readgroup\n");
        fprintf(stderr,"ids were: %s\n", ids);
        check_required_tag("ID");
        free(tags);
        tags=strdup("");
     }

tagvaluelist: tagvalue { fprintf(stderr, " one tagvaluelist\n"); }
  | tagvaluelist tagvalue { fprintf(stderr, " many tagvaluelist\n"); }
  ;

tagvalue: TAB TAG COLON VALUE {
        fprintf(stderr,"tagvalue:%s=%s\n", $2, $4);
        const char * tag=$2;
        const char * value=$4;

        if (strlen(tag)!=2)
        {
            fprintf(stderr,"tag '%s' must be 2 characters\n", tag);
        }

        if (islower(tag[0] &&
            islower(tag[1])))
        {
            fprintf(stderr,"optional tag\n");
        } else
        {
            validate(tag, value);
            tags=realloc(tags, strlen(tags) + strlen(tag) + 1 + 1);
            strcat(tags,tag); strcat(tags," ");

            if (!strcmp(tag,"SN"))
            {
                char * s=malloc(strlen(value)+2);
                strcpy(s,value);
                strcat(s," ");
                if (strstr(seqnames,s))
                {
                    fprintf(stderr,"error: duplicate sequence %s\n", value);
                }
                seqnames=realloc(seqnames,strlen(seqnames) + strlen(value) + 1 + 1);
                strcat(seqnames,s);
                free(s);
            }
            if (!strcmp(tag,"ID"))
            {
                char * s=malloc(strlen(value)+2);
                strcpy(s,value);
                strcat(s," ");
                if (strstr(ids,s))
                {
                    fprintf(stderr,"error: duplicate id %s\n", value);
                }
                ids=realloc(ids,strlen(ids) + strlen(value) + 1 + 1);
                strcat(ids,s);
                free(s);
            }
        }
        free($2);
        free($4);
        };
  | TAB TAB TAG COLON VALUE { fprintf(stderr,"two tabs\n"); }
  | TAB TAB EOL { fprintf(stderr,"empty tags\n"); }
  | TAB TAG TAG {
        const char * tag=$2;
        fprintf(stderr,"error: warning: malformed TAG:VALUE 'TAB %s(NOT COLON)...'\n", tag);
        }
  | TAB EOL { fprintf(stderr,"empty tags\n"); }

  ;

alignment:
    QNAME avlist
    {
        fprintf(stderr," avlist qname:%s fields=%zu\n", $1, alignfields);
        alignfields=2;
        free($1);
    }

avlist:
      av { fprintf(stderr," one av\n"); }
 |    avlist av {
           // fprintf(stderr,"bison: many avlist\n");
            }

av:
    TAB ALIGNVALUE
    {
        const char * field=$2;
        const char * opt="(required)";
        if (alignfields>=12) opt="(optional)";
        fprintf(stderr,"alignvalue #%zu%s: %s\n", alignfields, opt, field);
        switch (alignfields)
        {
            case 2: // FLAG
            {
                int flag;
                if (sscanf(field, "%d", &flag)!=1 ||
                    flag < 0 ||
                    flag > 4095)
                {
                    fprintf(stderr,"error parsing FLAG: %s\n", field);
                }
                fprintf(stderr,"flag is %d\n",flag);
                break;
            }
            case 3: // RNAME
            {
                const char * rname=field;
                if (!regexcheck("\\*|[!-)+-<>-~][!-~]*",rname))
                {
                    fprintf(stderr,"error parsing RNAME\n");
                }
                fprintf(stderr,"rname is %s\n",rname);
                break;
            }
            case 4: // POS
            {
                int pos;
                if (sscanf(field, "%d", &pos)!=1 ||
                    pos < 0 ||
                    pos > INT32_MAX)
                {
                    fprintf(stderr,"error parsing POS: %s\n", field);
                }
                fprintf(stderr,"pos is %d\n",pos);
                break;
            }
            case 5: // MAPQ
            {
                int mapq;
                if (sscanf(field, "%d", &mapq)!=1 ||
                    mapq < 0 ||
                    mapq > UINT8_MAX)
                {
                    fprintf(stderr,"error parsing MAPQ: %s\n", field);
                }
                fprintf(stderr,"mapq is %d\n", mapq);
                break;
            }
            case 6: // CIGAR
            {
                const char * cigar=field;
                if (!regexcheck("\\*|([0-9]+[MIDNSHPX=])+",cigar))
                {
                    fprintf(stderr,"error parsing cigar\n");
                }
                fprintf(stderr,"cigar is %s\n",cigar);
                break;
            }
            case 7: // RNEXT
            {
                const char * rnext=field;
                if (!regexcheck("\\*|=|[!-)+-<>-~][!-~]*",rnext))
                {
                    fprintf(stderr,"error parsing rnext\n");
                }
                fprintf(stderr,"rnext is %s\n",rnext);
                break;
            }
            case 8: // PNEXT
            {
                int pnext;
                if (sscanf(field, "%d", &pnext)!=1 ||
                    pnext < 0 ||
                    pnext > INT32_MAX)
                {
                    fprintf(stderr,"error parsing PNEXT: %s\n", field);
                }
                fprintf(stderr,"pnext is %d\n",pnext);
                break;
            }
            case 9: // TLEN
            {
                int tlen;
                if (sscanf(field, "%d", &tlen)!=1 ||
                    tlen < INT32_MIN ||
                    tlen > INT32_MAX)
                {
                    fprintf(stderr,"error parsing TLEN: %s\n", field);
                }
                fprintf(stderr,"tlen is %d\n", tlen);
                break;
            }
            case 10: // SEQ
            {
                const char * seq=field;
                if (!regexcheck("\\*|[A-Za-z=.]+",seq))
                {
                    fprintf(stderr,"error parsing seq\n");
                }
                fprintf(stderr,"seq is %s\n",seq);
                break;
            }
            case 11: // QUAL
            {
                const char * qual=field;
                if (!regexcheck("[!-~]+",qual))
                {
                    fprintf(stderr,"error parsing qual\n");
                }
                fprintf(stderr,"qual is %s\n", qual);
                break;
            }
            default: // Optional
            {
 //               /TT:t:
                if ((strlen(field)<5) ||
                  field[2]!=':' ||
                  field[4]!=':')
                  {
                    fprintf(stderr,"error: invald tagtypevalue:%s\n", field);
                  }
                const char type=field[3];
                if (!checkopttagtype(field))
                {
                    fprintf(stderr,"Optional field tag %s doesn't match type\n", field);

                }
                const char * value=&field[5];
                switch (type)
                {
                    case 'A':
                        if (!regexcheck("[!-~]", value))
                            fprintf(stderr,"error: value doesn't match A type:%s\n",value);
                        break;
                    case 'i':
                        if (!regexcheck("[-+]?[0-9]+", value))
                            fprintf(stderr,"error: value doesn't match i type:%s\n",value);
                        break;
                    case 'f':
                        if (!regexcheck("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", value))
                            fprintf(stderr,"error: value doesn't match f type:%s\n",value);
                        break;
                    case 'Z':
                        if (!regexcheck("[ !-~]*", value))
                            fprintf(stderr,"error: value doesn't match Z type:%s\n",value);
                        break;
                    case 'H':
                        if (!regexcheck("([0-9A-F][0-9A-F])*", value))
                            fprintf(stderr,"error: value doesn't match H type:%s\n",value);
                        break;
                    case 'B':
                        if
                        (!regexcheck("[cCsSiIf](,[-+]?[0-9]*\\.?[0-9]+(eE][-+]?[0-9]+)?)+", value))
                            fprintf(stderr,"error: value doesn't match B type:%s\n",value);
                        break;
                    default:
                        break;
                }
                fprintf(stderr,"optional field:%s\n", field);
                break;
            }
        }
        ++alignfields;
        free($2);
    }

%%


 /* Epilogue */

