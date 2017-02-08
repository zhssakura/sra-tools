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

   /*
      lex file is three sections:
      definitions: name definition
      %%
      rules: pattern action
      %%
      user code
   */

  /*
    @HD\tVN:1.4\tSO:coordinate
    @SQ\tSN:1\t...
    @SQ\tSN:16\t...
    @RG\tID:PM89\tPL:Illumina\tLB...
    HWI-...
  */

%top{
#include <stdio.h>
#include <string.h>
#include "sam.tab.h"

//#define YYSTYPE char
#define LEXDEBUG 0

#undef YY_INPUT
#define YY_INPUT(buf,numbytes,maxbytes) moredata(buf,&numbytes,maxbytes)
extern int moredata(char * buf,int * numbytes, int maxbytes);


    // do { if (LEXDEBUG) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__); } while (0)
#define lex_print(fmt, ...) \
    do { if (LEXDEBUG) fprintf(stderr, fmt, ##__VA_ARGS__); } while (0)

  }

/* definitions */
digits [[:digit:]]+
integer [-+]?[0-9]+
float  [-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?
eol    (\r\n|\r|\n)
controlchar [\x00-\x08\x0b-\x0c\x0e-\x1f]

tag [A-Za-z][A-Za-z0-9]

/* No @ in first, to distinguish from headers */
qname [!-?A-~][!-~]{0,253}

/* stop after first EOF */
%option prefix="SAM"
%option warn
%option noyywrap
%option nodefault
%option nounput
%option nodebug
/* Not a TTY */
%option never-interactive
%option yylineno
%option nounistd
/* Allow multitheading */
/* %option reentrant */
/* yylex takes yylval */
/* %option bison-bridge */
%x INALIGNMENT
%x AFTERCOLON

%%
 /* rules */
^@HD   { lex_print( "lex: Header\n"); return HEADER; }
^@SQ   { lex_print( "lex: SQ\n"); return SEQUENCE; }
^@RG   { lex_print( "lex: RG\n"); return READGROUP; }
^@PG   { lex_print( "lex: PG\n"); return PROGRAM; }
^@CO\t.*{eol} { lex_print( "lex: Comment\n"); return COMMENT; }

 /* TODO: replace strdups with pool/chunk allocator */
<INITIAL,INALIGNMENT>^{qname} { BEGIN INALIGNMENT; SAMlval.strval=strdup(yytext); lex_print("\nlex: alignment qname, INALIGNMENT\n"); return QNAME; }
<INITIAL,INALIGNMENT,AFTERCOLON>{controlchar} { lex_print("lex:CONTROLCHAR1\n"); return CONTROLCHAR; }
<INITIAL,INALIGNMENT>\t { lex_print("lex: tab\n"); return TAB;}
<AFTERCOLON>\t { BEGIN INITIAL; lex_print("lex: tab INITIAL\n"); return TAB;}

<INITIAL,INALIGNMENT,AFTERCOLON>{eol} { BEGIN INITIAL; lex_print("lex: eol INITIAL\n"); return EOL;}
<INALIGNMENT>[ -~]+ { SAMlval.strval=strdup(yytext); lex_print("lex: alignvalue\n"); return ALIGNVALUE; }
{tag}  { SAMlval.strval=strdup(yytext); lex_print("lex: Valid Tag:%s\n", yytext); return TAG; }
<INITIAL,INALIGNMENT>: { BEGIN AFTERCOLON; lex_print("lex: colon AFTERCOLON\n"); return COLON;}
<AFTERCOLON>[ -~]+ { SAMlval.strval=strdup(yytext); lex_print("lex: Value:%s\n", yytext); return VALUE; }

.|\n       lex_print("DEFAULT '%c' ", *yytext);

%%
 /* user code */
