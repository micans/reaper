

/* features:
      transpose tabular data.
      record and field separator can be changed.
      can read gzipped files.
      fast as long as data fits into memory.
      can handle headers lacking dummy lead field.
*/


/* todo
   -  encapsulate
       - input reading
       - add_fake_tab offset fiddling
       - removed_last_newline logic
   #  account for missing tab in header line .... (count mismatch)
   #  warn for count mismatches / output counts
   -  in check mode do not read in file in entirety.
   /  optionally do not print out trailing nl
   #  gzread
   ?  mmap file
   -  write row, col labels
   /  test with empty lines (rec_offset stress-test)
   /  test with all-empty fields
   -  optify na/nan/-/empty output fields
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>

#include "slib.h"


unsigned char g_tab = '\t';
unsigned char g_nl  = '\n';
int removed_last_newline = 0;

int g_debug = 0, g_check = 0;



struct offset
{  unsigned long this
;  unsigned long next
;  unsigned long n_sep
;
}  ;


int field_check_ok
(  unsigned char* data
,  unsigned long N
,  struct offset* offsets
,  unsigned long n_offsets
)
   {  int ret = 1
   ;  unsigned long i

   ;  for (i=0;i<n_offsets-1;i++)
      {  unsigned char* o = data+offsets[i].this, *p = o, *z = o + offsets[i].next
      ;  while (p < z && (p = memchr(p, g_tab, offsets[i].next - offsets[i].this - (p-o))))
            offsets[i].n_sep++
         ,  p++
      ;  if (!g_check && i == 3)
         break
   ;  }

                           /* note: the fake tab was already added (add_fake_tab), so computations
                            * below look weird: equality means deficit of one,
                            * number of separators == number of fields for leading record.

                            * in check mode we compare with the *second* record.
                           */
      if (g_check)
      {  unsigned long n_nok = 0
      ;  if (n_offsets > 1)
         arrr("second record has %lu fields", (offsets[1].n_sep))
      ;  for (i=2;i<n_offsets-1;i++)
         if (offsets[i].n_sep != offsets[1].n_sep)
            fprintf(stderr, " %lu/%lu", (i+1), offsets[i].n_sep)
         ,  n_nok++
      ;  if (n_nok)
         fputc('\n', stderr)
      ;  else
         arrr("file is clear")
      ;  exit(n_nok ? 1 : 0)
   ;  }

      else if (i >= 2 && offsets[0].n_sep == offsets[1].n_sep)
         arrr
         (  "first and second line have %lu and %lu separators; I'll add a leading separator"
         ,  offsets[0].n_sep-1
         ,  offsets[1].n_sep
         )
      ,  ret = 0

   ;  else if (i >= 2 && offsets[0].n_sep-1 != offsets[1].n_sep)
      arrr("first and second line have %lu and %lu separators; you probably know what your are doing", offsets[0].n_sep-1, offsets[1].n_sep)

   ;  return ret
;  }


void transpose
(  void* fpo
,  unsigned char* data
,  unsigned long N
,  struct offset* offsets
,  unsigned long n_offsets
,  int zippit
)
   {  int n_active
   ;  if (field_check_ok(data, N, offsets, n_offsets))
         offsets[0].this = 1                /* this is tied to add_fake_tab */
      ,  data[1] = g_nl
   ;  do 
      {  unsigned long i
      ;  n_active = 0
      ;  for (i=0; i<n_offsets-1;i++)
         {  unsigned long o = offsets[i].this
         ;  unsigned remaining = offsets[i].next - offsets[i].this
         ;  unsigned char* p = remaining ? memchr(data+o, g_tab, remaining) : NULL
         ;  if (remaining && !p)
            p = data+offsets[i].next
         ;  if (i) {
#if WE_USE_ZLIB
               if (zippit)
               gzputc(fpo, g_tab)
            ;  else
#endif
               fputc(g_tab, fpo)
         ;  }
            if (p)
            {  p[0] = '\0'
#if WE_USE_ZLIB
            ;  if (zippit)
               gzwrite(fpo, data+o+1, p-(data+o)-1)
            ;  else
#endif
               fwrite(data+o+1, 1, p-(data+o)-1,  fpo)
            ;  offsets[i].this = p-data
         ;  }
            if (offsets[i].this < offsets[i].next)
            n_active++
      ;  }
         if (n_active || g_nl == '\n') {
#if WE_USE_ZLIB
            if (zippit)
            gzputc(fpo, g_nl)
         ;  else
#endif
          putc(g_nl, fpo)
      ;  }
      }
      while (n_active > 0)
   ;  if (removed_last_newline) {
#if WE_USE_ZLIB
         if (zippit)
         gzputc(fpo, '\n')
      ;  else
#endif
         putc('\n', fpo)
   ;  }
   }


int main
(  int argc
,  char* argv[]
)
   {  unsigned long n_data_alloc = 1<<20
   ;  unsigned char* data = myalloc(n_data_alloc)  /* must be bigger than BUFSIZE + 1 */
   ;  ZFILE fpo = stdout
   ;  ZFILE input
#define BUFSIZE 4096
   ;  unsigned long N = 0, n_offsets = 0, n_offsets_alloc = 1<<10
   ;  long n = 0
   ;  struct offset *offsets = myalloc(n_offsets_alloc * sizeof offsets[0])
   ;  int i
   ;  int zippit = 1
   ;  const char* g_fnin = "-"
   ;  const char* g_fnout = "-"
   ;

/* enter macromagical option world */
      
      arg_switch()

      optarg("-i")  g_fnin = thearg();  endarg()
      optarg("-o")  g_fnout = thearg();  endarg()
      optarg("-t")  g_tab  = thearg()[0];  endarg()
      optarg("-n")  g_nl   = thearg()[0];  endarg()
      uniarg("--debug")  g_debug = 1; endarg()
#if WE_USE_ZLIB
      uniarg("--nozip")  zippit = 0;  endarg()
#endif
      uniarg("--check")  g_check = 1; endarg()
   uniarg("-h")
puts("-i <fname>        input stream (gzipped file allowed)");
puts("-o <fname>        output stream (gzipped file allowed)");
puts("-t <CHAR>         field delimiter (default TAB)");
puts("-n <CHAR>         record delimiter (default NEWLINE)");
puts("--debug           output record offset information");
puts("--check           check field count consistency");
#if WE_USE_ZLIB
puts("--nozip           do not gzip the output");
#endif
exit(0);
   endarg()

      failarg()
      arg_done()

/* exit macromagicalitaciousness */

   ;  if (!(input = myfopen(g_fnin, "r", 1)))
      exit(1)
   ;  if (!g_check && !(fpo = myfopen(g_fnout, "w", zippit)))
      exit(1)

   ;  offsets[0].this = 0
   ;  offsets[0].next = -1
   ;  n_offsets = 1

   ;  data[N++] = g_nl
   ;  data[N++] = g_tab  /* in case we need to add one when header line is one short: add_fake_tab */

#  if WE_USE_ZLIB
   ;  while ((n = gzread(input, data+N, BUFSIZE)) > 0)
#  else
   ;  while ((n = fread(data+N, 1, BUFSIZE, input)) > 0)
#endif
      {  if (N+n+2+BUFSIZE > n_data_alloc)     /* extra space to write '\0', optionally '\n\0' */
         {  n_data_alloc *= 1.141
         ;  if (n_data_alloc < N+n+2+BUFSIZE)
            n_data_alloc = N+n+2+10*BUFSIZE
         ;  data = myrealloc(data, n_data_alloc)
      ;  }

         data[N+n] = '\0'        /* Should not be necessary actually */

      ;  {  unsigned char* o = data+N, *p = o
         ;  while (p-o < n && (p = memchr(p, g_nl, n-(p-o))))
            {  if (n_offsets+1 > n_offsets_alloc)       /* extra space in case no trailing g_nl */
               {  n_offsets_alloc *= 1.141
               ;  offsets = myrealloc(offsets, n_offsets_alloc * sizeof offsets[0])
            ;  }
               offsets[n_offsets].this = p-data
            ;  offsets[n_offsets-1].next = p-data
            ;  n_offsets++
            ;  p++
         ;  }
         }
         N += n
   ;  }

      myfzclose(input, 1)

  ;   if (g_nl != '\n' && data[N-1] == '\n')
         data[--N] = '\0'
      ,  removed_last_newline = 1

  ;   if (N > 0 && data[N-1] != g_nl)
      {  data[N]  = g_nl
      ;  data[N+1] = '\0'
      ;  offsets[n_offsets-1].next = N
      ;  offsets[n_offsets].this = N
      ;  n_offsets++
      ;  N++
   ;  }

      offsets[n_offsets-1].next = -1

;if(g_debug)
 for (i=0;i<n_offsets;i++)
 arrr("%d record %d - %d", (int) n, (int) offsets[i].this, offsets[i].next)

;if(0)fprintf(stdout, "[%s], %d offsets\n", data, (int) n_offsets)

   ;  transpose(fpo, data, N, offsets, n_offsets, zippit)
   ;  myfzclose(fpo, zippit)
   ;  free(data)
   ;  free(offsets)
   ;  return 0
;  }





