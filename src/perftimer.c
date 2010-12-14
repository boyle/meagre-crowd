/* Meagre-Crowd: A sparse distributed matrix solver testbench for performance benchmarking.
 * Copyright (C) 2010 Alistair Boyle <alistair.js.boyle@gmail.com>
 *
 *     This file is part of Meagre-Crowd.
 *
 *     Meagre-Crowd program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "config.h"
#include <stdlib.h> // malloc/free
#include <errno.h>  // malloc error
#include <err.h>    // malloc error text
#include <string.h> // strnlen, strncpy, etc
#include "perftimer.h"

// local function
// allocate head block
// ph: ptr to ptr, a head block
// out: 0 - success, -1 allocation failed
static inline   int _perftimer_malloc_head(perftimer_t** ph);
static inline   int _perftimer_malloc_head(perftimer_t** ph) {
  if((*ph = malloc(sizeof(perftimer_t))) == NULL) {
    warn("perftimer_malloc(head)"); // shows errno string
    return -1; // allocation failure
  }
  (*ph)->head = NULL;
  (*ph)->tail = NULL;
  (*ph)->old  = NULL;
  (*ph)->current_depth = 0;
  return 0;
}



// perftimer_malloc()
// allocate a perftimer structure
// out: a new perftimer ptr, NULL does NOT indicate failure
//  (memory exhaustion will be indicated when you try to perftimer_inc)
inline perftimer_t* perftimer_malloc() {
  perftimer_t* ptr;
  _perftimer_malloc_head(&ptr); // ignores result code (don't care)
  return ptr;
}

// perftimer_free()
// release entire perftimer structure h
void perftimer_free(perftimer_t* h) {
  if(h == NULL) // nothing to do
    return;

  // deallocate the singly linked list, from the head down
  perftimer_tic_t* p = h->head;
  perftimer_tic_t* p_next;
  while(p != NULL) {
    p_next = p->next;
    free(p->desc); // free the copied string
    free(p);       // and free itself
    p = p_next;    // then move forward one
  }

  // finally, release the head of the structure
  perftimer_t* old = h->old;
  free(h);
  // and follow the 'old' link to release any repititions
  perftimer_free(old); // recursive!
}

// perftimer_inc
// increment the perftimer
// h: a perftimer structure ptr
// s: the description of the most recent event completed
// d: depth level, 0:min
// out: result - 0: success, -2: failure (memory exhaustion), -1: NULL ptr
int perftimer_inc(perftimer_t* h, char const*const s, const size_t n) {
  if(h == NULL)
    return -1; // bad ptr

  // allocate a new block into the selected location
  perftimer_tic_t* ptr;
  if((ptr = malloc(sizeof(perftimer_tic_t))) == NULL) {
    warn("perftimer_inc()"); // shows errno string
    return -2; // malloc failure
  }

  if(h->tail == NULL) { // first in list
    h->head = ptr;
    h->tail = ptr;
  }
  else { // append to list
    h->tail->next = ptr;
    h->tail = ptr;
  }

  // fill in the structure
  if(gettimeofday(&(ptr->now),NULL) != 0)
    warn("failed to store time"); // shouldn't be any reason to fail?
  ptr->depth = h->current_depth;
  ptr->next = NULL;
  { // safe string copy
    // this initial counting is a bit inefficient but worth it for the safety
    size_t nn = strnlen(s,n);
    if(nn == 0) { // no string
      ptr->desc = NULL;
    }
    else { // copy the string in, make sure its null terminated (safety first)
      // TODO use strndup instead?
      ptr->desc = malloc((nn+1)*sizeof(char)); // need an extra char for '\0'
      strncpy(ptr->desc,s,nn);
      ptr->desc[nn] = '\0'; // make sure the string is null terminated
    }
  }

  return 0; // success
}

// utility function: calculate difference between times, returns in seconds
// t2 is assumed to be later than t1
static inline   double   _calc_perftimer_diff(perftimer_tic_t const* const t1, perftimer_tic_t const* const t2);
static inline   double   _calc_perftimer_diff(perftimer_tic_t const* const t1, perftimer_tic_t const* const t2) {
  if((t1 == NULL) || (t2 == NULL))
    return 0.0;

  struct timeval t;
  timersub(&(t2->now), &(t1->now), &t);

  // convert to double
  return ((double) t.tv_sec) + ((double) t.tv_usec) * 1.0e-6;
}

// perftimer_printlen()
// determine the length of the string that would be generated without
// a constraint on length
// h: a perftimer structure ptr
// d: max depth, -1:unlimited, 0:min
// out: the required size of the string, excluding terminating '\0'
size_t perftimer_printlen(perftimer_t const * const h, const unsigned int d) {
  if((h == NULL) || (h->head == NULL))
    return 0;

  size_t c = 0;
  perftimer_tic_t const * ptr = h->head;

  // count desc
  while(ptr != NULL) {
    if((ptr->desc != NULL) && (ptr->depth <= d))
      // we know that all ptr->desc are zero terminated: we copied them
      // need to leave lots of room for calculated time
      c += strlen(ptr->desc) +60;
    // TODO test for overflow on 'c += ...'
    ptr = ptr->next;
  }

  return c+60-1; // extra padding for total
}


// find the longest set of links in the structure
static unsigned int _max_links(perftimer_t const * h);
static unsigned int _max_links(perftimer_t const * h) {
  unsigned int m = 0;
  while(h != NULL) {
    unsigned int c = 0;

    // find count of tic blocks -1
    perftimer_tic_t* p = h->head;
    while(p != NULL) {
      p = p->next;
      if(p != NULL)
        c++;
    }

    if(c > m) // update max if larger
      m = c;

    h = h->old; // try the next list
  }
  return m;
}

// perftimer_snprintf()
// create a string describing the events so far,
// including details upto a depth limit, and restricted to string length n
// h: a perftimer structure ptr
// s: the output string (must be preallocated)
// n: the string length limit
// d: max depth, -1:unlimited
// out: a count of bytes
//
// this code does NOT assume that the different rounds are the same length
// and it only averages across entries that match
int perftimer_snprintf(perftimer_t const * const h, char* s, const size_t n, const unsigned int d) {
  s[0] = '\0'; // string starts empty
  if((h==NULL) || (h->head == NULL)) // nothing to show
    return 0;

  unsigned int m_max = _max_links(h); // longest number of links
  // create list of times
  double* t = calloc(m_max,sizeof(double));
  unsigned int * r = calloc(m_max,sizeof(unsigned int)); // how many went into this count
  perftimer_t const * hs = h;
  perftimer_tic_t const * ptr = h->head;
  while(hs != NULL) {
    ptr = hs->head;
    unsigned int i; // current link
    for(i=0;i<m_max && (ptr != NULL) && (ptr->next != NULL);i++) {
      perftimer_tic_t const * stop = ptr->next;
      while((stop->next != NULL) && (stop->depth > ptr->depth)) {
        stop = stop->next;
      } // summarize lower depths
      t[i] +=  _calc_perftimer_diff(ptr, stop);
      r[i]++;
      ptr = ptr->next;
    }

    hs = hs->old;
  }
  // Note: this could be done without 'r'
  // i.e. t[i] = t[i]*(r-1)/r + new_val/r
  // but this wouldn't be as precise - not sure if it matters
  // but we are dealing with small numbers

  // now add the descriptions, if required

  // count desc
  static const int nw = 60;
  char tmp [nw];
  // start at the top, having already established that there is a first entry
  ptr = h->head;
  int i = 0;
  while(ptr->next != NULL) {
    if((ptr->desc != NULL) && (ptr->depth <= d)) {
      { // indent as appropriate for this depth
        int j;
        for(j=0;j<ptr->depth;j++)
          strcat(s, "  ");
      }
      strncat(s, ptr->desc, n-1);
      if(t[i] > 0.0) {
        if(r[i] > 1)
          snprintf(tmp,nw,":\t%0.3fms (av. %0.3fms in %d rounds)\n",
	    t[i]*1e3,
	    t[i]/((double)r[i])*1e3,
	    r[i]);
	else
          snprintf(tmp,nw,":\t%0.3fms\n",t[i]*1e3);
      }
      else {
        snprintf(tmp,nw,"\n");
      }
      strncat(s, tmp, n-1);
    }
    ptr = ptr->next;
    i++;
  }
  snprintf(tmp,nw,"total:\t%0.3fms",perftimer_wall(h)*1e3);
  strncat(s,tmp,n-1);
  unsigned int rnds = perftimer_rounds(h); // rounds
  if(rnds > 1) {
    snprintf(tmp,nw," (av. %0.3fms in %d rounds)",perftimer_wall_av(h)*1e3,rnds);
    strncat(s,tmp,n-1);
  }
  // make the string safe: zero terminate it
  s[n-1] = '\0'; // Note: man claims strncat appends \0 always, but it doesn't

  free(t); // release time list
  free(r);

  // this should come out to the same size as perftimer_printlen()
  return strnlen(s,n); // total bytes produced
}


// perftimer_snprintf_csv_header()
// create a string describing the events so far,
// including details upto a depth limit, and restricted to string length n
// h: a perftimer structure ptr
// s: the output string (must be preallocated)
// n: the string length limit
// d: max depth, -1:unlimited
// out: a count of bytes
//
// this code does NOT assume that the different rounds are the same length
// and it only averages across entries that match
int perftimer_snprintf_csv_header(perftimer_t const * const h, char* s, const size_t n, const unsigned int d) {
  s[0] = '\0'; // string starts empty
  if((h==NULL) || (h->head == NULL)) // nothing to show
    return 0;

  // count desc
  static const int nw = 60;
  char tmp [nw];
  // start at the top, having already established that there is a first entry
  perftimer_tic_t const * ptr = h->head;
  int i = 0;
  int first = 1;
  while(ptr->next != NULL) {
    if((ptr->desc != NULL) && (ptr->depth <= d)) {
      if(!first) // comma, except first number
        strncat(s, ", ", n-1);
      first = 0;
      strncat(s, ptr->desc, n-1);
    }
    ptr = ptr->next;
    i++;
  }
  if(!first) // comma, except first number
    strncat(s, ", ", n-1);
  first = 0;
  snprintf(tmp,nw,"total (ms)");
  strncat(s,tmp,n-1);
  if(perftimer_rounds(h) > 1) {
    snprintf(tmp,nw,", rounds");
    strncat(s,tmp,n-1);
  }
  // make the string safe: zero terminate it
  s[n-1] = '\0'; // Note: man claims strncat appends \0 always, but it doesn't

  // this is definately not the predicted length
  return strnlen(s,n); // total bytes produced
}


// perftimer_snprintf_csv_body()
// create a string describing the events so far,
// including details upto a depth limit, and restricted to string length n
// h: a perftimer structure ptr
// s: the output string (must be preallocated)
// n: the string length limit
// d: max depth, -1:unlimited
// out: a count of bytes
//
// this code does NOT assume that the different rounds are the same length
// and it only averages across entries that match
int perftimer_snprintf_csv_body(perftimer_t const * const h, char* s, const size_t n, const unsigned int d) {
  s[0] = '\0'; // string starts empty
  if((h==NULL) || (h->head == NULL)) // nothing to show
    return 0;

  unsigned int m_max = _max_links(h); // longest number of links
  // create list of times
  double* t = calloc(m_max,sizeof(double));
  unsigned int * r = calloc(m_max,sizeof(unsigned int)); // how many went into this count
  perftimer_t const * hs = h;
  perftimer_tic_t const * ptr = h->head;
  while(hs != NULL) {
    ptr = hs->head;
    unsigned int i; // current link
    for(i=0;i<m_max && (ptr != NULL) && (ptr->next != NULL);i++) {
      perftimer_tic_t const * stop = ptr->next;
      while((stop->next != NULL) && (stop->depth > ptr->depth)) {
        stop = stop->next;
      } // summarize lower depths
      t[i] +=  _calc_perftimer_diff(ptr, stop);
      r[i]++;
      ptr = ptr->next;
    }

    hs = hs->old;
  }
  // Note: this could be done without 'r'
  // i.e. t[i] = t[i]*(r-1)/r + new_val/r
  // but this wouldn't be as precise - not sure if it matters
  // but we are dealing with small numbers

  // count desc
  static const int nw = 60;
  char tmp [nw];
  // start at the top, having already established that there is a first entry
  ptr = h->head;
  int i = 0;
  int first = 1;
  while(ptr->next != NULL) {
    if((ptr->desc != NULL) && (ptr->depth <= d)) {
      if(!first) // comma, except first number
        strncat(s, ", ", n-1);
      first = 0;
      snprintf(tmp,nw,"%0.3f", t[i]/((double)r[i])*1e3);
      strncat(s, tmp, n-1);
    }
    ptr = ptr->next;
    i++;
  }
  if(!first) // comma, except first number
    strncat(s, ", ", n-1);
  first = 0;
  snprintf(tmp,nw,"%0.3f",perftimer_wall_av(h)*1e3);
  strncat(s,tmp,n-1);
  unsigned int rnds = perftimer_rounds(h); // rounds
  if(rnds > 1) {
    snprintf(tmp,nw,", %d",rnds);
    strncat(s,tmp,n-1);
  }
  // make the string safe: zero terminate it
  s[n-1] = '\0'; // Note: man claims strncat appends \0 always, but it doesn't

  free(t); // release time list
  free(r);

  // this won't come out to the expected length
  return strnlen(s,n); // total bytes produced
}



// perftimer_printf()
// h: a perftimer structure ptr
// d: max depth, 0:unlimited
// out: printed to stdout (printf)
void perftimer_printf(perftimer_t const * const h, const unsigned int d) {
  size_t n = perftimer_printlen(h,d);
  if(n == 0) // nothing to do
    return;
  char* s = malloc(n*sizeof(char));
  perftimer_snprintf(h,s,n,d);
  printf("%s\n",s);
  free(s);
  return;
}

// TODO refactor common code w/ perftimer_printf
void perftimer_printf_csv_header(perftimer_t const * const h, const unsigned int d) {
  size_t n = perftimer_printlen(h,d); // actual size will be smaller but that's okay
  if(n == 0) // nothing to do
    return;
  char* s = malloc(n*sizeof(char));
  perftimer_snprintf_csv_header(h,s,n,d);
  printf("%s\n",s);
  free(s);
  return;
}

// TODO refactor common code w/ perftimer_printf
void perftimer_printf_csv_body(perftimer_t const * const h, const unsigned int d) {
  size_t n = perftimer_printlen(h,d); // actual size will be smaller but that's okay
  if(n == 0) // nothing to do
    return;
  char* s = malloc(n*sizeof(char));
  perftimer_snprintf_csv_body(h,s,n,d);
  printf("%s\n",s);
  free(s);
  return;
}

// perftimer_wall()
// total time accounted for in the perftimer structure
// h: a perftimer structure ptr
// out: total time
double perftimer_wall(perftimer_t const * const h) {
  perftimer_t const * hh = h;
  while((hh != NULL) && (hh->old != NULL))
    hh = hh->old;
  if(h == NULL) // || hh == null
    return 0.0;
  return _calc_perftimer_diff(hh->head,h->tail);
}
// average over all runs
double perftimer_wall_av(perftimer_t const * h) {
  if(h == NULL)
    return 0.0;
  double sum = 0.0;
  double runs = ((double) perftimer_rounds(h)); // cast to double
  if(runs == 0)
    runs = 1; // avoid div by zero
  while(h != NULL) {
    sum += _calc_perftimer_diff(h->head,h->tail);
    h = h->old;
  }
  return sum / runs;
}

// perftimer_delta()
// the most recent time delta
// h: a perftimer structure ptr
// d: the depth, 0:min
double perftimer_delta(perftimer_t const * const h){
  if((h == NULL) || (h->head == NULL) || (h->head == h->tail))
    return 0.0;

  // We've established that there are at least two entries in
  // the singly linked list: safe to use ptr->next->next
  // Find the second last entry:
  perftimer_tic_t const * ptr = h->head;
  while(ptr->next->next != NULL)
    ptr = ptr->next;

  return _calc_perftimer_diff(ptr, h->tail);
}

// perftimer_rounds()
// h: a perftimer structure ptr
// out: number of timing iterations that have been done
//      effectively perftimer_restart()s +1 unless the next
//      timing round hasn't started yet
unsigned int perftimer_rounds(perftimer_t const* h) {
  unsigned int c = 0;
  while(h != NULL) {
    // if there are at least two tics recorded, count this round toward
    // the total
    if((h->head != NULL) && (h->head != h->tail))
      c++;
    h = h->old;
  }
  return c;
}

// perftimer_restart()
// start timing again from the beginning
// &h: perftimer structure ptr-to-ptr
void perftimer_restart(perftimer_t ** ph) {
  if(ph == NULL) // can't do anything here...
    return;
  perftimer_t* old = *ph;
  if(_perftimer_malloc_head(ph) != 0)
    perftimer_free(old); // release the old head if we failed to allocate
  else
    (*ph)->old = old;
}

// adjust_depth()
// set the current
// h: perftimer structure ptr
// d: the change in depth (int: +- n)
void perftimer_adjust_depth(perftimer_t * const h, int d) {
  if(h == NULL)
    return;
  else if((d < 0) && (-d > h->current_depth))
    h->current_depth = 0;
  else
    h->current_depth += d;
}
