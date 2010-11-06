#include <stdlib.h> // malloc/free
#include <errno.h>  // malloc error
#include <err.h>    // malloc error text
#include <string.h> // strnlen, strncpy, etc
#include "perftimer.h" 



// perftimer_malloc()
// allocate a perftimer structure
// out: a new perftimer ptr, NULL does NOT indicate failure
//  (memory exhaustion will be indicated when you try to perftimer_inc)
inline perftimer_t* perftimer_malloc() {
  perftimer_t* ptr;
  if((ptr = malloc(sizeof(perftimer_t))) == NULL) {
    warn("perftimer_malloc()"); // shows errno string
    return NULL;
  }
  ptr->head = NULL;
  ptr->tail = NULL;
  ptr->current_depth = 0;
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
  free(h);
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
  ptr->now  = time(NULL); // TODO higher resolution
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
static inline   double   calc_perftimer_diff(perftimer_tic_t const* const t1, perftimer_tic_t const* const t2);
static inline   double   calc_perftimer_diff(perftimer_tic_t const* const t1, perftimer_tic_t const* const t2) {
  if((t1 == NULL) || (t2 == NULL)) 
    return 0.0;
  return ((double)(t2->now - t1->now)) / CLOCKS_PER_SEC;
}

// perftimer_printlen()
// determine the length of the string that would be generated without 
// a constraint on length
// h: a perftimer structure ptr
// d: max depth, 0:unlimited
// out: the required size of the string, excluding terminating '\0'
size_t perftimer_printlen(perftimer_t const * const h, const unsigned int d) {
  if((h == NULL) || (h->head == NULL))
    return 0;

  size_t c = 0;
  perftimer_tic_t const * ptr = h->head;

  // count desc
  // TODO handle breadth (d)
  while(ptr != NULL) {
    if(ptr->desc != NULL)
      // we know that all ptr->desc are zero terminated: we copied them
      // need to leave lots of room for calculated time
      c += strlen(ptr->desc) +30;
    // TODO test for overflow on 'c += ...'
    ptr = ptr->next;
  }

  return c+30-1; // TODO extra padding?
}

// perftimer_snprintf()
// create a string describing the events so far,
// including details upto a depth limit, and restricted to string length n
// h: a perftimer structure ptr
// s: the output string (must be preallocated)
// n: the string length limit
// d: max depth, 0:unlimited
// out: a count of bytes 
int perftimer_snprintf(perftimer_t const * const h, char* s, const size_t n, const unsigned int d) {
  s[0] = '\0'; // string starts empty
  if((h==NULL) || (h->head == NULL)) // nothing to show
    return 0;

  // start at the top, having already established that there is a first entry
  perftimer_tic_t const * ptr = h->head;

  // count desc
  // TODO handle breadth (d)
  static const int nw = 30;
  char tmp [nw];
  double delta;
  while(ptr->next != NULL) {
    if(ptr->desc != NULL) {
      strncat(s, ptr->desc, n-1);
      if((delta = calc_perftimer_diff(ptr,ptr->next)) > 0.0) {
        snprintf(tmp,nw,": %0.3es\n",delta);
      }
      else {
        snprintf(tmp,nw,"\n");
      }
      strncat(s, tmp, n-1);
    }
    ptr = ptr->next;
  }
  snprintf(tmp,nw,"total: %0.3es",perftimer_wall(h));
  strncat(s,tmp,n-1);
  // make the string safe: zero terminate it
  s[n-1] = '\0'; // Note: man claims strncat appends \0 always, but it doesn't

  // this should come out to the same size as perftimer_printlen()
  // TODO add assertion?
  return strnlen(s,n); // total bytes produced
}


// perftimer_printf()
// h: a perftimer structure ptr
// d: max depth, 0:unlimited
// out: printed to stdout (printf)
void perftimer_printf(perftimer_t const * const h, const unsigned int d) {
  size_t n = perftimer_printlen(h,d);
  char* s = malloc(n*sizeof(char));
  perftimer_snprintf(h,s,n,d);
  printf("%s\n\n",s);
  free(s);
  return;
}

// perftimer_wall()
// total time accounted for in the perftimer structure
// h: a perftimer structure ptr
// out: total time
inline double perftimer_wall(perftimer_t const * const h) {
  return calc_perftimer_diff(h->head,h->tail);
}
double perftimer_wall_av(perftimer_t const * const h) {
  return perftimer_wall(h); // TODO
}

// perftimer_diff()
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

  return calc_perftimer_diff(ptr, h->tail);
}

// perftimer_rounds()
// h: a perftimer structure ptr
// out: number of timing iterations that have been done
//      effectively perftimer_restart()s +1 unless the next
//      timing round hasn't started yet
unsigned int perftimer_rounds(perftimer_t const* const h) {
  if(h == NULL)
    return 0;
  else
    return 1; // TODO or more...
}

// perftimer_restart()
// start timing again from the beginning
// &h: perftimer structure ptr-to-ptr
void perftimer_restart(perftimer_t ** ph) {
  perftimer_free(*ph); // TODO make a new list instead of destroying the old one
  *ph = perftimer_malloc();
}

