#include <stdlib.h> // malloc/free
#include <errno.h>  // malloc error
#include <err.h>    // malloc error text
#include <string.h> // strnlen, strncpy, etc
#include "perftimer.h" 

// helper functions (static = local only)
// go up as far as possible on this branch
static inline perftimer_t const * go_up(perftimer_t const * ptr);
static inline perftimer_t const * go_up(perftimer_t const * ptr) {
  if(ptr == NULL)
    return NULL;

  // go to the head of the list (upwards)
  while(ptr->prev != NULL) {
    ptr = ptr->prev;
  }
  return ptr;
}

// returns NULL if the ptr was NULL, head of the list otherwise
static inline   perftimer_t const *   find_head(perftimer_t const * ptr);
static inline   perftimer_t const *   find_head(perftimer_t const * ptr) {
  return go_up(ptr); // TODO go leftwards if needed (return to main branch)
}


// go down as far as possible on this branch
static inline   perftimer_t const *   go_down(perftimer_t const * ptr);
static inline   perftimer_t const *   go_down(perftimer_t const * ptr) {
  if(ptr == NULL)
    return NULL;

  // go to the head of the list (upwards)
  while(ptr->next != NULL) {
    ptr = ptr->next;
  }
  return ptr;
}

// go up but stop at the first branch
static inline   perftimer_t const *   go_up__nearest_branch(perftimer_t const * ptr);
static inline   perftimer_t const *   go_up__nearest_branch(perftimer_t const * ptr) {
  return go_up(ptr); // TODO look for branches instead
}



// perftimer_malloc()
// allocate a perftimer structure
// out: a new perftimer ptr, NULL does NOT indicate failure
//  (memory exhaustion will be indicated when you try to perftimer_inc)
inline perftimer_t* perftimer_malloc() {
  return NULL; // we start out with an empty list
}

// perftimer_free()
// release entire perftimer structure T
void perftimer_free(perftimer_t * ptr) {
  
  ptr = (perftimer_t*) find_head(ptr); // note: must cast it back to non-const
  if(ptr == NULL)
    return; // nothing to do

  // release memory
  // TODO deal with breadth as well
  perftimer_t* ptr_next;
  while(ptr != NULL) {
    ptr_next = ptr->next;
    free(ptr->desc); // free the copied string
    free(ptr);       // and free itself
    ptr = ptr_next;  // then move forward one
  }
}

// perftimer_inc
// increment the perftimer
// T: a perftimer structure ptr-to-a-ptr
// s: the description of the most recent event completed
// d: depth level, 0:min
// out: result - 0: success, -2: failure (memory exhaustion), -1: NULL ppT
int perftimer_inc(perftimer_t** ppT, char const*const s, const size_t n, const unsigned int d) {
  if(ppT == NULL)
    return -1; // bad ptr

  // find the bottom
  // TODO currently ignore d
  perftimer_t * ptr = (perftimer_t*) find_head(*ppT); // NOTE: this isn't efficient now but needed when looking at 'd' -- TODO
  // cast away const from search
  ptr = (perftimer_t*) go_down(ptr); // cast away const from search

  // allocate a new block
  if((*ppT = malloc(sizeof(perftimer_t))) == NULL) {
    err(1, NULL); // print an error message
    // must leave a ptr to the perftimer structure if there ever was one,
    // on failure, so that the structure doesn't get lost
    *ppT = ptr;
    return -2; // malloc failure
  }

  // fill in the structure
  perftimer_t * const pT = *ppT; // deref ppT since we're accessing it a lot
  pT->now  = time(NULL); // TODO higher resolution
  pT->prev = ptr; // cast away const from the search // TODO refactor const-ness?
  pT->next = NULL;
  { // safe string copy
    // this initial counting is a bit inefficient but worth it for the safety
    size_t nn = strnlen(s,n); 
    if(nn == 0) { // no string
      pT->desc = NULL;
    }
    else { // copy the string in, make sure its null terminated (safety first)
      // TODO use strndup instead?
      pT->desc = malloc((nn+1)*sizeof(char)); // need an extra char for '\0'
      strncpy(pT->desc,s,nn);
      pT->desc[nn] = '\0'; // make sure the string is null terminated
    }
  }
  // and link the old ptr back to the newly allocated object
  // (doubly linked list)
  if(ptr != NULL)
    ptr->next = pT;

  return 0; // success
}

// utility function: calculate difference between times, returns in seconds
// t2 is assumed to be later than t1
static inline   double   calc_perftimer_diff(perftimer_t const* const t1, perftimer_t const* const t2);
static inline   double   calc_perftimer_diff(perftimer_t const* const t1, perftimer_t const* const t2) {
  if((t1 == NULL) || (t2 == NULL)) 
    return 0.0;
  return ((double)(t2->now - t1->now)) / CLOCKS_PER_SEC;
}

// perftimer_wall()
// perftimer_snprintf()
// create a string describing the events so far,
// including details upto a depth limit, and restricted to string length n
// T: a perftimer structure ptr
// s: the output string (must be preallocated)
// n: the string length limit
// d: max depth, 0:unlimited
// out: a count of bytes 
int perftimer_snprintf(perftimer_t const * const pT, char* s, const size_t n, const unsigned int d) {
  s[0] = '\0'; // string starts empty
  perftimer_t const * ptr = find_head(pT); // start at the top

  // count desc
  // TODO handle breadth (d)
  static const int nw = 30;
  char tmp [nw];
  double delta;
  while(ptr != NULL) {
    if(ptr->desc != NULL) {
      strncat(s, ptr->desc, n-1);
      if((delta = calc_perftimer_diff(ptr->prev,ptr)) > 0.0) {
        snprintf(tmp,nw,": %0.3es\n",delta);
      }
      else {
        snprintf(tmp,nw,"\n");
      }
      strncat(s, tmp, n-1);
    }
    ptr = ptr->next;
  }
  snprintf(tmp,nw,"total: %0.3es",perftimer_wall(pT));
  strncat(s,tmp,n-1);
  // make the string safe: zero terminate it
  s[n-1] = '\0'; // Note: man claims strncat appends \0 always, but it doesn't

  // this should come out to the same size as perftimer_printlen()
  // TODO add assertion?
  return strnlen(s,n); // total bytes produced
}

// perftimer_printlen()
// determine the length of the string that would be generated without 
// a constraint on length
// T: a perftimer structure ptr
// d: max depth, 0:unlimited
// out: the required size of the string, excluding terminating '\0'
size_t perftimer_printlen(perftimer_t const * const pT, const unsigned int d) {
  size_t c = 0;
  perftimer_t const * ptr = find_head(pT);

  // count desc
  // TODO handle breadth (d)
  while(ptr != NULL) {
    if(ptr->desc != NULL)
      // we know that all ptr->desc are zero terminated: we copied them
      // need to leave lots of room for times
      c += strlen(ptr->desc) +30;
    // TODO test for overflow on 'c += ...'
    ptr = ptr->next;
  }

  return c+30-1;
}

// perftimer_printf()
// T: a perftimer structure ptr
// d: max depth, 0:unlimited
// out: printed to stdout (printf)
void perftimer_printf(perftimer_t const * const pT, const unsigned int d) {
  size_t n = perftimer_printlen(pT,d);
  char* s = malloc((n+1+500)*sizeof(char)); // TODO rm 500
  perftimer_snprintf(pT,s,n,d);
  printf("%s\n\n",s);
  free(s);
  return;
}

// total time accounted for in the perftimer structure
// T: a perftimer structure ptr
// out: total time
double perftimer_wall(perftimer_t const * const pT) {
  perftimer_t const * const head = find_head(pT);
  perftimer_t const * const tail = go_down(head);

  // TODO look for branches and go down those too
  
  // TODO use higher resolution perftimers? time.h -> sys/time.h
  return calc_perftimer_diff(head,tail);
}
double perftimer_wall_av(perftimer_t const * const pT) {
  return perftimer_wall(pT); // TODO
}

// perftimer_diff()
// the most recent time delta at this depth
// T: a perftimer structure ptr
// d: the depth, 0:min
double perftimer_diff(perftimer_t const * const pT, const unsigned int d){
  perftimer_t const * const last = go_down(find_head(pT));
  if((last == NULL) || (last->prev == NULL))
    return 0.0;

  return calc_perftimer_diff(last->prev, last);
}
double perftimer_diff_av(perftimer_t const * const pT, const unsigned int d) {
  return perftimer_diff(pT,d);
}

// perftimer_rounds()
// T: a perftimer structure ptr
// out: number of timing iterations that have been done
//      effectively perftimer_restart()s +1 unless the next
//      timing round hasn't started yet
unsigned int perftimer_rounds(perftimer_t const* const pT) {
  if(pT == NULL)
    return 0;
  else
    return 1; // TODO or more...
}

// perftimer_restart()
// start timing again from the beginning
// T: perftimer structure ptr-to-ptr
void perftimer_restart(perftimer_t ** ppT) {
  perftimer_free(*ppT); // TODO make a new list instead of destroying the old one
  *ppT = perftimer_malloc();
}

