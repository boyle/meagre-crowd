#ifndef _PERFTIMER_H_
#define _PERFTIMER_H_

#include <time.h>
#include <stdio.h>

// doubly linked structure
// TODO switch this to a forward declaration (internal structure is private)
typedef struct perftimer_t {
  struct perftimer_tic_t *head;
  struct perftimer_tic_t *tail;
  unsigned int current_depth;
} perftimer_t;

typedef struct perftimer_tic_t {
  time_t now;
  unsigned int depth;
  struct perftimer_tic_t* next;
  char* desc;
} perftimer_tic_t;

// perftimer_malloc()
// allocate a perftimer structure
// out: a new perftimer ptr, NULL does NOT indicate failure
//  (memory exhaustion will be indicated when you try to perftimer_inc)
perftimer_t* perftimer_malloc();

// perftimer_free()
// release entire perftimer structure T
void perftimer_free(perftimer_t * pT);

// perftimer_inc
// increment the perftimer
// T: a perftimer structure ptr-to-a-ptr
// s: the description of the most recent event completed
// n: max length of string
// d: depth level, 0:min
// out: result - 0: success, -2:failure (memory exhaustion), -1:null ppT
int perftimer_inc(perftimer_t* pT, char const*const s, const size_t n);

// perftimer_snprintf()
// create a string describing the events so far,
// including details upto a depth limit, and restricted to string length n
// T: a perftimer structure ptr
// s: the output string (must be preallocated)
// n: the string length limit
// d: max depth, 0:unlimited
// out: the output string length, <0 means error (TODO what sorts of errors?)
int perftimer_snprintf(perftimer_t const * const pT, char* s, const size_t n, const unsigned int d);

// perftimer_printlen()
// determine the length of the string that would be generated without 
// a constraint on length
// T: a perftimer structure ptr
// d: max depth, 0:unlimited
// out: the required size of the string, not including terminating '\0'
size_t perftimer_printlen(perftimer_t const * const pT, const unsigned int d);

// perftimer_printf()
// T: a perftimer structure ptr
// d: max depth, 0:unlimited
// out: printed to stdout (printf)
void perftimer_printf(perftimer_t const * const pT, const unsigned int d);

// perftimer_wall()
// total time accounted for in the perftimer structure
// T: a perftimer structure ptr
// out: total time
double perftimer_wall(perftimer_t const * const pT);
double perftimer_wall_av(perftimer_t const * const pT);

// perftimer_diff()
// the most recent time delta
// T: a perftimer structure ptr
double perftimer_delta(perftimer_t const * const pT);

// perftimer_rounds()
// T: a perftimer structure ptr
// out: number of timing iterations that have been done
//      effectively perftimer_restart()s +1 unless the next
//      timing round hasn't started yet
unsigned int perftimer_rounds(perftimer_t const* const pT);


// perftimer_restart()
// start timing again from the beginning
// T: perftimer structure ptr-to-ptr
void perftimer_restart(perftimer_t ** ppT);

#endif //_PERFTIMER_H_
