/* Meagre-Crowd: A sparse distributed matrix solver testbench for performance benchmarking.
 * Copyright (C) 2010 Alistair Boyle <alistair.js.boyle@gmail.com>
 *
 *     This file is part of Meagre-Crowd.
 *
 *     Meagre-Crowd program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _PERFTIMER_H_
#define _PERFTIMER_H_

#include "config.h"
#include <sys/time.h>
#include <stdio.h>

// doubly linked structure
// TODO switch this to a forward declaration (internal structure is private)
typedef struct perftimer_t {
  struct perftimer_tic_t *head;
  struct perftimer_tic_t *tail;
  struct perftimer_t* old;
  unsigned int current_depth;
} perftimer_t;

typedef struct perftimer_tic_t {
  struct timeval now;
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
void perftimer_free( perftimer_t * pT );

// perftimer_inc
// increment the perftimer
// T: a perftimer structure ptr-to-a-ptr
// s: the description of the most recent event completed
// n: max length of string
// d: depth level, 0:min
// out: result - 0: success, -2:failure (memory exhaustion), -1:null ppT
int perftimer_inc( perftimer_t* pT, char const*const s, const size_t n );

// perftimer_snprintf()
// create a string describing the events so far,
// including details upto a depth limit, and restricted to string length n
// T: a perftimer structure ptr
// s: the output string (must be preallocated)
// n: the string length limit
// d: max depth, 0:unlimited
// out: the output string length
int perftimer_snprintf( perftimer_t const * const pT, char* s, const size_t n, const unsigned int d );
int perftimer_snprintf_csv_header( perftimer_t const * const pT, char* s, const size_t n, const unsigned int d );
int perftimer_snprintf_csv_body( perftimer_t const * const pT, char* s, const size_t n, const unsigned int d );

// perftimer_printlen()
// determine the length of the string that would be generated without
// a constraint on length
// T: a perftimer structure ptr
// d: max depth, 0:unlimited
// out: the required size of the string, not including terminating '\0'
size_t perftimer_printlen( perftimer_t const * const pT, const unsigned int d );

// perftimer_printf()
// T: a perftimer structure ptr
// d: max depth, 0:unlimited
// out: printed to stdout (printf)
void perftimer_printf( perftimer_t const * const pT, const unsigned int d );
void perftimer_printf_csv_header( perftimer_t const * const pT, const unsigned int d );
void perftimer_printf_csv_body( perftimer_t const * const pT, const unsigned int d );

// perftimer_wall()
// total time accounted for in the perftimer structure
// T: a perftimer structure ptr
// out: total time
double perftimer_wall( perftimer_t const * const pT );
double perftimer_wall_av( perftimer_t const * pT );

// perftimer_diff()
// the most recent time delta
// T: a perftimer structure ptr
double perftimer_delta( perftimer_t const * const pT );

// perftimer_rounds()
// T: a perftimer structure ptr
// out: number of timing iterations that have been done
//      effectively perftimer_restart()s +1 unless the next
//      timing round hasn't started yet
unsigned int perftimer_rounds( perftimer_t const* pT );


// perftimer_restart()
// start timing again from the beginning
// T: perftimer structure ptr-to-ptr
void perftimer_restart( perftimer_t ** ppT );

// adjust_depth()
// set the current
// h: perftimer structure ptr
// d: the change in depth (int: +- n)
void perftimer_adjust_depth( perftimer_t * const h, int d );


#endif //_PERFTIMER_H_
