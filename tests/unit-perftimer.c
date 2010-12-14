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
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "perftimer.h"

// do something to kill time
void do_nothing(int i);
void do_nothing(int i) {
  int a[10] = {0};
  int j;
  for(j=0;j<i*1000;j++)
    a[j%10] = i+j*(i-2);
}

void test_basic();
void test_basic() {
  unsigned int i;
  
  // check we exit properly when given a NULL ptr
  assert(perftimer_inc(NULL,"",10) == -1);
  perftimer_free(NULL);
  char s [10];
  assert(perftimer_snprintf(NULL,s,10,-1) == 0);
  assert(perftimer_snprintf_csv_header(NULL,s,10,-1) == 0);
  assert(perftimer_snprintf_csv_body(NULL,s,10,-1) == 0);
  assert(perftimer_printlen(NULL,-1) == 0);
  perftimer_printf(NULL,-1);
  perftimer_printf_csv_header(NULL,-1);
  perftimer_printf_csv_body(NULL,-1);
  assert(perftimer_wall(NULL) == 0.0);
  assert(perftimer_wall_av(NULL) == 0.0);
  assert(perftimer_delta(NULL) == 0.0);
  assert(perftimer_rounds(NULL) == 0);
  perftimer_restart(NULL);
  perftimer_t* t = NULL;
  perftimer_restart(&t);
  perftimer_free(t); // we now have a valid timer after restart, need to free it
  perftimer_adjust_depth(NULL,+10);
  perftimer_adjust_depth(NULL,-10);
  perftimer_adjust_depth(NULL,0);

  perftimer_t* T = perftimer_malloc();
  assert(perftimer_rounds(T) == 0);

  // check that the perftimer advances cleanly
  assert(perftimer_inc(T,"start",10) == 0);
  assert(perftimer_rounds(T) == 0);
  perftimer_printf(T,0);
  printf("csv:\n");
  perftimer_printf_csv_header(T,0);
  perftimer_printf_csv_body(T,0);
  assert(perftimer_wall(T) == 0.0);
  assert(perftimer_delta(T) == 0.0);
  
  perftimer_printf(T,0);
  perftimer_printf_csv_header(T,0);
  perftimer_printf_csv_body(T,0);

  do_nothing(10000);

  assert(perftimer_inc(T,"s1",10) == 0);
  assert(perftimer_rounds(T) == 1);
  assert(perftimer_printlen(T,0) > 5);
  assert(perftimer_wall(T) >= 0.0);
  assert(perftimer_delta(T) >= 0.0);

  perftimer_printf(T,0);
  perftimer_printf_csv_header(T,0);
  perftimer_printf_csv_body(T,0);
  
  do_nothing(200000);
  
  perftimer_adjust_depth(T,+1);
  assert(perftimer_inc(T,"ss1",10) == 0);
  perftimer_adjust_depth(T,+1);
  assert(perftimer_inc(T,"ss2",10) == 0);
  perftimer_adjust_depth(T,-2);
  assert(perftimer_inc(T,"s3",10) == 0);
  assert(perftimer_inc(T,"sss4",10) == 0);
  perftimer_adjust_depth(T,+6);
  assert(perftimer_inc(T,"sss5",10) == 0);
  assert(perftimer_inc(T,"ss6",10) == 0);
  for(i=0;i<10;i++) {
    printf("depth %d\n",i);
    perftimer_printf(T,i);
    perftimer_printf_csv_header(T,i);
    perftimer_printf_csv_body(T,i);
  }

  // try a restart
  assert(perftimer_rounds(T) == 1);
  perftimer_restart(&T);
  assert(perftimer_rounds(T) == 1);
  assert(perftimer_inc(T,"ss1",10) == 0);
  assert(perftimer_rounds(T) == 1);
  assert(perftimer_inc(T,"ss2",10) == 0);
  assert(perftimer_rounds(T) == 2);
  assert(perftimer_inc(T,"ss3",10) == 0);
  perftimer_printf(T,0);
  perftimer_printf_csv_header(T,0);
  perftimer_printf_csv_body(T,0);

  // try wall_av(): average wall time
  assert(perftimer_rounds(T) == 2);
  printf("wall_av() = %fs in %d reps\n",perftimer_wall_av(T),perftimer_rounds(T));
  assert(perftimer_wall_av(T) != 0.0);
  assert(perftimer_wall_av(T) < perftimer_wall(T));
  assert(perftimer_wall_av(T) < 1.0e6); // not infinity

  perftimer_free(T);
}

int main(int argc, char **argv) {
  test_basic();
  return 0;
}
