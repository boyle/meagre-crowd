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
  assert(perftimer_printlen(NULL,-1) == 0);
  perftimer_printf(NULL,-1);
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
  assert(perftimer_wall(T) == 0.0);
  assert(perftimer_delta(T) == 0.0);
  
  perftimer_printf(T,0);

  do_nothing(10000);

  assert(perftimer_inc(T,"s1",10) == 0);
  assert(perftimer_rounds(T) == 1);
  assert(perftimer_printlen(T,0) > 5);
  assert(perftimer_wall(T) >= 0.0);
  assert(perftimer_delta(T) >= 0.0);

  perftimer_printf(T,0);
  
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
