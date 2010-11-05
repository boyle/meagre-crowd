#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "../src/perftimer.h"

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
  perftimer_t* T = perftimer_malloc();

  // check we exit properly when given a NULL ptr-to-ptr
  assert(perftimer_inc(NULL,"",10,0) == -1);
  assert(perftimer_rounds(NULL) == 0);

  // check that the perftimer advances cleanly
  assert(perftimer_inc(&T,"start",10,0) == 0);
// TODO  char* s[100];
  perftimer_printf(T,0);
// TODO  assert(perftimer_printlen(T,0) == strlen()));
  assert(perftimer_wall(T) == 0.0);
  assert(perftimer_diff(T,0) == 0.0);
  
  perftimer_printf(T,0);

  do_nothing(1);

  assert(perftimer_inc(&T,"s1",10,0) == 0);
  assert(perftimer_printlen(T,0) > 5);
  assert(perftimer_wall(T) >= 0.0);
  assert(perftimer_diff(T,0) >= 0.0);

  perftimer_printf(T,0);
  
  do_nothing(200000);
  
  assert(perftimer_inc(&T,"ss1",10,1) == 0);
  assert(perftimer_inc(&T,"ss2",10,1) == 0);
  assert(perftimer_inc(&T,"s3",10,0) == 0);
  assert(perftimer_inc(&T,"sss4",10,2) == 0);
  assert(perftimer_inc(&T,"sss5",10,2) == 0);
  assert(perftimer_inc(&T,"ss6",10,1) == 0);
  for(i=0;i<10;i++)
    perftimer_printf(T,i);

  perftimer_free(T);
}

int main(int argc, char **argv) {
  test_basic();
  return 0;
}
