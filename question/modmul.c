/* Copyright (C) 2017 Daniel Page <csdsp@bristol.ac.uk>
 *
 * Use of this source code is restricted per the CC BY-NC-ND license, a copy of
 * which can be found via http://creativecommons.org (and should be included as
 * LICENSE.txt within the associated archive or repository).
 */

#include "modmul.h"
#include <time.h>

/* Perform stage 1:
 *
 * - read each 3-tuple of N, e and m from stdin,
 * - compute the RSA encryption c, then
 * - write the ciphertext c to stdout.
 */

void windExp(mpz_t t, mpz_t x, mpz_t y2, mpz_t n, int k) {
  mpz_t y,temp;
  mpz_set_ui(t,1);
  mpz_init(y);
  mpz_init(temp);
  mpz_mul(y,y2,t);

  //gmp_printf("23: %Zd\n",y);

  int j = 2;
  j = j << (k-1);

  mpz_t T [j + 2];
  mpz_init_set_ui(T[0],1);

  for (mp_bitcnt_t i = 0; i < (j + 2); i ++) {
    mpz_init(T[i+1]);
    mpz_mul(T[i+1],T[i],x);
    mpz_mod(T[i+1],T[i+1],n);
  }


  mp_bitcnt_t i = (mpz_sizeinbase(y2,2));

  mp_bitcnt_t l;
  mp_bitcnt_t u;
  mp_bitcnt_t temp;

  while (i >= 0) {
    if (mpz_tstbit(y,i) == 0) {
      l = i;
      u = 0;
    }
    else {
      if (i-k+1 > 0) {l = i - k + 1;}
      else {l = 0;}
      u = 0;
      while (mpz_tstbit(y,l) == 0) {l = l + 1;}
      //mpz_mul(temp,y,1);
      //for (int a = 0, a < l, a) {mpz};
      // GRAB THE REQUIRED NUMBER
      //while (temp >= l) {
      //  u = u*2;
      //  if (mpz_tstbit(y,temp) == 1) {u = u + 1;}
      //  temp = temp - 1;
      //}
    }

    for (int j = i - l + 1; j > 0; j --) {mpz_mul(t,t,t);}

    if (u != 0) {mpz_add(t,t,T[(int)((u - 1)/2)]);}
    i = l - 1;
  }
  mpz_clear(y);
}

void stage1(mpz_t c, mpz_t n, mpz_t e, mpz_t m) {
  mpz_powm_sec(c,m,e,n);
  // fill in this function with solution

}

/* Perform stage 2:
 *
 * - read each 9-tuple of N, d, p, q, d_p, d_q, i_p, i_q and c from stdin,
 * - compute the RSA decryption m, then
 * - write the plaintext m to stdout.
 */

void stage2(mpz_t m, mpz_t n, mpz_t d, mpz_t p, mpz_t q, mpz_t d_p, mpz_t d_q, mpz_t i_p, mpz_t i_q, mpz_t c) {
  mpz_powm_sec(m,c,d,n);
  // fill in this function with solution

}

/* Perform stage 3:
 *
 * - read each 5-tuple of p, q, g, h and m from stdin,
 * - compute the ElGamal encryption c = (c_1,c_2), then
 * - write the ciphertext c to stdout.
 */

void stage3(mpz_t c1, mpz_t c2, mpz_t p, mpz_t q, mpz_t g, mpz_t h, mpz_t m) {
  int start,end,mytime;
  start = clock();
  mpz_t r;
  gmp_randstate_t rState;
  mpz_init(r);
  end = clock();
  end = start - end;
  printf("%d\n",end);
  // MUST CHANGE TO RANDOMIZER
  gmp_randinit_default(rState);
  mytime = time(NULL);
  gmp_randseed_ui(rState,mytime);
  mpz_urandomm(r,rState,q);
  gmp_printf("%ZX from %ZX\n",r,q);

  mpz_mod(r,r,q);
  mpz_powm_sec(c1,g,r,p);
  mpz_powm_sec(c2,h,r,p);
  mpz_mul(c2,c2,m);
  mpz_mod(c2,c2,p);
  // fill in this function with solution
  mpz_clear(r);
}

/* Perform stage 4:
 *
 * - read each 5-tuple of p, q, g, x and c = (c_1,c_2) from stdin,
 * - compute the ElGamal decryption m, then
 * - write the plaintext m to stdout.
 */

void stage4(mpz_t m, mpz_t p, mpz_t q, mpz_t g, mpz_t x, mpz_t c1, mpz_t c2) {

  mpz_t e;
  mpz_init(e);
  mpz_mul_si(e,x,-1);
  mpz_mod(e,e,q);
  mpz_powm_sec(m,c1,e,p);
  mpz_mul(m,m,c2);
  mpz_mod(m,m,p);
  // fill in this function with solution

}

/* The main function acts as a driver for the assignment by simply invoking the
 * correct function for the requested stage.
 */

int main( int argc, char* argv[] ) {
  if( 2 != argc ) {
    abort();
  }

  if     ( !strcmp( argv[ 1 ], "stage1" ) ) {
    mpz_t n, e, m, c;

    mpz_init(n);
    mpz_init(e);
    mpz_init(m);
    mpz_init(c);

    if (1 != gmp_scanf("%ZX\n",n)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",e)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",m)) {
      abort();
    }

    stage1(c,n,e,m);

    gmp_printf("%ZX\n",c);

    mpz_clear(n);
    mpz_clear(e);
    mpz_clear(m);
    mpz_clear(c);
  }
  else if( !strcmp( argv[ 1 ], "stage2" ) ) {
    mpz_t n, d, p, q, d_p, d_q, i_p, i_q, c, m;

    mpz_init(n);
    mpz_init(d);
    mpz_init(p);
    mpz_init(q);
    mpz_init(d_p);
    mpz_init(d_q);
    mpz_init(i_p);
    mpz_init(i_q);
    mpz_init(c);
    mpz_init(m);

    if (1 != gmp_scanf("%ZX\n",n)){
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",d)){
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",p)){
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",q)){
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",d_p)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",d_q)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",i_p)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",i_q)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",c)) {
      abort();
    }

    stage2(m,n,d,p,q,d_p,d_q,i_p,i_q,c);

    gmp_printf("%ZX\n",m);

    mpz_clear(n);
    mpz_clear(d);
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(d_p);
    mpz_clear(d_q);
    mpz_clear(i_p);
    mpz_clear(c);
    mpz_clear(m);

  }
  else if( !strcmp( argv[ 1 ], "stage3" ) ) {
// read each 5-tuple of p, q, g, h and m from stdin,
    mpz_t p, q, g, h, m, c1, c2;

    mpz_init(p);
    mpz_init(q);
    mpz_init(g);
    mpz_init(h);
    mpz_init(m);
    mpz_init(c1);
    mpz_init(c2);

    if (1 != gmp_scanf("%ZX\n",p)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",q)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",g)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",h)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",m)) {
      abort();
    }

    stage3(c1,c2,p,q,g,h,m);

    gmp_printf("%ZX\n%ZX\n",c1,c2);

    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(g);
    mpz_clear(h);
    mpz_clear(m);
    mpz_clear(c1);
    mpz_clear(c2);
  }
  else if( !strcmp( argv[ 1 ], "stage4" ) ) {
    mpz_t m, p, q, g, x, c1, c2;
    mpz_init(m);
    mpz_init(p);
    mpz_init(q);
    mpz_init(g);
    mpz_init(x);
    mpz_init(c1);
    mpz_init(c2);

    if (1 != gmp_scanf("%ZX\n",p)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",q)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",g)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",x)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",c1)) {
      abort();
    }
    if (1 != gmp_scanf("%ZX\n",c2)) {
      abort();
    }

    stage4(m,p,q,g,x,c1,c2);

    gmp_printf("%ZX\n",m);

    mpz_clear(m);
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(g);
    mpz_clear(x);
    mpz_clear(c1);
    mpz_clear(c2);
  }
  else if( !strcmp( argv[ 1 ], "test" ) ) {
    mpz_t t,x,y,n;
    mpz_init_set_ui(t,1);
    mpz_init_set_ui(x,1234);
    mpz_init_set_ui(y,3094);
    mpz_init_set_ui(n,503*509);
    int k = 4;

    windExp(t,x,y,n,k);

    mpz_clear(t);
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(n);
  }
  else {
    abort();
  }

  return 0;
}
