/* Copyright (C) 2017 Daniel Page <csdsp@bristol.ac.uk>
 *
 * Use of this source code is restricted per the CC BY-NC-ND license, a copy of
 * which can be found via http://creativecommons.org (and should be included as
 * LICENSE.txt within the associated archive or repository).
 */

#include "modmul.h"
#include <time.h>
#include <stdio.h>
#include <stdbool.h>

const int windExpK = 8;
const int montK = 5;



/* Grabs relevant r_i from a base 2^k:
 *
 * - Takes as input r_i, r, k and i
 * - sets r_i <- r_i from r in base 2^k
 */

void grabBits(mpz_t r_i, mpz_t r, int k, int i) {
  int r_j = 0;
  int l = i*k;
  for (int j = (k - 1); j >= 0; j --) {
    r_j = r_j << 1;
    r_j = r_j + mpz_tstbit(r,(j+l));
  }
  mpz_set_ui(r_i,r_j);
}

/* Conducts Montgomery multiplication:
 *
 * - Takes as input r, x, y, p, w, n and k
 * - sets r <- x*y mod n
 */

void montMul(mpz_t r, mpz_t x, mpz_t y, mpz_t p, mpz_t w, mpz_t n, int k) {
  mpz_set_ui(r,0);
  int b = 1 << k;
  int j = (strlen(mpz_get_str(NULL,b,n)));

  mpz_t u, r_0, y_i, x_0, t0, t1, bZ;
  mpz_init_set_ui(u,0);
  mpz_init_set_ui(r_0,0);
  mpz_init_set_ui(y_i,0);
  mpz_init_set_ui(x_0,0);
  mpz_init_set_ui(t1,0);
  mpz_init_set_ui(t0,0);
  mpz_init_set_ui(bZ,b);
  grabBits(x_0,x,k,0);

  for (int i = 0; i < j; i++) {
    grabBits(r_0,r,k,0);
    grabBits(x_0,x,k,0);
    grabBits(y_i,y,k,i);

    mpz_mul(u,y_i,x_0);
    mpz_add(u,r_0,u);
    mpz_mul(u,u,w);
    mpz_mod(u,u,bZ);

    mpz_mul(t0,y_i,x);
    mpz_mul(t1,u,n);
    mpz_add(t1,t0,t1);
    mpz_add(t1,r,t1);
    mpz_fdiv_q_2exp(r,t1,k);
  }

  if (mpz_cmp(r,n) >= 0) {
    mpz_sub(r,r,n);
  }
  mpz_clear(u);
  mpz_clear(r_0);
  mpz_clear(y_i);
  mpz_clear(x_0);
  mpz_clear(t0);
  mpz_clear(t1);
  mpz_clear(bZ);
}

/* Sets up variables for Montgomery multiplication:
 *
 * - Takes as input w, n and rho
 * - sets omega <- -n^-1 mod rho
 */

void montOmega(mpz_t w, mpz_t n, mpz_t p) {
  mpz_sub(w,p,n);
  mpz_invert(w,w,p);
}

/* Sets up xHat for Montgomery multiplication:
 *
 * - Takes as input xHat, p, w, x, n and k
 * - sets xHat <- x*p mod n
 */

void montHat(mpz_t xHat, mpz_t p, mpz_t w, mpz_t x, mpz_t n, int k) {
  mpz_t p2;
  mpz_init(p2);
  mpz_mul(p2,p,p);
  mpz_mod(p2,p2,n);
  montMul(xHat,x,p2,p,w,n,k);
  mpz_clear(p2);
}

/* Gets back x value from xHat that was used for Montgomery multiplication:
 *
 * - Takes as input x, xHat, p, w, n and k
 * - sets x <- xHat in normal representation
 */

void deMont(mpz_t x, mpz_t xHat, mpz_t p, mpz_t w, mpz_t n, int k) {
  mpz_t one;
  mpz_init_set_ui(one,1);
  montMul(x,xHat,one,p,w,n,k);
  mpz_clear(one);
}

/* Perform windowed exponentiation without Montgomery:
 *
 * - Takes as input t, x, y, n and k
 * - sets t <- x^y mod n
 */

void windExpNoMont(mpz_t t, mpz_t x, mpz_t y, mpz_t n, int k) {
  mpz_set_ui(t,1);
  //mpz_init(y);
  //mpz_mul(y,y2,t);

  int j = 2;
  j = j << (k-2);

  mpz_t T [j + 2];
  mpz_init_set(T[0],x);
  for (mp_bitcnt_t i = 0; i < (j + 2); i ++) {
    mpz_init(T[i+1]);
    mpz_mul(T[i+1],T[i],x);
    mpz_mul(T[i+1],T[i+1],x);
    mpz_mod(T[i+1],T[i+1],n);
  }

  mp_bitcnt_t i = (strlen(mpz_get_str(NULL,2,y)));

  mp_bitcnt_t l;
  mp_bitcnt_t u;

  while ((i + 1) > 0)  {

    if (mpz_tstbit(y,i) == 0) {
      l = i;
      u = 0;
    }
    else {
      j = i-k+1;
      if (j > 0) {l = i - k + 1;}
      else {l = 0;}
      u = 0;
      while (mpz_tstbit(y,l) == 0) {
        l = l + 1;
      }
      //mpz_mul(temp,y,1);
      //for (int a = 0, a < l, a) {mpz};
      for (j = i; (j + 1) > l; j --) {
        u = u << 1 ;
        u = u + mpz_tstbit(y,j);
      }
    }

    for (int j = i - l + 1; j > 0; j --) {
      mpz_mul(t,t,t);
      mpz_mod(t,t,n);
    }

    if (u != 0) {
      mpz_mul(t,t,T[(int)((u - 1)/2)]);
      mpz_mod(t,t,n);
    }
    i = l - 1;
  }
  for (mp_bitcnt_t i = 0; i < (j + 3); i ++) {
    mpz_clear(T[i]);
  }
}

/* Sets up rho for Montgomery multiplication:
 *
 * - Takes as input p, n and k
 * - sets p <- required rho value
 */

void montRho(mpz_t p, mpz_t n, int k) {

  mpz_t b,two,kz;
  mpz_init(b);
  mpz_init_set_ui(two,2);
  mpz_init_set_ui(kz,k);
  windExpNoMont(b,two,kz,n,windExpK);
  mpz_set(p,b);
  while (mpz_cmp(p,n) < 0) {
    mpz_mul(p,p,b);
  }
}

/* Sets up variables for Montgomery multiplication:
 *
 * - Takes as input xHat, p, w, x, n and k
 * - sets p <- required rho value
 * - sets omega <- -n^-1 mod rho
 * - sets xHat <- x*p mod n
 */

void montSetup(mpz_t xHat, mpz_t p, mpz_t w, mpz_t x, mpz_t n, int k) {
  montRho(p,n,k);
  montOmega(w,n,p);
  montHat(xHat,p,w,x,n,k);
}

/* Perform windowed exponentiation with external Montgomery:
 *
 * - Takes as input tHat = t*p (mod n), xHat = x*p (mod n), y, p , w, n and k
 * - sets tHat <- xHat^y mod n
 */

void windExp(mpz_t tH, mpz_t xH, mpz_t y, mpz_t p, mpz_t w, mpz_t n, int k) {
  mpz_t one, zH;
  mpz_init_set_ui(one,1);
  mpz_init(zH);
  montHat(tH,p,w,one,n,k);
  int j = 2;
  j = j << (k-2);
  mpz_t T [j + 2];
  mpz_init_set(T[0],xH);
  for (mp_bitcnt_t i = 0; i < (j + 2); i ++) {
    mpz_init(T[i+1]);
    montMul(T[i+1],T[i],xH,p,w,n,montK);
    montMul(zH,T[i+1],xH,p,w,n,montK);
    mpz_set(T[i+1],zH);
  }
  mp_bitcnt_t i = (strlen(mpz_get_str(NULL,2,y)));

  mp_bitcnt_t l;
  mp_bitcnt_t u;

  while ((i + 1) > 0)  {
    if (mpz_tstbit(y,i) == 0) {
      l = i;
      u = 0;
    }
    else {
      j = i-k+1;
      if (j > 0) {l = i - k + 1;}
      else {l = 0;}
      u = 0;
      while (mpz_tstbit(y,l) == 0) {
        l = l + 1;
      }
      //mpz_mul(temp,y,1);
      //for (int a = 0, a < l, a) {mpz};
      for (j = i; (j + 1) > l; j --) {
        u = u << 1 ;
        u = u + mpz_tstbit(y,j);
      }
    }
    for (int j = i - l + 1; j > 0; j --) {
      montMul(zH,tH,tH,p,w,n,montK);
      mpz_set(tH,zH);
    }
    if (u != 0) {
      montMul(zH,tH,T[(int)((u - 1)/2)],p,w,n,montK);
      mpz_set(tH,zH);
    }
    i = l - 1;
  }
  for (mp_bitcnt_t i = 0; i < (j + 3); i ++) {
    mpz_clear(T[i]);
  }
  mpz_clear(zH);
  mpz_clear(one);
}

/* Perform windowed exponentiation with internal Montgomery:
 *
 * - Takes as input t, x, y, n and k
 * - sets t <- x^y mod n
 */

void windExpInternalMont(mpz_t t, mpz_t x, mpz_t y, mpz_t n, int k) {
  mpz_set_ui(t,1);
  mpz_t w, p, xH, tH, zH;
  mpz_init(w);
  mpz_init(p);
  mpz_init(xH);
  mpz_init(tH);
  mpz_init(zH);
  int montK = 2;

  montSetup(xH,p,w,x,n,montK);
  montHat(tH,p,w,t,n,montK);

  int j = 2;
  j = j << (k-2);

  mpz_t T [j + 2];
  mpz_init_set(T[0],xH);
  for (mp_bitcnt_t i = 0; i < (j + 2); i ++) {
    mpz_init(T[i+1]);
    montMul(T[i+1],T[i],xH,p,w,n,montK);
    montMul(zH,T[i+1],xH,p,w,n,montK);
    mpz_set(T[i+1],zH);
  }
  mp_bitcnt_t i = (strlen(mpz_get_str(NULL,2,y)));
  mp_bitcnt_t l;
  mp_bitcnt_t u;

  while ((i + 1) > 0)  {
    if (mpz_tstbit(y,i) == 0) {
      l = i;
      u = 0;
    }
    else {
      j = i-k+1;
      if (j > 0) {l = i - k + 1;}
      else {l = 0;}
      u = 0;
      while (mpz_tstbit(y,l) == 0) {
        l = l + 1;
      }
      //mpz_mul(temp,y,1);
      //for (int a = 0, a < l, a) {mpz};
      for (j = i; (j + 1) > l; j --) {
        u = u << 1 ;
        u = u + mpz_tstbit(y,j);
      }
    }
    for (int j = i - l + 1; j > 0; j --) {
      montMul(zH,tH,tH,p,w,n,montK);
      mpz_set(tH,zH);
    }
    if (u != 0) {
      montMul(zH,tH,T[(int)((u - 1)/2)],p,w,n,montK);
      mpz_set(tH,zH);
    }
    i = l - 1;
  }
  for (mp_bitcnt_t i = 0; i < (j + 3); i ++) {
    mpz_clear(T[i]);
  }
  deMont(t,tH,p,w,n,montK);
  mpz_clear(w);
  mpz_clear(p);
  mpz_clear(tH);
  mpz_clear(xH);
  mpz_clear(zH);
}

//////////////////////////////////////

/* Perform stage 1:
 *
 * - read each 3-tuple of N, e and m from stdin,
 * - compute the RSA encryption c, then
 * - write the ciphertext c to stdout.
 */

void stage1(mpz_t c, mpz_t n, mpz_t e, mpz_t m) {
  //mpz_t cH,p,w,mH;
  //mpz_init(cH);
  //mpz_init(p);
  //mpz_init(w);
  //mpz_init(mH);
  //montSetup(mH,p,w,m,n,montK);
  //windExp(cH,mH,e,p,w,n,windExpK);
  //deMont(c,cH,p,w,n,montK);

  if (1 != gmp_scanf("%ZX\n",n)) {
    abort();
  }

  if (1 != gmp_scanf("%ZX\n",e)) {
    abort();
  }

  if (1 != gmp_scanf("%ZX\n",m)) {
    abort();
  }

  windExpNoMont(c,m,e,n,windExpK);
  // fill in this function with solution
  gmp_printf("%ZX\n",c);
  //mpz_clear(cH);
  //mpz_clear(p);
  //mpz_clear(w);
  //mpz_clear(mH);
}

/* Perform stage 2:
 *
 * - read each 9-tuple of N, d, p, q, d_p, d_q, i_p, i_q and c from stdin,
 * - compute the RSA decryption m, then
 * - write the plaintext m to stdout.
 */

void stage2(mpz_t m, mpz_t n, mpz_t d, mpz_t p, mpz_t q, mpz_t d_p, mpz_t d_q, mpz_t i_p, mpz_t i_q, mpz_t c) {
  mpz_t m1,m2/*,m1H,m2H,mH,omega,rho,cH*/;
  mpz_init(m1);
  mpz_init(m2);
  //mpz_init(m1H);
  //mpz_init(m2H);
  //mpz_init(omega);
  //mpz_init(rho);
  //mpz_init(cH);
  //mpz_init(mH);

/*  montSetup(cH,rho,omega,c,p,montK);

  windExp(m1H,cH,d_p,rho,omega,p,windExpK);

  deMont(m1,m1H,rho,omega,p,montK);

  montSetup(cH,rho,omega,c,q,montK);

  windExp(m2H,cH,d_q,rho,omega,q,windExpK);

  deMont(m2,m2H,rho,omega,q,montK);

  montSetup(m2H,rho,omega,m2,n,montK);
  montHat(m1H,rho,omega,m1,n,montK);
*/
  //i_p = p^-1 (mod q)
  //i_q = q^-1 (mod p)
  //deMont(m1,m1H,rho,omega,n,montK);
  //deMont(m2,m2H,rho,omega,n,montK);
  //if (mpz_cmp(p,q) < 0) {
/*    mpz_sub(mH,m2H,m1H);
    montHat(m2H,rho,omega,i_p,n,montK);
    montMul(mH,mH,m2H,rho,omega,n,montK);
    montHat(m2H,rho,omega,p,n,montK);
    montMul(mH,mH,m2H,rho,omega,n,montK);
    mpz_add(mH,mH,m1H);
    deMont(m,mH,rho,omega,n,montK); */
  //}
  //else {
  //  mpz_sub(m,m1,m2);
  //  mpz_mul(m,m,i_q);
  //  mpz_mod(m,m,n);
  //  mpz_mul(m,m,q);
  //  mpz_add(m,m,m2);
  //  mpz_mod(m,m,n);
  //}

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

  windExpNoMont(m1,c,d_p,p,windExpK);
  windExpNoMont(m2,c,d_q,q,windExpK);

  //i_p = p^-1 (mod q)
  //i_q = q^-1 (mod p)

  //if (mpz_cmp(p,q) < 0) {
    mpz_sub(m,m2,m1);
    mpz_mul(m,m,i_p);
    mpz_mod(m,m,n);
    mpz_mul(m,m,p);
    mpz_add(m,m,m1);
    mpz_mod(m,m,n);

  gmp_printf("%ZX\n",m);

  //windExp(m,c,d,n,windExpK);

  mpz_clear(m1);
  mpz_clear(m2);
  //mpz_clear(m1H);
  //mpz_clear(m2H);
  //mpz_clear(omega);
  //mpz_clear(rho);
  //mpz_clear(cH);
  //mpz_clear(mH);
}

/* Perform stage 3:
 *
 * - read each 5-tuple of p, q, g, h and m from stdin,
 * - compute the ElGamal encryption c = (c_1,c_2), then
 * - write the ciphertext c to stdout.
 */

void stage3(mpz_t c1, mpz_t c2, mpz_t p, mpz_t q, mpz_t g, mpz_t h, mpz_t m) {
  mpz_t r/*,c1H,c2H,ghat,hhat,mH, rho,omega*/;
  gmp_randstate_t rState;
  mpz_init(r);
  //mpz_init(c2H);
  //mpz_init(c1H);
  //mpz_init(rho);
  //mpz_init(omega);
  //mpz_init(ghat);
  //mpz_init(hhat);
  //mpz_init(mH);

  //montSetup(ghat,rho,omega,g,p,montK);
  //montHat(hhat,rho,omega,h,p,montK);
  //montHat(mH,rho,omega,m,p,montK);

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

  gmp_randinit_default(rState);

  int randomSeed = 0;

  setbuf(stdout, NULL);
  FILE *randomData = fopen("/dev/urandom", "r");
  char randInp[50];
  fgets(randInp, sizeof randInp, randomData);

  fclose(randomData);
  for (int ij = 0; ij < 50; ij++){
    randomSeed = (randomSeed*10) + (int)(randInp[ij] - '0');
  }

  gmp_randseed_ui(rState,randomSeed);
  mpz_sub_ui(q,q,1);
  mpz_urandomm(r,rState,q);

  /*
  mpz_mod(r,r,q);
  windExp(c1H,ghat,r,rho,omega,p,windExpK);
  windExp(c2H,hhat,r,rho,omega,p,windExpK);
  montMul(hhat,c2H,mH,rho,omega,p,montK);
  deMont(c2,hhat,rho,omega,p,montK);
  deMont(c1,c1H,rho,omega,p,montK);
  // fill in this function with solution
  mpz_clear(r);
  mpz_clear(c2H);
  mpz_clear(c1H);
  mpz_clear(rho);
  mpz_clear(omega);
  mpz_clear(ghat);
  mpz_clear(hhat);
  mpz_clear(mH); */
  mpz_add_ui(r,r,1);
  windExpNoMont(c1,g,r,p,windExpK);
  windExpNoMont(c2,h,r,p,windExpK);
  mpz_mul(c2,c2,m);
  mpz_mod(c2,c2,p);
  // fill in this function with solution
  gmp_printf("%ZX\n%ZX\n",c1,c2);
  mpz_clear(r);

}

/* Perform stage 4:
 *
 * - read each 5-tuple of p, q, g, x and c = (c_1,c_2) from stdin,
 * - compute the ElGamal decryption m, then
 * - write the plaintext m to stdout.
 */

void stage4(mpz_t m, mpz_t p, mpz_t q, mpz_t g, mpz_t x, mpz_t c1, mpz_t c2) {

  mpz_t e/*, c1H, c2H, mH, zH, rho, omega*/;
  mpz_init(e);
  //mpz_init(c1H);
  //mpz_init(c2H);
  //mpz_init(mH);
  //mpz_init(rho);
  //mpz_init(omega);
  //mpz_init(zH);

  //montSetup(c1H,rho,omega,c1,p,montK);
  //montHat(c2H,rho,omega,c2,p,montK);
  //montHat(mH,rho,omega,m,p,montK);

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

  mpz_neg(e,x);
  mpz_mod(e,e,q);

/*  windExp(mH,c1H,e,rho,omega,p,windExpK);
  montMul(zH,mH,c2,rho,omega,p,montK);
  deMont(m,zH,rho,omega,p,montK); */

  windExpNoMont(m,c1,e,p,windExpK);
  mpz_mul(m,m,c2);
  mpz_mod(m,m,p);
  // fill in this function with solution
  gmp_printf("%ZX\n",m);

  mpz_clear(e);
  //mpz_clear(c1H);
  //mpz_clear(c2H);
  //mpz_clear(mH);
  //mpz_clear(rho);
  //mpz_clear(omega);
  //mpz_clear(zH);
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

    stage1(c,n,e,m);

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


    stage2(m,n,d,p,q,d_p,d_q,i_p,i_q,c);

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

    stage3(c1,c2,p,q,g,h,m);

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

    stage4(m,p,q,g,x,c1,c2);

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
    mpz_init(t);
    mpz_init(x);
    mpz_init(y);
    mpz_init(n);
    mpz_set_ui(t,1);
    mpz_set_ui(x,16224235);
    mpz_set_ui(y,1121214);
    int intN = 2190323525;
    mpz_set_ui(n,intN);


    if (1 != gmp_scanf("%ZX\n",x)) {
      abort();
    }

    if (1 != gmp_scanf("%ZX\n",y)) {
      abort();
    }

    if (1 != gmp_scanf("%ZX\n",n)) {
      abort();
    }

    //if (1 != gmp_scanf("%ZX\n",x)) {
    //  abort();
    //}
    //if (1 != gmp_scanf("%ZX\n",y)) {
    //  abort();
    //}
    //if (1 != gmp_scanf("%ZX\n",n)) {
    //  abort();
    //}

    gmp_printf("Finished Reading in\n");

    int k = 8;


    mpz_t w,p,xH,yH,zH,z;
    mpz_init_set_ui(p,1);
    mpz_init_set_ui(w,1);
    mpz_init_set_ui(xH,1);
    mpz_init_set_ui(yH,1);
    mpz_init_set_ui(zH,1);
    mpz_init_set_ui(z,1);
    k = montK;

    montSetup(xH,p,w,x,n,k);
    montHat(yH,p,w,y,n,k);

    windExp(zH,xH,y,p,w,n,k);
    deMont(t,zH,p,w,n,k);
    mpz_mod(t,t,n);
    gmp_printf("Full Mont Wind Exp gives: %Zd\n", t);
    windExpInternalMont(t,x,y,n,k);
    mpz_mod(t,t,n);
    gmp_printf("Half Mont Wind Exp gives: %Zd\n", t);
    windExpNoMont(t,x,y,n,k);
    gmp_printf("Result should be: %Zd\n", t);
    mpz_clear(t);

    montMul(zH,xH,yH,p,w,n,k);
    deMont(z,zH,p,w,n,k);
    gmp_printf("Mont gives: %Zd\n", z);

    mpz_mul(z,x,y);
    mpz_mod(z,z,n);
    gmp_printf("Result should be: %Zd\n", z);

    mpz_sub(zH,xH,yH);
    deMont(z,zH,p,w,n,k);
    gmp_printf("Add in mont = %Zd\n", z);

    mpz_sub(z,x,y);
    mpz_mod(z,z,n);
    gmp_printf("Add in mont = %Zd\n", z);

    //for (int ij = 0; ij < intN; ij++) {
    //  montHat(xH,p,w,x,n,k);
    //  deMont(y,xH,p,w,n,k);
    //  if (mpz_cmp(x,y) != 0) {
    //    gmp_printf(".");
    //  }
    //}


    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(n);
  }
  else {
    abort();
  }

  return 0;
}
