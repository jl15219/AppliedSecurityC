/* Copyright (C) 2017 Daniel Page <csdsp@bristol.ac.uk>
 *
 * Use of this source code is restricted per the CC BY-NC-ND license, a copy of
 * which can be found via http://creativecommons.org (and should be included as
 * LICENSE.txt within the associated archive or repository).
 */

#include "modmul.h"
#include <time.h>

// NEEDS TO BE BIGGER THAN 4 FOR NO REASON
const int windExpK = 8;

/* Perform stage 1:
 *
 * - read each 3-tuple of N, e and m from stdin,
 * - compute the RSA encryption c, then
 * - write the ciphertext c to stdout.
 */

int rdrand64( int* r ) {
 int success;
 asm( "rdrand %0 ; setc %1"
      : "=r" (*r), "=qm" (success) );
 return success;
}


void montRed(mpz_t r, mpz_t t, mpz_t w, mpz_t p, mpz_t n, int k) {
  mpz_set(r,t);
  int b = 1 << k;
  int j = (strlen(mpz_get_str(NULL,b,n)));
  mpz_t u;
  mpz_init(u);
  for (int i = 0; i < j; i++) {

  }
}

void grabBits(mpz_t r_i, mpz_t r, int k, int i) {
  int r_j = 0;
  int l = i*k;
  for (int j = (k - 1); j >= 0; j --) {
    r_j = r_j << 1;
    r_j = r_j + mpz_tstbit(r,(j+l));
  }
  mpz_set_ui(r_i,r_j);
}


void montMul(mpz_t r, mpz_t x, mpz_t y, mpz_t p, mpz_t w, mpz_t n, int k) {
  mpz_set_ui(r,0);
  int b = 1 << k;
  int j = (strlen(mpz_get_str(NULL,b,n)));
  //gmp_printf("j = %d\n",j);

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
    //gmp_printf("%d)a)t1 = %Zd <- (u = %Zd)*(n = %Zd)\n",i,t1, u,n);
    mpz_add(u,r_0,u);
    mpz_mul(u,u,w);
    mpz_mod(u,u,bZ);

    mpz_mul(t0,y_i,x);
    mpz_mul(t1,u,n);
    //gmp_printf("%d)a)t1 = %Zd <- (u = %Zd)*(n = %Zd)\n",i,t1, u,n);
    mpz_add(t1,t0,t1);
    //gmp_printf("%d)a)t1 = %Zd <- (t0 = %Zd)+t1\n",i,t1, t0);
    mpz_add(t1,r,t1);
    //gmp_printf("%d)a)t1 = %Zd <- (r = %Zd)*(t1)\n",i,t1, r);
    mpz_fdiv_q_2exp(r,t1,k);
    //gmp_printf("r = %Zd <- (t1 = %Zd)/(2^%d)", r ,t1, k);
    //gmp_printf("%d) r = %Zd, r_0 = %Zd, u = %Zd, x_0 = %Zd, y_i = %Zd, x = %Zo, y = %Zo, t0 = %Zd, t1 = %Zd, n = %Zd\n", i, r, r_0, u, x_0, y_i, x, y, t0, t1,n);
  }
  //gmp_printf("\n");
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

void montOmega(mpz_t w, mpz_t n, mpz_t p) {
  mpz_sub(w,p,n);
  mpz_invert(w,w,p);
}

void montHat(mpz_t xHat, mpz_t p, mpz_t w, mpz_t x, mpz_t n, int k) {
  mpz_t p2;
  mpz_init(p2);
  mpz_mul(p2,p,p);
  mpz_mod(p2,p2,n);
  montMul(xHat,x,p2,p,w,n,k);
  mpz_clear(p2);
}

void deMont(mpz_t x, mpz_t xHat, mpz_t p, mpz_t w, mpz_t n, int k) {
  mpz_t one;
  mpz_init_set_ui(one,1);
  montMul(x,xHat,one,p,w,n,k);
  mpz_clear(one);
}

void windExpTest(mpz_t t, mpz_t x, mpz_t y, mpz_t n, int k) {
  mpz_set_ui(t,1);
  //mpz_init(y);
  //mpz_mul(y,y2,t);
  //gmp_printf("23: %Zd\n",y);

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

  mp_bitcnt_t i = (mpz_sizeinbase(y,2));

  mp_bitcnt_t l;
  mp_bitcnt_t u;

  while ((i + 1) > 0)  {
    //gmp_printf("%Zd\n", t);
    if (mpz_tstbit(y,i) == 0) {
      l = i;
      u = 0;
    }
    else {
      j = i-k+1;
      if (j > 0) {l = i - k + 1;}
      else {l = 0;}
      //gmp_printf("%d",l);
      u = 0;
      while (mpz_tstbit(y,l) == 0) {
        l = l + 1;
        //gmp_printf(", %d", l);
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

void montRho(mpz_t p, mpz_t n, int k) {

  mpz_t b,two,kz;
  mpz_init(b);
  mpz_init_set_ui(two,2);
  mpz_init_set_ui(kz,k);
  windExpTest(b,two,kz,n,windExpK);
  mpz_set(p,b);
  while (mpz_cmp(p,n) < 0) {
    mpz_mul(p,p,b);
  }
}

void montSetup(mpz_t xHat, mpz_t p, mpz_t w, mpz_t x, mpz_t n, int k) {
  //gmp_printf("Creating Rho\n");
  montRho(p,n,k);
  //gmp_printf("Rho made\nCreating Omega\n");
  montOmega(w,n,p);
  //gmp_printf("Omega made\nCreating xHat\n");
  montHat(xHat,p,w,x,n,k);
  //gmp_printf("xHat made\n");
}

void windExp(mpz_t t, mpz_t x, mpz_t y, mpz_t n, int k) {
  mpz_set_ui(t,1);
  //mpz_init(y);
  //mpz_mul(y,y2,t);
  //gmp_printf("23: %Zd\n",y);
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
  mp_bitcnt_t i = (mpz_sizeinbase(y,2));

  mp_bitcnt_t l;
  mp_bitcnt_t u;

  while ((i + 1) > 0)  {
    //gmp_printf("%Zd\n", t);
    if (mpz_tstbit(y,i) == 0) {
      l = i;
      u = 0;
    }
    else {
      j = i-k+1;
      if (j > 0) {l = i - k + 1;}
      else {l = 0;}
      //gmp_printf("%d",l);
      u = 0;
      while (mpz_tstbit(y,l) == 0) {
        l = l + 1;
        //gmp_printf(", %d", l);
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

void stage1(mpz_t c, mpz_t n, mpz_t e, mpz_t m) {
  windExp(c,m,e,n,windExpK);
  //mpz_powm_sec(c,m,e,n);
  // fill in this function with solution
  gmp_printf("%ZX\n",c);
}

/* Perform stage 2:
 *
 * - read each 9-tuple of N, d, p, q, d_p, d_q, i_p, i_q and c from stdin,
 * - compute the RSA decryption m, then
 * - write the plaintext m to stdout.
 */

void stage2(mpz_t m, mpz_t n, mpz_t d, mpz_t p, mpz_t q, mpz_t d_p, mpz_t d_q, mpz_t i_p, mpz_t i_q, mpz_t c) {
  mpz_t m1,m2;
  mpz_init(m1);
  mpz_init(m2);

  windExp(m1,c,d_p,p,windExpK);
  windExp(m2,c,d_q,q,windExpK);

  //i_p = p^-1 (mod q)
  //i_q = q^-1 (mod p)

  //if (mpz_cmp(p,q) < 0) {
    mpz_sub(m,m2,m1);
    mpz_mul(m,m,i_p);
    mpz_mod(m,m,n);
    mpz_mul(m,m,p);
    mpz_add(m,m,m1);
    mpz_mod(m,m,n);
  //}
  //else {
  //  mpz_sub(m,m1,m2);
  //  mpz_mul(m,m,i_q);
  //  mpz_mod(m,m,n);
  //  mpz_mul(m,m,q);
  //  mpz_add(m,m,m2);
  //  mpz_mod(m,m,n);
  //}


  gmp_printf("%ZX\n",m);

  //windExp(m,c,d,n,windExpK);
  // fill in this function with solution
  //gmp_printf("%ZX\n",m);
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

  int randomData;
  int * randomPointer;
  randomPointer = &randomData;
  rdrand64(randomPointer);

  //int randomData = open("/dev/urandom", O_RDONLY);
  //if (randomData < 0) {
  //  myclock = time(NULL);
  //}
  //else
  //{
  //  char myRandomData[50];
  //  ssize_t result = read(randomData, myRandomData, sizeof myRandomData);
  //  if (result < 0)
  //  {
  //      // something went wrong
  //  }
  //}
  gmp_randseed_ui(rState,4);
  mpz_urandomm(r,rState,q);
  gmp_printf("%ZX from %ZX\n",r,q);

  mpz_mod(r,r,q);
  windExp(c1,g,r,p,windExpK);
  windExp(c2,h,r,p,windExpK);
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

  mpz_t e;
  mpz_init(e);
  mpz_mul_si(e,x,-1);
  mpz_mod(e,e,q);
  windExp(m,c1,e,p,windExpK);
  mpz_mul(m,m,c2);
  mpz_mod(m,m,p);
  // fill in this function with solution
  gmp_printf("%ZX\n",m);
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
    mpz_set_ui(x,16224);
    mpz_set_ui(y,1124);
    int intN = 21903;
    mpz_set_ui(n,intN);

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
    windExp(t,x,y,n,k);
    gmp_printf("Mont Wind Exp gives: %Zd\n", t);
    windExpTest(t,x,y,n,k);
    gmp_printf("Result should be: %Zd\n", t);
    mpz_clear(t);

    mpz_t w,p,xH,yH,zH,z;
    mpz_init_set_ui(p,1);
    mpz_init_set_ui(w,1);
    mpz_init_set_ui(xH,1);
    mpz_init_set_ui(yH,1);
    mpz_init_set_ui(zH,1);
    mpz_init_set_ui(z,1);
    k = 3;

    montSetup(xH,p,w,x,n,k);
    montHat(yH,p,w,y,n,k);
    //gmp_printf("\nN = %Zd, b = %d, Rho = %Zd, Omega = %Zd, xHat = %Zd, yHat = %Zd\n",n,1 << k,p,w,xH,yH);


    montMul(zH,xH,yH,p,w,n,k);
    deMont(z,zH,p,w,n,k);
    gmp_printf("Mont gives: %Zd\n", z);

    mpz_mul(z,x,y);
    mpz_mod(z,z,n);
    gmp_printf("Result should be: %Zd\n", z);

    for (int ij = 0; ij < intN; ij++) {
      montHat(xH,p,w,x,n,k);
      deMont(y,xH,p,w,n,k);
      if (mpz_cmp(x,y) != 0) {
        gmp_printf("NOOOO\n");
      }
    }


    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(n);
  }
  else {
    abort();
  }

  return 0;
}
