#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <signal.h>
#include <sys/sysinfo.h>
#include <pthread.h>

#define MAXBits 300

void handl(int signum)
{
    exit(0);
}

typedef struct
{
    mpz_t x;
    mpz_t y;
} point;
void clearPoint(point *p)
{
    mpz_clear(p->x);
    mpz_clear(p->y);
}
typedef struct
{
    mpz_t a;
    mpz_t b;
    mpz_t n;
    mpz_t x0;
    mpz_t y0;
} courbeElliptique;
void clearCourbe(courbeElliptique *c)
{
    mpz_clear(c->a);
    mpz_clear(c->b);
    mpz_clear(c->n);
    mpz_clear(c->x0);
    mpz_clear(c->y0);
}

void trouveBorne(mpz_t newborne, mpz_t borneinit)
{
    mpz_t a;
    mpz_init(a);
    mpz_set_si(a, 2);
    mpz_t prod;
    mpz_init(prod);
    mpz_set_si(prod, 1);
    while (mpz_cmp(prod, borneinit) < 0)
    {
        mpz_mul(prod, prod, a);
        mpz_nextprime(a, a);
    }
    // if (mpz_cmp(prod, borneinit)>0){
    //    mpz_sub_ui(a,a,1);
    //   mpz_div(prod, prod, a);
    //}
    mpz_set(newborne, a);
    mpz_clear(a);
    mpz_clear(prod);
}
void addition(point *res2, point p1, point p2, courbeElliptique c, mpz_t *valide)
{

    if (mpz_cmp_si(p1.x, 0) == 0 && mpz_cmp_si(p1.y, 0) == 0)
    {

        mpz_set(res2->x, p2.x);
        mpz_set(res2->y, p2.y);
        return;
    }
    if (mpz_cmp_si(p2.x, 0) == 0 && mpz_cmp_si(p2.y, 0) == 0)
    {
        mpz_set(res2->x, p1.x);
        mpz_set(res2->y, p1.y);
        return;
    }

    mpz_t pgcd;
    mpz_init(pgcd);
    mpz_t sub;
    mpz_init(sub);
    mpz_sub(sub, p1.x, p2.x);
    mpz_mod(sub, sub, c.n);

    mpz_gcd(pgcd, sub, c.n);
    // gmp_printf("pgcd 1 : %Zd\n", pgcd);
    // gmp_printf("pgcd : %Zd\n", pgcd);
    if (mpz_cmp_d(pgcd, 1) == 0)
    {
        mpz_t lambda;
        mpz_init(lambda);
        mpz_t x3;
        mpz_t y3;
        mpz_init(x3);
        mpz_init(y3);

        mpz_invert(lambda, sub, c.n);
        mpz_sub(sub, p1.y, p2.y);
        mpz_mod(sub, sub, c.n);
        mpz_mul(lambda, sub, lambda);
        mpz_mod(lambda, lambda, c.n);
        mpz_mul(sub, lambda, lambda);
        mpz_mod(sub, sub, c.n);
        mpz_sub(sub, sub, p1.x);
        mpz_mod(sub, sub, c.n);
        mpz_sub(x3, sub, p2.x);
        mpz_mod(x3, x3, c.n);
        mpz_sub(sub, p1.x, x3);
        mpz_mod(sub, sub, c.n);
        mpz_mul(sub, sub, lambda);
        mpz_mod(sub, sub, c.n);
        mpz_sub(y3, sub, p1.y);
        mpz_mod(y3, y3, c.n);
        mpz_set(res2->x, x3);
        mpz_set(res2->y, y3);
        mpz_clear(sub);
        mpz_clear(pgcd);

        mpz_clear(lambda);
        mpz_clear(x3);
        mpz_clear(y3);

        return;
    }
    else if (mpz_cmp(pgcd, c.n) == 0)
    {

        mpz_add(sub, p1.y, p2.y);

        mpz_mod(sub, sub, c.n);
        mpz_gcd(pgcd, sub, c.n);

        if (mpz_cmp_d(pgcd, 1) != 0 && mpz_cmp(pgcd, c.n) != 0)
        {
            printf("FACTEUR TROUVE\n");
            mpz_set(*valide, pgcd);
            mpz_clear(pgcd);
            // gmp_printf("valide = %Zd", *valide);
            mpz_clear(sub);
            mpz_set(res2->x, p1.x);
            mpz_set(res2->y, p1.y);
            return;
        }
        if (mpz_cmp(pgcd, c.n) == 0)
        {
            // printf("pgcd(y) egal a n\n ")

            mpz_clear(pgcd);
            mpz_clear(sub);
            mpz_set_si(res2->x, 0);
            mpz_set_si(res2->y, 0);
            return;
        }
        if (mpz_cmp_d(pgcd, 1) == 0)
        {

            // printf("pgcd(y) egal a 1\n");
            mpz_t lambda;
            mpz_t lambdacarre;
            mpz_init(lambdacarre);
            mpz_t x3;
            mpz_t y3;
            mpz_init(y3);
            mpz_init(x3);
            mpz_init(lambda);
            mpz_mul(lambda, p1.x, p1.x);   // x^2
            mpz_mod(lambda, lambda, c.n);  // mod n
            mpz_mul_ui(lambda, lambda, 3); // 3x^2

            mpz_mod(lambda, lambda, c.n);

            mpz_add(lambda, lambda, c.a);

            mpz_invert(sub, sub, c.n);

            mpz_mul(lambda, sub, lambda);

            mpz_mod(lambda, lambda, c.n);
            mpz_mul(lambdacarre, lambda, lambda);

            mpz_sub(sub, lambdacarre, p1.x);
            mpz_sub(sub, sub, p2.x);
            mpz_mod(sub, sub, c.n);
            mpz_set(x3, sub);

            mpz_sub(sub, p1.x, x3);

            mpz_mul(sub, sub, lambda);

            mpz_mod(sub, sub, c.n);

            mpz_sub(sub, sub, p1.y);

            mpz_mod(sub, sub, c.n);
            mpz_set(res2->x, x3);
            mpz_set(res2->y, sub);
            mpz_clear(lambda);
            mpz_clear(sub);
            mpz_clear(x3);
            mpz_clear(y3);
            mpz_clear(pgcd);
            mpz_clear(lambdacarre);
            return;
        }
    }
    else
    {
        mpz_set(*valide, pgcd);
        printf("FACTEUR 2 ");
        mpz_clear(pgcd);
        mpz_clear(sub);
        mpz_set(res2->x, p1.x);
        mpz_set(res2->y, p1.y);
        return;
    }
}
void buildCurve(courbeElliptique *c, mpz_t n, gmp_randstate_t state)
{

    mpz_t a;
    mpz_init(a);
    mpz_urandomm(a, state, n);
    mpz_urandomm(a, state, n);
    mpz_mod(a, a, n);
    mpz_init(c->a);
    mpz_set(c->a, a);
    mpz_t x0;
    mpz_init(x0);
    mpz_urandomm(x0, state, n);
    mpz_mod(x0, x0, n);
    mpz_t y0;
    mpz_init(y0);
    mpz_urandomm(y0, state, n);
    mpz_mod(y0, y0, n);
    mpz_t b;
    mpz_t xcube;
    mpz_t ax0;
    mpz_init(ax0);
    mpz_init(xcube);
    mpz_init(b);
    mpz_mul(b, y0, y0);
    mpz_mul(xcube, x0, x0);
    mpz_mul(xcube, xcube, x0);
    mpz_mul(ax0, a, x0);
    mpz_sub(b, b, xcube);
    mpz_sub(b, b, ax0);
    mpz_mod(b, b, n);

    mpz_init(c->b);
    mpz_set(c->b, b);
    mpz_init(c->n);
    mpz_set(c->n, n);
    mpz_init(c->x0);
    mpz_init(c->y0);
    mpz_set(c->x0, x0);
    mpz_set(c->y0, y0);
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(x0);
    mpz_clear(y0);
    mpz_clear(xcube);
    mpz_clear(ax0);
}
int racinecaree(int nb)
{
    int i = 0;
    while (i * i < nb)
    {
        i++;
    }
    return i;
}
int expbinaire(int nb, int puissance)
{
    int res = 1;
    while (puissance > 0)
    {
        if (puissance % 2 == 1)
        {
            res = res * nb;
        }
        nb = nb * nb;
        puissance = puissance / 2;
    }
    return res;
}
void bonnepuissance(mpz_t *res, mpz_t nb, mpz_t borne)
{
    mpz_t anciennb;
    mpz_init(anciennb);
    mpz_set(anciennb, nb);
    // gmp_printf("appel de bonnepuissance sur nb = %Zd, borne = %Zd\n", nb, borne);
    mpz_t nbmul;
    mpz_init(nbmul);
    mpz_set(nbmul, nb);
    mpz_t nb2;
    mpz_init(nb2);
    mpz_sqrt(nb2, borne);
    mpz_mul_ui(nb2, nb2, 2);
    mpz_add_ui(nb2, nb2, 1);
    mpz_add(nb2, borne, nb2);
    while (mpz_cmp(nb, nb2) < 0)
    {
        mpz_mul(nb, nbmul, nb);
        //  gmp_printf("nombre = %Zd\n, borne = %Zd\n", nb, nb2);
    }
    if (mpz_cmp(nb, borne) > 0)
    {
        mpz_div(nb, nb, nbmul);
    }
    mpz_set(*res, nb);
    mpz_set(nb, anciennb);
    mpz_clear(anciennb);
    mpz_clear(nb2);
    mpz_clear(nbmul);
    // gmp_printf("res = %Zd\n", *res);
}
int degreeBin(mpz_t nb)
{
    return mpz_sizeinbase(nb, 2);
}
void multiplication(point *res1, point p1, mpz_t k, courbeElliptique c, mpz_t *valide)
{ // Multiplication de p1 par k, meme principe que l'exponentiation binaire

    point res;
    mpz_init(res.x);
    mpz_init(res.y);
    mpz_set_d(res.x, 0);
    mpz_set_d(res.y, 0);
    point ad;
    mpz_init(ad.x);
    mpz_init(ad.y);
    mpz_set(ad.x, p1.x);
    mpz_set(ad.y, p1.y);
    for (size_t i = 0; i < mpz_sizeinbase(k, 2); i++)
    {

        if (mpz_tstbit(k, i) == 1)
        {
            if (i == 0)
            {

                addition(&res, res, ad, c, valide);
                // gmp_printf("point : x = %Zd\n y = %Zd\n", p1.x, p1.y);
                if (mpz_cmp_si(*valide, 1) != 0)
                {
                    return;
                }
            }
            else
            {

                addition(&res, res, ad, c, valide);

                if (mpz_cmp_si(*valide, 1) != 0)
                {
                    return;
                }
                addition(&ad, ad, ad, c, valide);

                if (mpz_cmp_si(*valide, 1) != 0)
                {
                    return;
                }
            }
        }
        else
        {

            addition(&ad, ad, ad, c, valide);

            if (mpz_cmp_si(*valide, 1) != 0)
            {
                return;
            }
        }
    }

    mpz_set(res1->x, res.x);
    mpz_set(res1->y, res.y);
    clearPoint(&res);
    clearPoint(&ad);

    /*point res;
    mpz_init(res.x);
    mpz_init(res.y);
    mpz_set_d(res.x, 0);
    mpz_set_d(res.y, 0);
    mpz_t i;
    mpz_init(i);
    mpz_set(i, k);
    gmp_printf("%Zd\n", k);
    while(mpz_cmp_ui(i,0)!=0){
        //gmp_printf(" i = %Zd\n",i);
        addition(&res, res, p1, c, valide);
        //gmp_printf("p1.y = %Zd\n",p1.y );
        //gmp_printf("res.x = %Zd\n", res.x);
        //gmp_printf("res.y = %Zd\n", res.y);
        if(mpz_cmp_ui(*valide, 1)!=0){
                    return;
                }
        mpz_sub_ui(i,i,1);
    }
    */
}
void ECM(mpz_t n, gmp_randstate_t state, int B2, mpz_t borne)
{
    time_t start = time(NULL);
    int j = 1;
    courbeElliptique c;

    while (1)
    {
        // printf("j = %d\n",j);
        mpz_t verif1;
        mpz_t verif2;
        mpz_t pgcd;

        mpz_init(verif2);

        do
        {
            buildCurve(&c, n, state);

            mpz_init(pgcd);
            mpz_init(verif1);
            mpz_set_si(verif1, 4);
            mpz_set_si(verif2, 27);
            mpz_mul(verif1, verif1, c.a);
            mpz_mul(verif1, verif1, c.a);
            mpz_mul(verif1, verif1, c.a);
            mpz_mul(verif2, verif2, c.b);
            mpz_mul(verif2, verif2, c.b);
            mpz_add(verif1, verif1, verif2);
            mpz_gcd(pgcd, verif1, n);
            mpz_clear(verif1);
            mpz_clear(verif2);
        } while (mpz_cmp(pgcd, n) == 0);
        //printf("Courbe : \n");
        //gmp_printf("a : %Zd\n", c.a);
        //gmp_printf("b : %Zd\n", c.b);
        //gmp_printf("n : %Zd\n", c.n);

        if (mpz_cmp_d(pgcd, 1) != 0 && mpz_cmp(pgcd, n) != 0)
        {
            gmp_printf("facteur trouvé à la %d eme itération : %Zd\n", j, pgcd);
            clearCourbe(&c);
            mpz_clear(pgcd);
            _exit(0);
        }

        else
        { // On choisit une borne de friabilité
            mpz_clear(pgcd);

            mpz_t pow;
            mpz_init(pow);
            point p1;
            mpz_init(p1.x);
            mpz_init(p1.y);
            mpz_t valide;
            mpz_init(valide);
            mpz_set_si(valide, 1);
            mpz_set(p1.x, c.x0);
            mpz_set(p1.y, c.y0);
            // gmp_printf("point : x = %Zd\n y = %Zd\n", p1.x, p1.y);
            mpz_t i;
            mpz_init(i);
            mpz_set_si(i, 2);

            while (mpz_cmp(borne, i) > 0)
            {
                bonnepuissance(&pow, i, borne);
                // gmp_printf("multiplication de p1: x =%Zd, y = %Zd par %Zd\n",p1.x, p1.y, pow);
                multiplication(&p1, p1, pow, c, &valide);
                //  gmp_printf("point : x = %Zd\n y = %Zd\n", p1.x, p1.y);
                if (mpz_cmp_si(valide, 1) != 0)
                {
                    time_t end = time(NULL);
                    gmp_printf("facteur trouvé en %d secondes à l'itération %d : %Zd\n", end - start, j, valide);
                    mpz_clear(valide);
                    clearPoint(&p1);
                    clearCourbe(&c);
                    _exit(0);
                }
                // gmp_printf("%Zd\n", i);
                mpz_nextprime(i, i);

                // mpz_mul(prod, i, prod);
            }
            if (B2 != 0)
            {
                //printf("entrée dans la phase 2\n");
                // B2 PHASE :
                mpz_t newborne;
                mpz_t diff;
                mpz_init(newborne);
                mpz_set(newborne, borne);
                mpz_mul_si(newborne, newborne, 100);
                mpz_init(diff);
                mpz_set(diff, i);
                point ad;
                point p2;
                mpz_init(p2.x);
                mpz_init(p2.y);
                mpz_set(p2.x, p1.x);
                mpz_set(p2.y, p1.y);

                mpz_init(ad.x);
                mpz_init(ad.y);
                point tab[100];
                for (int k = 0; k < 100; k++)
                {
                    mpz_t k2;
                    mpz_init(k2);
                    mpz_set_si(k2, 2 + 2 * k);
                    mpz_init(tab[k].x);
                    mpz_init(tab[k].y);
                    point temp;
                    mpz_init(temp.x);
                    mpz_init(temp.y);
                    multiplication(&temp, p1, k2, c, &valide);
                    mpz_set(tab[k].x, temp.x);
                    mpz_set(tab[k].y, temp.y);

                    mpz_clear(k2);
                    clearPoint(&temp);
                }
                multiplication(&p2, p1, diff, c, &valide);
                if (mpz_cmp_si(valide, 1) != 0)
                {
                    time_t end = time(NULL);
                    gmp_printf("facteur trouvé en %d secondes à l'itération %d : %Zd\n", end - start, j, valide);
                    mpz_clear(valide);
                    clearPoint(&p1);
                    clearCourbe(&c);
                    _exit(0);
                }
                else
                {
                    while (mpz_cmp(newborne, i) > 0)
                    {
                        mpz_nextprime(i, i);
                        mpz_sub(diff, i, diff);

                        mpz_div_ui(diff, diff, 2);
                        mpz_sub_ui(diff, diff, 1);
                        int indice = mpz_get_ui(diff);
                        // gmp_printf("i = %Zd, newborne = %Zd, diff = %Zd, x : %Zd, y:%Zd\n",i, newborne, diff, p1.x, p1.y);

                        if (indice < 100)
                        {
                            addition(&p2, p2, tab[indice], c, &valide);
                            if (mpz_cmp_si(valide, 1) != 0)
                            {
                                time_t end = time(NULL);
                                gmp_printf("B2 : facteur trouvé en %d secondes à l'itération %d : %Zd\n", end - start, j, valide);
                                mpz_clear(valide);
                                clearPoint(&p1);
                                clearCourbe(&c);
                                _exit(0);
                            }
                        }
                        mpz_set(diff, i);
                    }
                }
            }
            mpz_clear(pow);
            // mpz_clear(newborne);
            mpz_clear(i);

            // gmp_printf("prod = :%Zd ", prod);
            // mpz_clear(prod);
            mpz_clear(valide);
            clearPoint(&p1);
        }
        // printf("Aucun facteur trouvé pour la %d eme courbe\n", i);

        j += 1;
        time_t tour = time(NULL);
        if (tour - start > 180)
        {
            clearCourbe(&c);
            printf("Aucun facteur trouvé\n");
            return;
        }
    }
    clearCourbe(&c);
    printf("Aucun facteur trouvé\n");
}

typedef struct thread_data
{
    mpz_t *n;
    gmp_randstate_t *state;
    mpz_t *borne;
    mpz_t *p;
    mpz_t *q;
} thread_data_t;

void *find_factor(void *threadarg)
{
    thread_data_t *my_data = (thread_data_t *)threadarg;
    ECM(*my_data->n, *my_data->state, 1, *my_data->borne);
    gmp_randclear(*my_data->state);
    mpz_clear(*my_data->n);
    mpz_clear(*my_data->p);
    mpz_clear(*my_data->q);
    mpz_clear(*my_data->borne);
    return NULL;
}

int main()
{
    mpz_t n;
    mpz_init(n);
    // mpz_set_str(n, "216703300309849365159018393570821",10);
    mpz_t p;
    mpz_init(p);
    mpz_t q;
    mpz_init(q);
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, time(NULL));
    mpz_urandomb(q, state, 70);
    mpz_nextprime(q, q);
    mpz_urandomb(p, state, 70);
    mpz_nextprime(p, p);

    mpz_mul(n, q, p);
    gmp_printf(" n = %Zd = %Zd*%Zd\n", n, p, q);

    mpz_t borne;
    mpz_init(borne);
    mpz_set_str(borne, "50000", 10);

    pid_t fork1 = fork();
    if (fork1 == 0)
    {
        wait(NULL);
        gmp_randclear(state);
        mpz_clear(n);
        mpz_clear(p);
        mpz_clear(q);
        mpz_clear(borne);
        return 0;
    }
    else
    {
        pid_t fork2 = fork();
        if (fork2 > 0)
        {
            const int nb_threads = get_nprocs();
            pthread_t threads[nb_threads];

            thread_data_t data;
            data.n = &n;
            data.state = &state;
            data.borne = &borne;
            data.p = &p;
            data.q = &q;

            for (int i = 0; i < nb_threads; i++)
            {
                pthread_create(&threads[i], NULL, find_factor, &data);
            }

            for (int i = 0; i < nb_threads; i++)
            {
                pthread_join(threads[i], NULL);
            }
        } else if (fork2 == 0) {
            exit(0);
        }
    }

    return 0;
}
