/* Untimed independent all-entry verification for the final lattice kernel.
 * This file never calls the benchmark generator, a floating reference dot,
 * the candidate's recurrence, or any timer. Its primary data are int32_t.
 * Four fixed held-out cases reproduce lattice_benchmark.c's xorshift64 stream
 * and construction order using integer coefficients; conversion to double is
 * exclusively scalbn((double)X,e). The oracle is a direct signed int64 dot.
 * For these cases n<=8192 and |X|<=32768, every partial absolute sum is <=2^43,
 * below the independently checked 2^46 limit and int64 overflow by a wide gap.
 * Build/run is delegated to algebra_lattice_full_entry_verify.py --run.
 */
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>
#include <fenv.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int lattice_sequential_gram(const double *,size_t,size_t,size_t,
                            double *,size_t,int,size_t);

struct test_case {
    const char *family;
    size_t m,n,q;
    int e,dist;
    uint64_t seed;
};
static const struct test_case cases[] = {
    {"windows",512,4096,2,-15,0,UINT64_C(20261201)},
    {"innov",768,8192,4,0,0,UINT64_C(20261201)},
    {"innov",512,4096,2,-15,1,UINT64_C(20261201)},
    {"innov",384,2048,2,-15,2,UINT64_C(20261201)}
};

/* Independent PRNG state; arithmetic is defined modulo 2^64. */
static uint64_t random_state;
static uint64_t next_random(void) {
    random_state ^= random_state << 13;
    random_state ^= random_state >> 7;
    random_state ^= random_state << 17;
    return random_state;
}
static int32_t integer_coefficient(int distribution,size_t index) {
    switch (distribution) {
        case 0: return (int32_t)(next_random() & UINT64_C(65535))-32768;
        case 1: return 16384+(int32_t)(next_random()%UINT64_C(17))-8;
        case 2: return (index&1u)?-32768:32767;
        case 3: return index%17u?0:((index&1u)?-32768:32767);
        default:return (int32_t)(next_random()%UINT64_C(3))-1;
    }
}
static int make_integer_input(const struct test_case *c,int32_t *x) {
    const size_t m=c->m,n=c->n;
    random_state=c->seed?c->seed:UINT64_C(1);
    if (!strcmp(c->family,"innov")) {
        for (size_t k=0;k<n;++k) x[k]=integer_coefficient(c->dist,k);
        for (size_t i=1;i<m;++i) {
            memcpy(x+i*n,x+(i-1)*n+1,(n-1)*sizeof(*x));
            x[i*n+n-1]=integer_coefficient(c->dist,i+n-1);
            for (size_t h=0;h<c->q;++h) {
                const size_t k=(size_t)(next_random()%(uint64_t)(n-1));
                x[i*n+k]=integer_coefficient(c->dist,i*n+k);
            }
        }
    } else {
        int32_t *stream=malloc((m+n-1)*sizeof(*stream));
        if (!stream) return -1;
        for (size_t k=0;k<m+n-1;++k)
            stream[k]=integer_coefficient(c->dist,k);
        for (size_t i=0;i<m;++i)
            memcpy(x+i*n,stream+i,n*sizeof(*x));
        free(stream);
        if (c->q) {
            for (size_t i=0;i<m;++i) x[i*n+(17*i+7)%n]=0;
        }
    }
    return 0;
}

/* The checker is deliberately direct, O(m^2 n), and uses no Gram update. */
static int direct_integer_oracle(const struct test_case *c,const int32_t *x,
                                  int64_t *gram,uint64_t *max_abs_dot) {
    const uint64_t limit=UINT64_C(1)<<46;
    for (size_t p=0;p<c->m*c->n;++p)
        if (x[p]<-32768 || x[p]>32768) return -1;
    if ((uint64_t)c->n*UINT64_C(32768)*UINT64_C(32768)>limit) return -1;
    *max_abs_dot=0;
    for (size_t i=0;i<c->m;++i) {
        for (size_t j=0;j<=i;++j) {
            int64_t sum=0;
            for (size_t k=0;k<c->n;++k)
                sum+=(int64_t)x[i*c->n+k]*(int64_t)x[j*c->n+k];
            gram[i*c->m+j]=sum;
            const uint64_t magnitude=(uint64_t)(sum<0?-sum:sum);
            if (magnitude>limit) return -1;
            if (magnitude>*max_abs_dot) *max_abs_dot=magnitude;
        }
    }
    return 0;
}

struct method_result {
    const char *name;
    int status;
    size_t checked,nonfinite,mismatches,upper_or_padding_changed,input_changed;
    double maximum_absolute_error,maximum_normalized_error;
    size_t zero_norm_errors;
};
static int check_and_record(const struct test_case *c,size_t case_id,
                            const char *layout,struct method_result *r,
                            const int64_t *integer_gram,const double *g,
                            size_t ldg,FILE *events) {
    const int full=!strcmp(r->name,"accelerate_gemm");
    for (size_t i=0;i<c->m;++i) {
        for (size_t j=0;j<ldg;++j) {
            if (j>=c->m || (!full && j>i)) {
                if (!isnan(g[i*ldg+j])) ++r->upper_or_padding_changed;
                continue;
            }
            const size_t a=i>=j?i:j,b=i>=j?j:i;
            const int64_t integer=integer_gram[a*c->m+b];
            const int64_t norm_i=integer_gram[i*c->m+i];
            const int64_t norm_j=integer_gram[j*c->m+j];
            const double exact=scalbn((double)integer,2*c->e);
            const double got=g[i*ldg+j];
            ++r->checked;
            if (!isfinite(got)) {++r->nonfinite;continue;}
            if (got!=exact) {
                ++r->mismatches;
                /* All nonzero errors are exported losslessly. The Python
                 * receipt checks the preset tolerance with exact Fractions. */
                fprintf(events,"%zu\t%s\t%s\t%zu\t%zu\t%a\t%" PRId64
                               "\t%" PRId64 "\t%" PRId64 "\t%d\n",
                        case_id,layout,r->name,i,j,got,integer,norm_i,norm_j,c->e);
                const double error=fabs(got-exact);
                if (error>r->maximum_absolute_error) r->maximum_absolute_error=error;
                const double denominator=scalbn(sqrt((double)norm_i)*
                                                  sqrt((double)norm_j),2*c->e);
                if (denominator>0.0) {
                    const double ratio=error/denominator;
                    if (ratio>r->maximum_normalized_error) r->maximum_normalized_error=ratio;
                } else {++r->zero_norm_errors;}
            }
        }
    }
    return 0;
}
static void json_number(double x) {
    if (isfinite(x)) printf("%.17g",x); else printf("null");
}
static void print_method(const struct method_result *r,int first) {
    printf("%s{\"method\":\"%s\",\"status\":%d,\"entries_checked\":%zu,"
           "\"nonfinite\":%zu,\"exact_mismatches\":%zu,"
           "\"upper_or_padding_changed\":%zu,\"input_changed\":%zu,"
           "\"zero_norm_errors\":%zu,\"maximum_absolute_error\":",
           first?"":",",r->name,r->status,r->checked,r->nonfinite,r->mismatches,
           r->upper_or_padding_changed,r->input_changed,r->zero_norm_errors);
    json_number(r->maximum_absolute_error);
    printf(",\"maximum_normalized_error_approx\":");
    json_number(r->maximum_normalized_error);
    printf("}");
}
static int verify_layout(const struct test_case *c,size_t case_id,
                          const int32_t *x,const int64_t *oracle,
                          uint64_t maximum_dot,int padded,FILE *events) {
    const size_t lda=c->n+(padded?3u:0u),ldg=c->m+(padded?5u:0u);
    const char *layout=padded?"padded":"contiguous";
    double *a=malloc(c->m*lda*sizeof(*a)),*g=malloc(c->m*ldg*sizeof(*g));
    if (!a || !g) {free(a);free(g);return 1;}
    for (size_t i=0;i<c->m;++i) {
        for (size_t k=0;k<lda;++k)
            a[i*lda+k]=k<c->n?scalbn((double)x[i*c->n+k],c->e):NAN;
    }
    const char *names[3]={"lattice_sequential","accelerate_syrk","accelerate_gemm"};
    struct method_result results[3];
    int failed=0;
    for (size_t method=0;method<3;++method) {
        struct method_result *r=&results[method];
        memset(r,0,sizeof(*r));r->name=names[method];
        for (size_t p=0;p<c->m*ldg;++p) g[p]=NAN;
        if (method==0) {
            r->status=lattice_sequential_gram(a,c->m,c->n,lda,g,ldg,c->e,c->q);
        } else if (method==1) {
            cblas_dsyrk(CblasRowMajor,CblasLower,CblasNoTrans,(int)c->m,(int)c->n,
                        1.0,a,(int)lda,0.0,g,(int)ldg);
        } else {
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        (int)c->m,(int)c->m,(int)c->n,
                        1.0,a,(int)lda,a,(int)lda,0.0,g,(int)ldg);
        }
        check_and_record(c,case_id,layout,r,oracle,g,ldg,events);
        for (size_t i=0;i<c->m;++i) {
            for (size_t k=0;k<lda;++k) {
                if (k<c->n) {
                    if (a[i*lda+k]!=scalbn((double)x[i*c->n+k],c->e)) ++r->input_changed;
                } else if (!isnan(a[i*lda+k])) ++r->input_changed;
            }
        }
        if (r->status || r->nonfinite || r->upper_or_padding_changed ||
            r->input_changed || (!method && r->mismatches)) failed=1;
    }
    printf("{\"case_id\":%zu,\"family\":\"%s\",\"m\":%zu,\"n\":%zu,"
           "\"q_cap\":%zu,\"scale_exp\":%d,\"seed\":%" PRIu64
           ",\"dist\":%d,\"layout\":\"%s\",\"lda\":%zu,\"ldg\":%zu,"
           "\"oracle_max_abs_integer_dot\":%" PRIu64 ",\"methods\":[",
           case_id,c->family,c->m,c->n,c->q,c->e,c->seed,c->dist,layout,lda,ldg,maximum_dot);
    for (size_t method=0;method<3;++method) print_method(&results[method],method==0);
    printf("]}\n");fflush(stdout);
    free(a);free(g);return failed;
}
int main(int argc,char **argv) {
    if (argc!=2) {fprintf(stderr,"usage: algebra_lattice_full_entry EVENTS.tsv\n");return 2;}
    if (fegetround()!=FE_TONEAREST ||
        BLASSetThreading(BLAS_THREADING_SINGLE_THREADED) ||
        BLASGetThreading()!=BLAS_THREADING_SINGLE_THREADED) return 3;
    FILE *events=fopen(argv[1],"w");if (!events) return 4;
    fprintf(events,"case_id\tlayout\tmethod\ti\tj\tgot_hex\tinteger_target\tinteger_norm_i\tinteger_norm_j\tscale_exp\n");
    int failed=0;
    for (size_t id=0;id<sizeof(cases)/sizeof(cases[0]);++id) {
        const struct test_case *c=&cases[id];
        int32_t *x=malloc(c->m*c->n*sizeof(*x));
        int64_t *oracle=malloc(c->m*c->m*sizeof(*oracle));
        uint64_t maximum_dot=0;
        if (!x || !oracle || make_integer_input(c,x) ||
            direct_integer_oracle(c,x,oracle,&maximum_dot)) {
            free(x);free(oracle);fclose(events);return 5;
        }
        for (int padded=0;padded<=1;++padded)
            failed|=verify_layout(c,id,x,oracle,maximum_dot,padded,events);
        free(x);free(oracle);
    }
    if (fclose(events)) return 6;
    return failed?1:0;
}
