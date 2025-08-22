
/*
 * flex_demux_mtx.c (v2)
 *
 * - Assignments use **1-based** indices into allowed_features.tsv (better for MTX/AnnData parity).
 * - Doublets are written to a separate file (<prefix>.doublets.txt).
 *
 * Build:
 *   gcc -O3 -march=native -o flex_demux_mtx flex_demux_mtx.c -lm
 *   # Optional parallel per-cell loop:
 *   # gcc -O3 -march=native -fopenmp -DUSE_OPENMP -o flex_demux_mtx flex_demux_mtx.c -lm
 */

#define _POSIX_C_SOURCE 200809L   // expose strdup, getline, etc.
#define _DEFAULT_SOURCE          // (glibc) fallback for older POSIX defs
#include <sys/types.h>
#include <sys/stat.h>            // stat, mkdir
#include <getopt.h>              // <── NEW
#include <errno.h>               // errno

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

static void die(const char *msg){ fprintf(stderr,"ERROR: %s\n", msg); exit(1); }
static char *trim(char *s){ if(!s) return s; while(isspace((unsigned char)*s)) s++; if(*s==0) return s; char *e=s+strlen(s)-1; while(e>s && isspace((unsigned char)*e)) *e--=0; return s; }

typedef struct { char **data; int n, cap; } VecS;
typedef struct { int *data; int n, cap; } VecI;
typedef struct { double *data; int n, cap; } VecD;
static void vecs_init(VecS *v){ v->data=NULL; v->n=0; v->cap=0; }
static void veci_init(VecI *v){ v->data=NULL; v->n=0; v->cap=0; }
static void vecd_init(VecD *v){ v->data=NULL; v->n=0; v->cap=0; }
static void vecs_push(VecS *v,const char*s){ if(v->n==v->cap){ v->cap=v->cap? v->cap*2:256; v->data=(char**)realloc(v->data,v->cap*sizeof(char*)); if(!v->data) die("OOM vecs"); } v->data[v->n++]=strdup(s?s:""); }
static void veci_push(VecI *v,int x){ if(v->n==v->cap){ v->cap=v->cap? v->cap*2:256; v->data=(int*)realloc(v->data,v->cap*sizeof(int)); if(!v->data) die("OOM veci"); } v->data[v->n++]=x; }
static void vecd_push(VecD *v,double x){ if(v->n==v->cap){ v->cap=v->cap? v->cap*2:256; v->data=(double*)realloc(v->data,v->cap*sizeof(double)); if(!v->data) die("OOM vecd"); } v->data[v->n++]=x; }
static int cmp_int(const void*a,const void*b){ int x=*(const int*)a,y=*(const int*)b; return (x<y)?-1:(x>y); }

/* String -> int map (open addressing) */
typedef struct { char **keys; int *vals; int cap; int n; } MapSI;
static unsigned long hash_str(const char *s){ unsigned long h=1469598103934665603ULL; while(*s){ h^=(unsigned char)(*s++); h*=1099511628211ULL; } return h; }
static void mapsi_init(MapSI *m,int cap){ m->cap=1; while(m->cap<cap*2) m->cap<<=1; m->keys=(char**)calloc(m->cap,sizeof(char*)); m->vals=(int*)malloc(m->cap*sizeof(int)); m->n=0; }
static void mapsi_put(MapSI *m,const char*key,int val){ if((m->n+1)*2>m->cap){ int oc=m->cap; char **ok=m->keys; int *ov=m->vals; m->cap<<=1; m->keys=(char**)calloc(m->cap,sizeof(char*)); m->vals=(int*)malloc(m->cap*sizeof(int)); m->n=0; for(int i=0;i<oc;++i){ if(ok[i]){ unsigned long h=hash_str(ok[i]); int j=(int)(h&(m->cap-1)); while(m->keys[j]) j=(j+1)&(m->cap-1); m->keys[j]=ok[i]; m->vals[j]=ov[i]; m->n++; } } free(ok); free(ov); }
    unsigned long h=hash_str(key); int i=(int)(h&(m->cap-1));
    while(m->keys[i]){ if(strcmp(m->keys[i],key)==0){ m->vals[i]=val; return; } i=(i+1)&(m->cap-1); }
    m->keys[i]=strdup(key); m->vals[i]=val; m->n++;
}
static int mapsi_get(const MapSI*m,const char*key){ unsigned long h=hash_str(key); int i=(int)(h&(m->cap-1)); while(m->keys[i]){ if(strcmp(m->keys[i],key)==0) return m->vals[i]; i=(i+1)&(m->cap-1);} return -1; }
static void mapsi_free(MapSI *m){ for(int i=0;i<m->cap;++i) if(m->keys[i]) free(m->keys[i]); free(m->keys); free(m->vals); }

/* Tails + BH */
static double poisson_tail_ge(int k,double lam){ if(k<=0) return 1.0; double pmf=exp(-lam),cdf=pmf; for(int i=1;i<k;++i){ pmf*=lam/(double)i; cdf+=pmf; if(1.0-cdf<1e-16) return 0.0; } double t=1.0-cdf; return (t<0?0:t); }
static double nb_tail_ge(int k,double r,double p){ if(k<=0) return 1.0; double pmf=pow(p,r),cdf=pmf; for(int x=0;x<k-1;++x){ pmf*=(r+x)*(1.0-p)/(double)(x+1); cdf+=pmf; if(1.0-cdf<1e-16) return 0.0; } double t=1.0-cdf; return (t<0?0:t); }
static double binom_tail_ge(int k,int n,double p){ if(k<=0) return 1.0; if(k>n) return 0.0; double q=1.0-p, pmf=pow(q,n), cdf=pmf; for(int i=0;i<k-1;++i){ pmf*=(double)(n-i)/(double)(i+1)*(p/q); cdf+=pmf; if(1.0-cdf<1e-16) return 0.0; } double t=1.0-cdf; return (t<0?0:t); }
typedef struct { int idx; double p; } IdxP;
static int cmp_idxp(const void*a,const void*b){ double x=((const IdxP*)a)->p,y=((const IdxP*)b)->p; return (x<y)?-1:(x>y); }
static void bh_fdr(const double *p,int n,double*q){ if(n<=0) return; IdxP*arr=(IdxP*)malloc(n*sizeof(IdxP)); for(int i=0;i<n;++i){ arr[i].idx=i; arr[i].p=p[i]; } qsort(arr,n,sizeof(IdxP),cmp_idxp); double prev=1.0; for(int r=1;r<=n;++r){ int i=arr[r-1].idx; double val=arr[r-1].p*n/(double)r; if(val<prev) prev=val; q[i]=prev<1.0?prev:1.0; } free(arr); }

/* Readers */
static void read_lines(const char*path,VecS*out){ FILE*f=fopen(path,"r"); if(!f){ fprintf(stderr,"Cannot open %s\n",path); die("open"); } char*line=NULL; size_t len=0; ssize_t r; while((r=getline(&line,&len,f))>=0){ char*s=trim(line); if(*s==0) continue; vecs_push(out,s); } if(line) free(line); fclose(f); }
static void read_features_firstcol(const char*path,VecS*feat_ids){ FILE*f=fopen(path,"r"); if(!f){ fprintf(stderr,"Cannot open %s\n",path); die("open"); } char*line=NULL; size_t len=0; ssize_t r; while((r=getline(&line,&len,f))>=0){ char*s=line; char*tab=strchr(s,'\t'); if(tab) *tab=0; s=trim(s); if(*s==0) continue; vecs_push(feat_ids,s);} if(line) free(line); fclose(f); }

/* Helper: load barcode list with better diagnostics */
static int load_barcode_list(const char *path, VecS *dest, const char *label) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "INFO: Cannot open %s for %s barcodes\n", path, label);
        return 0;  /* not found */
    }
    
    char *line = NULL;
    size_t len = 0;
    ssize_t r;
    int count = 0;
    
    while ((r = getline(&line, &len, f)) >= 0) {
        char *s = trim(line);
        if (*s == 0) continue;
        vecs_push(dest, s);
        count++;
    }
    
    if (line) free(line);
    fclose(f);
    
    fprintf(stderr, "INFO: Loaded %d %s barcodes from %s\n", count, label, path);
    return 1;  /* success */
}

/* Helper: try to load filtered barcodes from MTX directory */
static int try_load_filtered_from_mtx(const char *mtx_dir, VecS *out) {
    char path[1024];
    
    /* Try filtered_barcodes.tsv first */
    snprintf(path, sizeof(path), "%s/filtered_barcodes.tsv", mtx_dir);
    if (load_barcode_list(path, out, "MTX-derived filtered")) {
        return 1;
    }
    
    /* Try filtered_barcodes.txt second */
    snprintf(path, sizeof(path), "%s/filtered_barcodes.txt", mtx_dir);
    if (load_barcode_list(path, out, "MTX-derived filtered")) {
        return 1;
    }
    
    return 0;  /* not found */
}

/* Helper: mkdir -p equivalent - create directory and all parent directories */
static int mkdir_p(const char *path) {
    char tmp[1024];
    char *p = NULL;
    size_t len;
    
    snprintf(tmp, sizeof(tmp), "%s", path);
    len = strlen(tmp);
    if (tmp[len - 1] == '/') {
        tmp[len - 1] = 0;
    }
    
    for (p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = 0;
            if (mkdir(tmp, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0) {
                if (errno != EEXIST) {
                    fprintf(stderr, "WARNING: Cannot create directory %s: %s\n", tmp, strerror(errno));
                    return -1;
                }
            }
            *p = '/';
        }
    }
    
    if (mkdir(tmp, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0) {
        if (errno != EEXIST) {
            fprintf(stderr, "WARNING: Cannot create directory %s: %s\n", tmp, strerror(errno));
            return -1;
        }
    }
    
    return 0;
}

/* Helper: build output path - handles directory vs prefix mode */
static void build_output_path(const char *out_prefix, const char *suffix, char *buf, size_t buf_size) {
    struct stat st;
    int is_directory = 0;
    
    /* Check if out_prefix ends with '/' or is an existing directory */
    size_t len = strlen(out_prefix);
    if (len > 0 && out_prefix[len - 1] == '/') {
        is_directory = 1;
    } else if (stat(out_prefix, &st) == 0 && S_ISDIR(st.st_mode)) {
        is_directory = 1;
    }
    
    if (is_directory) {
        /* Directory mode: create directory if needed and place files inside */
        if (mkdir_p(out_prefix) != 0) {
            fprintf(stderr, "WARNING: Failed to create output directory %s\n", out_prefix);
        }
        
        /* Combine directory + suffix without leading dot */
        if (len > 0 && out_prefix[len - 1] == '/') {
            snprintf(buf, buf_size, "%s%s", out_prefix, suffix);
        } else {
            snprintf(buf, buf_size, "%s/%s", out_prefix, suffix);
        }
    } else {
        /* Prefix mode: traditional behavior with dot separator */
        snprintf(buf, buf_size, "%s.%s", out_prefix, suffix);
    }
}

typedef struct { char*barcode; int*counts; int total; } CellRec;

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    const char *mtx_dir = NULL, *solo_dir = NULL, *out_prefix = NULL;
    const char *cell_list_path = NULL;          /* <-- NEW */
    double tau = 0.8, delta = 0.4, gamma = 0.9, alpha = 1e-4, ambient_q = 0.999;
    int floor_ = 12, use_fdr = 1;

    /* ---------- getopt_long definitions ---------- */
    static const struct option long_opts[] = {
        {"mtx-dir",       required_argument, 0, 'm'},
        {"starsolo-dir",  required_argument, 0, 's'},
        {"out-prefix",    required_argument, 0, 'o'},
        {"cell-list",     required_argument, 0, 'c'},   /* NEW */
        {"tau",           required_argument, 0, 't'},
        {"delta",         required_argument, 0, 'd'},
        {"gamma",         required_argument, 0, 'g'},
        {"alpha",         required_argument, 0, 'a'},
        {"floor",         required_argument, 0, 'f'},
        {"ambient-q",     required_argument, 0, 'q'},
        {"no-fdr",        no_argument,       0, 'n'},
        {0, 0, 0, 0}
    };
    const char *short_opts = "m:s:o:c:t:d:g:a:f:q:n";    /* NEW ‘c:’ */

    int opt;
    while ((opt = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
        switch (opt) {
            case 'm': mtx_dir   = optarg;      break;
            case 's': solo_dir  = optarg;      break;
            case 'o': out_prefix= optarg;      break;
            case 'c': cell_list_path = optarg;      break;   /* NEW */
            case 't': tau       = atof(optarg);break;
            case 'd': delta     = atof(optarg);break;
            case 'g': gamma     = atof(optarg);break;
            case 'a': alpha     = atof(optarg);break;
            case 'f': floor_    = atoi(optarg);break;
            case 'q': ambient_q = atof(optarg);break;
            case 'n': use_fdr   = 0;           break;
            default:
                fprintf(stderr,
                    "Usage: %s --mtx-dir DIR [--starsolo-dir DIR] --out-prefix PRE "
                    "[--cell-list FILE] "
                    "[--tau X] [--delta X] [--gamma X] [--alpha X] [--floor N] "
                    "[--ambient-q X] [--no-fdr]\n", argv[0]);
                return 1;
        }
    }
    /* ------------------------------------------------ */

    if (!mtx_dir || !out_prefix) {
        fprintf(stderr,
            "ERROR: --mtx-dir and --out-prefix are required\n");
        return 1;
    }

    char path_barcodes[1024],path_features[1024],path_mtx[1024],path_allowed[1024];
    snprintf(path_barcodes,sizeof(path_barcodes), "%s/barcodes.tsv", mtx_dir);
    snprintf(path_features,sizeof(path_features), "%s/features.tsv", mtx_dir);
    snprintf(path_mtx,sizeof(path_mtx), "%s/matrix.mtx", mtx_dir);
    snprintf(path_allowed,sizeof(path_allowed), "%s/allowed_features.tsv", mtx_dir);
    char path_raw[1024],path_filt[1024];
    snprintf(path_raw,sizeof(path_raw), "%s/raw/barcodes.tsv", solo_dir);
    snprintf(path_filt,sizeof(path_filt), "%s/filtered/barcodes.tsv", solo_dir);

    /* MTX barcodes */
    VecS mtx_barcodes={0}; vecs_init(&mtx_barcodes);
    read_lines(path_barcodes,&mtx_barcodes);
    int n_cols=mtx_barcodes.n; if(n_cols<=0) die("MTX barcodes.tsv empty");
    MapSI mtx_col_idx; mapsi_init(&mtx_col_idx,n_cols*2+1);
    for(int j=0;j<n_cols;++j) mapsi_put(&mtx_col_idx, mtx_barcodes.data[j], j);

    /* Barcode resolution with priority order */
    VecS raw_barcodes={0}, filt_barcodes={0};
    vecs_init(&raw_barcodes); vecs_init(&filt_barcodes);

    /* 1. Load filtered barcodes with priority order */
    int filtered_source = 0;  /* 0=none, 1=cell-list, 2=mtx-derived, 3=starsolo */
    
    if (cell_list_path) {
        /* Priority 1: --cell-list provided by user */
        if (load_barcode_list(cell_list_path, &filt_barcodes, "user-provided filtered")) {
            filtered_source = 1;
        }
    }
    
    if (filtered_source == 0) {
        /* Priority 2: MTX-derived filtered list */
        if (try_load_filtered_from_mtx(mtx_dir, &filt_barcodes)) {
            filtered_source = 2;
        }
    }
    
    if (filtered_source == 0 && solo_dir) {
        /* Priority 3: STARsolo fallback */
        if (load_barcode_list(path_filt, &filt_barcodes, "STARsolo filtered")) {
            filtered_source = 3;
        }
    }
    
    if (filtered_source == 0) {
        fprintf(stderr, "ERROR: No filtered barcode list found. Please either:\n");
        fprintf(stderr, "  1. Provide --cell-list FILE, or\n");
        fprintf(stderr, "  2. Add filtered_barcodes.tsv to your MTX directory, or\n");
        fprintf(stderr, "  3. Supply --starsolo-dir with filtered/barcodes.tsv\n");
        return 1;
    }

    /* 2. Load raw barcodes (for ambient estimation) */
    if (solo_dir) {
        /* Use STARsolo raw list when available */
        load_barcode_list(path_raw, &raw_barcodes, "STARsolo raw");
    } else {
        /* Fallback: use all MTX barcodes as raw */
        fprintf(stderr, "INFO: No --starsolo-dir provided; using all MTX barcodes as raw set\n");
        for (int i = 0; i < mtx_barcodes.n; ++i) {
            vecs_push(&raw_barcodes, mtx_barcodes.data[i]);
        }
    }
    int n_filtered = filt_barcodes.n;

    /* sets */
    VecI cells_cols={0}; veci_init(&cells_cols);
    VecI negs_cols={0}; veci_init(&negs_cols);
    int missing_cells=0;
    MapSI filt_map; mapsi_init(&filt_map, filt_barcodes.n*2+1);
    for(int i=0;i<filt_barcodes.n;++i) mapsi_put(&filt_map, filt_barcodes.data[i], 1);
    for(int i=0;i<filt_barcodes.n;++i){
        int col=mapsi_get(&mtx_col_idx, filt_barcodes.data[i]);
        if(col>=0) veci_push(&cells_cols,col); else missing_cells++;
    }
    for(int i=0;i<raw_barcodes.n;++i){
        if(mapsi_get(&filt_map, raw_barcodes.data[i])==1) continue;
        int col=mapsi_get(&mtx_col_idx, raw_barcodes.data[i]);
        if(col>=0) veci_push(&negs_cols,col);
    }

    /* features + allowlist */
    VecS feat_ids={0}; vecs_init(&feat_ids); read_features_firstcol(path_features,&feat_ids);
    VecS allowed={0}; vecs_init(&allowed); read_lines(path_allowed,&allowed);
    int K=0; int *row_to_allowed=(int*)malloc(feat_ids.n*sizeof(int));
    for(int i=0;i<feat_ids.n;++i) row_to_allowed[i]=-1;
    MapSI allowed_map; mapsi_init(&allowed_map, allowed.n*2+1);
    for(int i=0;i<allowed.n;++i) mapsi_put(&allowed_map, allowed.data[i], i+1); /* 1-based */
    for(int r=0;r<feat_ids.n;++r){ int idx=mapsi_get(&allowed_map, feat_ids.data[r]); if(idx>=1){ row_to_allowed[r]=idx; if(idx>K) K=idx; } }
    if(K==0) die("No allowed features found");

    int *col_status=(int*)calloc(n_cols,sizeof(int)); /* 0 ignore, 1 cell, 2 neg */
    for(int t=0;t<cells_cols.n;++t) col_status[cells_cols.data[t]]=1;
    for(int t=0;t<negs_cols.n;++t) col_status[negs_cols.data[t]]=2;

    double *ambient_colsum=(double*)calloc(K+1,sizeof(double)); /* 1..K */
    int *neg_totals=(int*)calloc(n_cols,sizeof(int));

    int num_cells=cells_cols.n;
    typedef struct { char*barcode; int *counts; int total; } CellRec;
    CellRec *cells=(CellRec*)calloc(num_cells,sizeof(CellRec));
    int *col_to_cellidx=(int*)malloc(n_cols*sizeof(int));
    for(int j=0;j<n_cols;++j) col_to_cellidx[j]=-1;
    for(int i=0;i<num_cells;++i){
        int col=cells_cols.data[i];
        cells[i].barcode=strdup(mtx_barcodes.data[col]);
        cells[i].counts=(int*)calloc(K+1,sizeof(int)); /* 1..K */
        cells[i].total=0; col_to_cellidx[col]=i;
    }

    /* read matrix.mtx */
    FILE *fm=fopen(path_mtx,"r"); if(!fm){ fprintf(stderr,"Cannot open %s\n",path_mtx); return 1; }
    char *line=NULL; size_t len=0; ssize_t rr;
    if((rr=getline(&line,&len,fm))<0) die("matrix.mtx empty");
    if(strncmp(line,"%%MatrixMarket",14)!=0){
        fprintf(stderr, "%s is missing %%MatrixMarket\n", path_mtx);
        return 1;
    }
    while((rr=getline(&line,&len,fm))>=0){
        char*s=trim(line); if(*s=='%') continue;
        int n_rows=0,n_cols2=0,nnz=0;
        if(sscanf(s,"%d %d %d",&n_rows,&n_cols2,&nnz)!=3) die("Bad size line");
        if(n_cols2!=n_cols) fprintf(stderr,"WARNING: barcodes.tsv columns (%d) != mtx columns (%d)\n", n_cols, n_cols2);
        int r_i,c_j; double v;
        for(int t=0;t<nnz;++t){
            if((rr=getline(&line,&len,fm))<0) die("Unexpected EOF");
            if(sscanf(line,"%d %d %lf",&r_i,&c_j,&v)!=3) die("Bad triplet");
            int row0=r_i-1, col0=c_j-1;
            if(row0<0||row0>=feat_ids.n||col0<0||col0>=n_cols) continue;
            int aidx=row_to_allowed[row0]; if(aidx<1) continue; /* not allowed */
            int status=col_status[col0];
            if(status==2){ ambient_colsum[aidx]+=v; neg_totals[col0]+= (int)v; }
            else if(status==1){ int ci=col_to_cellidx[col0]; if(ci>=0){ cells[ci].counts[aidx]+= (int)v; cells[ci].total+= (int)v; } }
        }
        break;
    }
    /* close input */
    if (line)
        free(line);
    fclose(fm);

    /* ambient vector r */
    VecD rvec={0}; vecd_init(&rvec); double sum_amb=0.0;
    for(int j=1;j<=K;++j) sum_amb+=ambient_colsum[j];
    if(sum_amb<=0){ for(int j=1;j<=K;++j) vecd_push(&rvec, 1.0/(double)K); }
    else { for(int j=1;j<=K;++j) vecd_push(&rvec, ambient_colsum[j]/sum_amb); }

    /* negative totals vector */
    VecI neg_vec={0}; veci_init(&neg_vec);
    for(int j=0;j<n_cols;++j) if(col_status[j]==2) veci_push(&neg_vec, neg_totals[j]);

    /* thresholds */
    int Nneg=neg_vec.n;
    int *neg_sorted=(int*)malloc(Nneg*sizeof(int)); memcpy(neg_sorted, neg_vec.data, Nneg*sizeof(int)); qsort(neg_sorted,Nneg,sizeof(int),cmp_int);
    int cut=(int)floor(0.7*Nneg); if(cut<1) cut=1;
    double lam_hat=0.0; for(int i=0;i<cut;++i) lam_hat+=neg_sorted[i]; lam_hat/=(double)cut;
    double p_cut=fmin(1e-6, 0.1/fmax(1,Nneg));
    VecI kept={0}; veci_init(&kept);
    for(int i=0;i<Nneg;++i){ double p=poisson_tail_ge(neg_vec.data[i], lam_hat); if(p>=p_cut) veci_push(&kept, neg_vec.data[i]); }
    int Kk=kept.n; int *ks=(int*)malloc(Kk*sizeof(int)); memcpy(ks, kept.data, Kk*sizeof(int)); qsort(ks,Kk,sizeof(int),cmp_int);
    int lo_i=(int)floor(0.10*Kk), hi_i=(int)floor(0.80*Kk); if(hi_i<=lo_i){ lo_i=0; hi_i=Kk; }
    double mu = 0.0, var = 0.0;
    int len2 = hi_i - lo_i;
    if (len2 <= 0) { lo_i = 0; hi_i = Kk; len2 = Kk; }

    for (int i = lo_i; i < hi_i; ++i)
        mu += ks[i];
    mu /= (double)fmax(1, len2);

    for(int i=lo_i;i<hi_i;++i){ double d=ks[i]-mu; var+=d*d; } var/=(double)fmax(1,len2-1);
    int over=(var>1.2*mu);
    int m_emp=0; if(Kk>0){ int idx=(int)floor((Kk-1)*ambient_q); if(idx<0) idx=0; if(idx>=Kk) idx=Kk-1; m_emp=ks[idx]; }
    double alpha_amb=1e-3; int m_mod=0;
    if(!over || var<=mu){ int m=(int)fmax(0.0,floor(mu)); while(poisson_tail_ge(m,mu)<=alpha_amb && m>0) m--; while(poisson_tail_ge(m,mu)>alpha_amb) m++; m_mod=m; }
    else { double k_nb=(mu*mu)/fmax(1e-12,(var-mu)); double p_nb=k_nb/(k_nb+mu); int m=(int)fmax(0.0,floor(mu)); while(nb_tail_ge(m,k_nb,p_nb)<=alpha_amb && m>0) m--; while(nb_tail_ge(m,k_nb,p_nb)>alpha_amb) m++; m_mod=m; }
    int M_ambient=(m_emp>m_mod? m_emp:m_mod);
    double rmax=0.0; for(int j=0;j<rvec.n;++j) if(rvec.data[j]>rmax) rmax=rvec.data[j];
    int M_stat=1; for(int m=1;m<=1000;++m){ int kth=(int)floor(tau*m); if(kth<1) kth=1; double pt=binom_tail_ge(kth,m,rmax); if(pt<alpha){ M_stat=m; break; } M_stat=m; }
    int M_min=M_ambient; if(M_stat>M_min) M_min=M_stat; if(floor_>M_min) M_min=floor_;
    M_min=6;
    M_ambient=10;
    printf("M_min: %d\n", M_min);
    printf("M_ambient: %d\n", M_ambient);
    printf("M_stat: %d\n", M_stat);
    printf("floor_: %d\n", floor_);
    printf("over: %d\n", over);
    printf("var: %f\n", var);
    printf("mu: %f\n", mu);
    
    /* per-cell stats */
    int C=num_cells;
    double *p1=(double*)malloc(C*sizeof(double)), *p2=(double*)malloc(C*sizeof(double));
    double *q1=(double*)malloc(C*sizeof(double)), *q2=(double*)malloc(C*sizeof(double));
    int *j1=(int*)malloc(C*sizeof(int)), *j2=(int*)malloc(C*sizeof(int));
    double *f1=(double*)malloc(C*sizeof(double)), *f2=(double*)malloc(C*sizeof(double));
    int *m_tot=(int*)malloc(C*sizeof(int)), *c1=(int*)malloc(C*sizeof(int)), *c2=(int*)malloc(C*sizeof(int));

    for(int i=0;i<C;++i){
        m_tot[i]=cells[i].total;
        if(m_tot[i]<M_min){ p1[i]=1.0; p2[i]=1.0; j1[i]=j2[i]=-1; f1[i]=f2[i]=0.0; c1[i]=c2[i]=0; continue; }
        int t1=-1,t2=-1;
        for(int kx=1;kx<=K;++kx){
            if(t1<0 || cells[i].counts[kx] > cells[i].counts[t1]){ t2=t1; t1=kx; }
            else if(t2<0 || cells[i].counts[kx] > cells[i].counts[t2]){ t2=kx; }
        }
        j1[i]=t1; j2[i]=t2;
        c1[i]=(t1>=1)? cells[i].counts[t1]:0;
        c2[i]=(t2>=1)? cells[i].counts[t2]:0;
        f1[i]=(m_tot[i]>0)? (double)c1[i]/(double)m_tot[i] : 0.0;
        f2[i]=(m_tot[i]>0)? (double)c2[i]/(double)m_tot[i] : 0.0;
        p1[i]=(t1>=1)? binom_tail_ge(c1[i], m_tot[i], rvec.data[t1-1]) : 1.0;
        p2[i]=(t2>=1 && c2[i]>0)? binom_tail_ge(c2[i], m_tot[i], rvec.data[t2-1]) : 1.0;
    }
    if(use_fdr){ bh_fdr(p1,C,q1); bh_fdr(p2,C,q2); } else { for(int i=0;i<C;++i){ q1[i]=p1[i]; q2[i]=p2[i]; } }

    /* outputs */
    char path_assign[1024],path_dbl[1024],path_amb[1024],path_unas[1024],path_miss[1024];
    build_output_path(out_prefix, "assignments.tsv", path_assign, sizeof(path_assign));
    build_output_path(out_prefix, "doublets.txt", path_dbl, sizeof(path_dbl));
    build_output_path(out_prefix, "ambiguous.txt", path_amb, sizeof(path_amb));
    build_output_path(out_prefix, "unassignable.txt", path_unas, sizeof(path_unas));
    build_output_path(out_prefix, "missing_cells.txt", path_miss, sizeof(path_miss));
    FILE *fa=fopen(path_assign,"w"); if(!fa) die("open assignments");
    FILE *fd=fopen(path_dbl,"w"); if(!fd) die("open doublets");
    FILE *fb=fopen(path_amb,"w"); if(!fb) die("open ambiguous");
    FILE *fu=fopen(path_unas,"w"); if(!fu) die("open unassignable");
    FILE *fmz=fopen(path_miss,"w"); if(!fmz) die("open missing");

    for(int i=0;i<filt_barcodes.n;++i){
        if(mapsi_get(&mtx_col_idx, filt_barcodes.data[i]) < 0) fprintf(fmz,"%s\n",filt_barcodes.data[i]);
    }
    fclose(fmz);

    int n_singlet=0,n_doublet=0,n_ambig=0,n_low=0;
    for(int i=0;i<C;++i){
        if(m_tot[i] < M_min){ fprintf(fu,"%s\n", cells[i].barcode); n_low++; continue; }
        double top_sum=f1[i]+f2[i];
        double frac1=(top_sum>0)? f1[i]/top_sum : 0.0;
        int is_singlet=(f1[i]>=tau) && ((f1[i]-f2[i])>=delta) && (q1[i] < alpha);
        int is_doublet=(top_sum>=gamma) && (frac1>=0.2 && frac1<=0.8) && (q1[i]<alpha) && (q2[i]<alpha);

        if(is_singlet && j1[i]>=1){ fprintf(fa, "%s\t%d\n", cells[i].barcode, j1[i]); n_singlet++; }
        else if(is_doublet){ fprintf(fd, "%s\n", cells[i].barcode); n_doublet++; }
        else { fprintf(fb, "%s\n", cells[i].barcode); n_ambig++; }
    }
    fclose(fa); fclose(fd); fclose(fb); fclose(fu);

    /* QC */
    int cells_in_mtx=cells_cols.n;
    printf("Filtered STARsolo cells: %d\n", n_filtered);
    printf("Cells present in MTX: %d  (missing: %d)\n", cells_in_mtx, missing_cells);
    int total_eval=n_singlet+n_doublet+n_ambig+n_low;
    printf("Evaluated (filtered∩MTX): %d\n", total_eval);
    printf("Singlet: %d (%.3f)\n", n_singlet, (total_eval? (double)n_singlet/total_eval : 0.0));
    printf("Doublet: %d (%.3f)\n", n_doublet, (total_eval? (double)n_doublet/total_eval : 0.0));
    printf("Ambiguous: %d (%.3f)\n", n_ambig, (total_eval? (double)n_ambig/total_eval : 0.0));
    printf("Low support: %d (%.3f)\n", n_low, (total_eval? (double)n_low/total_eval : 0.0));
    printf("M_min=%d (M_ambient=%d, M_stat=%d, floor=%d), overdispersed=%s, r_max=%.4f\n",
           M_min, M_ambient, M_stat, floor_, (over? "yes":"no"), rmax);

    /* cleanup */
    mapsi_free(&mtx_col_idx);
    mapsi_free(&filt_map);
    mapsi_free(&allowed_map);

    for (int i = 0; i < mtx_barcodes.n; ++i)
        free(mtx_barcodes.data[i]);
    free(mtx_barcodes.data);

    for (int i = 0; i < raw_barcodes.n; ++i)
        free(raw_barcodes.data[i]);
    free(raw_barcodes.data);

    for (int i = 0; i < filt_barcodes.n; ++i)
        free(filt_barcodes.data[i]);
    free(filt_barcodes.data);

    for (int i = 0; i < feat_ids.n; ++i)
        free(feat_ids.data[i]);
    free(feat_ids.data);

    for (int i = 0; i < allowed.n; ++i)
        free(allowed.data[i]);
    free(allowed.data);

    free(row_to_allowed);
    free(col_status);
    free(ambient_colsum);
    free(neg_totals);
    free(neg_sorted);
    free(ks);

    for (int i = 0; i < num_cells; ++i) {
        free(cells[i].barcode);
        free(cells[i].counts);
    }
    free(cells);

    free(p1); free(p2); free(q1); free(q2);
    free(j1); free(j2); free(f1); free(f2);
    free(m_tot); free(c1); free(c2);
    free(rvec.data);

    return 0;
}
