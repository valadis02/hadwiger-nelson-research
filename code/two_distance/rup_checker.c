/* Minimal forward DRUP checker (RUP lemmas only; deletions honoured).
   usage: rup formula.cnf proof.drup
   Every lemma must be RUP w.r.t. the current clause database; success iff the empty clause
   is derived (or the DB becomes conflicting under unit propagation). Deleting clauses only
   weakens the database, so honouring deletions is sound. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { int *lit; int n; int alive; } Clause;
static Clause *C; static int nC = 0, capC = 0;
static int nvars = 0;
static int **W; static int *Wn, *Wcap;      /* watch lists indexed by literal code */
static signed char *val;                    /* value per var: 0 unassigned, 1 true, -1 false */
static int *trail; static int tn = 0;

static int code(int l) { return l > 0 ? 2 * l : 2 * (-l) + 1; }
static int litval(int l) { int v = val[abs(l)]; return l > 0 ? v : -v; }
static void watch(int l, int ci) {
  int c = code(l);
  if (Wn[c] == Wcap[c]) { Wcap[c] = Wcap[c] ? 2 * Wcap[c] : 4; W[c] = realloc(W[c], sizeof(int) * Wcap[c]); }
  W[c][Wn[c]++] = ci;
}
static void ensure_vars(int v) {
  if (v <= nvars) return;
  int old = nvars; nvars = v * 2;
  val = realloc(val, nvars + 1); memset(val + old + 1, 0, nvars - old);
  trail = realloc(trail, sizeof(int) * (nvars + 1));
  W = realloc(W, sizeof(int *) * (2 * nvars + 2)); Wn = realloc(Wn, sizeof(int) * (2 * nvars + 2)); Wcap = realloc(Wcap, sizeof(int) * (2 * nvars + 2));
  for (int i = 2 * old + 2; i < 2 * nvars + 2; i++) { W[i] = NULL; Wn[i] = 0; Wcap[i] = 0; }
}
static int cmpint(const void *a, const void *b) { int x = *(int *)a, y = *(int *)b; return (x > y) - (x < y); }

/* hash of sorted clause for deletion lookup */
#define HB (1 << 22)
static int *hhead, *hnext;
static unsigned hash(int *l, int n) { unsigned h = 2166136261u; for (int i = 0; i < n; i++) { h ^= (unsigned)l[i]; h *= 16777619u; } return h & (HB - 1); }

static int add_clause(int *lits, int n) {
  if (nC == capC) { capC = capC ? 2 * capC : 1024; C = realloc(C, sizeof(Clause) * capC); hnext = realloc(hnext, sizeof(int) * capC); }
  int *l = malloc(sizeof(int) * (n ? n : 1)); memcpy(l, lits, sizeof(int) * n);
  qsort(l, n, sizeof(int), cmpint);
  C[nC].lit = l; C[nC].n = n; C[nC].alive = 1;
  unsigned h = hash(l, n); hnext[nC] = hhead[h]; hhead[h] = nC;
  if (n >= 2) { watch(l[0], nC); watch(l[1], nC); }
  return nC++;
}

/* unit propagation from scratch over current alive clauses; returns 1 on conflict.
   assumptions: neg of lemma literals already on trail. Units are handled by scanning unit clauses. */
static int *units; static int nunits = 0, capunits = 0;
static int propagate(int start) {
  int qh = start;
  while (qh < tn) {
    int l = trail[qh++];          /* l is true; visit clauses watching -l */
    int c = code(-l);
    for (int k = 0; k < Wn[c]; ) {
      int ci = W[c][k]; Clause *cl = &C[ci];
      if (!cl->alive) { W[c][k] = W[c][--Wn[c]]; continue; }
      int *L = cl->lit;
      if (L[0] == -l) { int t = L[0]; L[0] = L[1]; L[1] = t; }
      if (litval(L[0]) == 1) { k++; continue; }
      int found = 0;
      for (int m = 2; m < cl->n; m++) if (litval(L[m]) != -1) {
        int t = L[1]; L[1] = L[m]; L[m] = t; watch(L[1], ci); W[c][k] = W[c][--Wn[c]]; found = 1; break; }
      if (found) continue;
      if (litval(L[0]) == -1) return 1;
      val[abs(L[0])] = L[0] > 0 ? 1 : -1; trail[tn++] = L[0]; k++;
    }
  }
  return 0;
}
static void backtrack(int to) { while (tn > to) { int l = trail[--tn]; val[abs(l)] = 0; } }

static int base_units() {        /* assign all alive unit clauses at level 0; 1 if conflict */
  for (int u = 0; u < nunits; u++) {
    Clause *cl = &C[units[u]]; if (!cl->alive) continue;
    int l = cl->lit[0];
    if (litval(l) == -1) return 1;
    if (litval(l) == 0) { val[abs(l)] = l > 0 ? 1 : -1; trail[tn++] = l; }
  }
  return propagate(0);
}

static int read_clause(FILE *f, int *buf, int *n, int *isdel) {
  int x, k = 0; *isdel = 0; int c;
  while ((c = fgetc(f)) == ' ' || c == '\n' || c == '\r' || c == '\t');
  if (c == EOF) return 0;
  if (c == 'd') *isdel = 1; else ungetc(c, f);
  while (fscanf(f, "%d", &x) == 1 && x != 0) { buf[k++] = x; ensure_vars(abs(x)); }
  *n = k; return 1;
}

int main(int argc, char **argv) {
  hhead = malloc(sizeof(int) * HB); for (int i = 0; i < HB; i++) hhead[i] = -1;
  ensure_vars(16);
  FILE *f = fopen(argv[1], "r"); char line[256];
  int buf[100000], n, isdel;
  /* skip header */
  int c;
  while ((c = fgetc(f)) == 'c' || c == 'p') { if (!fgets(line, sizeof line, f)) break; }
  ungetc(c, f);
  long nform = 0;
  while (read_clause(f, buf, &n, &isdel)) {
    int ci = add_clause(buf, n); nform++;
    if (n == 1) { if (nunits == capunits) { capunits = capunits ? 2 * capunits : 64; units = realloc(units, sizeof(int) * capunits); } units[nunits++] = ci; }
    if (n == 0) { printf("formula contains empty clause\n"); return 0; }
  }
  fclose(f);
  printf("formula clauses: %ld\n", nform);
  FILE *p = fopen(argv[2], "r");
  long lem = 0, dels = 0, fail = 0;
  while (read_clause(p, buf, &n, &isdel)) {
    if (isdel) {
      int tmp[100000]; memcpy(tmp, buf, sizeof(int) * n); qsort(tmp, n, sizeof(int), cmpint);
      unsigned h = hash(tmp, n);
      for (int ci = hhead[h]; ci != -1; ci = hnext[ci]) {
        if (C[ci].alive && C[ci].n == n) {
          int *s = malloc(sizeof(int) * (n ? n : 1)); memcpy(s, C[ci].lit, sizeof(int) * n); qsort(s, n, sizeof(int), cmpint);
          int eq = !memcmp(s, tmp, sizeof(int) * n); free(s);
          if (eq && n != 1) { C[ci].alive = 0; dels++; break; }   /* never delete units (conservative, still sound) */
        }
      }
      continue;
    }
    lem++;
    backtrack(0);
    int conflict = base_units();
    if (!conflict) {
      int t0 = tn, triv = 0;
      for (int i = 0; i < n; i++) {
        int l = buf[i];
        if (litval(l) == 1) { triv = 1; break; }
        if (litval(l) == 0) { val[abs(l)] = l > 0 ? -1 : 1; trail[tn++] = -l; }
      }
      if (!triv) conflict = propagate(t0);
      else conflict = 1;
    }
    if (!conflict) { fail++; printf("lemma %ld NOT RUP\n", lem); return 1; }
    backtrack(0);
    int ci = add_clause(buf, n);
    if (n == 1) { if (nunits == capunits) { capunits = capunits ? 2 * capunits : 64; units = realloc(units, sizeof(int) * capunits); } units[nunits++] = ci; }
    if (n == 0) { printf("empty clause derived at lemma %ld (deletions applied: %ld)\nVERIFIED UNSAT\n", lem, dels); return 0; }
  }
  backtrack(0);
  if (base_units()) { printf("database conflicting after %ld lemmas (deletions: %ld)\nVERIFIED UNSAT\n", lem, dels); return 0; }
  printf("proof ended without contradiction (%ld lemmas)\nNOT VERIFIED\n", lem);
  return 1;
}
