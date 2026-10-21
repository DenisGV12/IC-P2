/*
 * solver3d_mpi.c — solver de fluidos 3D de Jos Stam, versión MPI (Práctica 3, IC, UA 2026-27)
 *
 * Basado en solver3d.c (práctica 2). Estrategia:
 *   - Reparto de la malla en planos k contiguos entre procesos (slabs, ver slab.h).
 *   - lin_solve: Gauss-Seidel LEXICOGRÁFICO (original) -> ROJO-NEGRO. El original actualiza
 *     in-place con dependencia (i-1,j-1,k-1) ya actualizados, imposible de repartir. En rojo-negro
 *     las celdas de un color solo dependen del otro, así que cada medio barrido es paralelo y solo
 *     requiere intercambiar los planos frontera (halo) entre vecinos. El resultado es idéntico
 *     para cualquier número de procesos (cada celda depende únicamente de sus vecinos).
 *   - advect: semi-lagrangiano, lee d0 en posiciones arbitrarias -> se replica d0 con MPI_Allgatherv
 *     y cada proceso calcula solo sus planos.
 *   - project/add_source/set_bnd: locales a cada slab (set_bnd sin comunicación; ver comentarios).
 *
 * Cambios respecto a solver3d.c marcados con [mod P3].
 */
#include <stdlib.h>
#include "slab.h"

#define IX(i,j,k) ((i)+(M+2)*(j) + (M+2)*(N+2)*(k))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define MAX(a,b)            (((a) > (b)) ? (a) : (b))
#ifndef LINEARSOLVERTIMES
#define LINEARSOLVERTIMES 10
#endif

Slab sl;

/* ------------------------------------------------------------------ */
/* Infraestructura de reparto y comunicación                          */
/* ------------------------------------------------------------------ */

int slab_init ( int M, int N, int O )
{
	MPI_Comm_rank ( MPI_COMM_WORLD, &sl.rank );
	MPI_Comm_size ( MPI_COMM_WORLD, &sl.nprocs );
	if ( O < sl.nprocs ) return 0;                       /* al menos un plano interior por proceso */

	int base = O / sl.nprocs, rem = O % sl.nprocs;       /* reparto equilibrado: ±1 plano */
	sl.klo = 1 + sl.rank*base + (sl.rank < rem ? sl.rank : rem);
	sl.khi = sl.klo + base + (sl.rank < rem ? 1 : 0) - 1;
	sl.own_lo = sl.rank == 0              ? 0   : sl.klo;
	sl.own_hi = sl.rank == sl.nprocs - 1  ? O+1 : sl.khi;
	sl.plane  = (size_t)(M+2)*(N+2);

	sl.counts = malloc ( sl.nprocs * sizeof(int) );
	sl.displs = malloc ( sl.nprocs * sizeof(int) );
	if ( !sl.counts || !sl.displs ) return 0;
	for ( int r = 0 ; r < sl.nprocs ; r++ ) {
		int lo = 1 + r*base + (r < rem ? r : rem);
		int hi = lo + base + (r < rem ? 1 : 0) - 1;
		int olo = r == 0 ? 0 : lo, ohi = r == sl.nprocs-1 ? O+1 : hi;
		sl.displs[r] = (int)(olo * sl.plane);
		sl.counts[r] = (int)((ohi - olo + 1) * sl.plane);
	}
	sl.t_halo = sl.t_gather = sl.bytes_halo = sl.bytes_gather = 0.0;
	sl.n_halo = sl.n_gather = 0;
	return 1;
}

void slab_free ( void ) { free ( sl.counts ); free ( sl.displs ); }

/*
 * Intercambio de halos: cada vecino envía su plano frontera interior y recibe el del vecino
 * en el plano fantasma. Dos MPI_Sendrecv (hacia arriba y hacia abajo); MPI_PROC_NULL en los
 * extremos, donde el plano fantasma es el contorno del dominio y lo gestiona set_bnd.
 */
void halo_exchange ( float * x )
{
	if ( sl.nprocs == 1 ) return;
	double t0 = MPI_Wtime ();
	int up   = sl.rank < sl.nprocs-1 ? sl.rank+1 : MPI_PROC_NULL;
	int down = sl.rank > 0           ? sl.rank-1 : MPI_PROC_NULL;
	int n = (int)sl.plane;
	/* hacia arriba: mi último plano -> ghost inferior del vecino de arriba */
	MPI_Sendrecv ( x + sl.khi*sl.plane,     n, MPI_FLOAT, up,   0,
	               x + (sl.klo-1)*sl.plane, n, MPI_FLOAT, down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	/* hacia abajo: mi primer plano -> ghost superior del vecino de abajo */
	MPI_Sendrecv ( x + sl.klo*sl.plane,     n, MPI_FLOAT, down, 1,
	               x + (sl.khi+1)*sl.plane, n, MPI_FLOAT, up,   1, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	sl.t_halo += MPI_Wtime () - t0;
	sl.n_halo++;
	sl.bytes_halo += 2.0 * n * sizeof(float);   /* enviados por este proceso */
}

/* Estrategia B (overlap=1): intercambio no bloqueante; el cálculo de los planos interiores se solapa con la comunicación. */
void halo_start ( float * x, MPI_Request rq[4] )
{
	double t0 = MPI_Wtime ();
	int up   = sl.rank < sl.nprocs-1 ? sl.rank+1 : MPI_PROC_NULL;
	int down = sl.rank > 0           ? sl.rank-1 : MPI_PROC_NULL;
	int n = (int)sl.plane;
	MPI_Irecv ( x + (sl.klo-1)*sl.plane, n, MPI_FLOAT, down, 0, MPI_COMM_WORLD, &rq[0] );
	MPI_Irecv ( x + (sl.khi+1)*sl.plane, n, MPI_FLOAT, up,   1, MPI_COMM_WORLD, &rq[1] );
	MPI_Isend ( x + sl.khi*sl.plane,     n, MPI_FLOAT, up,   0, MPI_COMM_WORLD, &rq[2] );
	MPI_Isend ( x + sl.klo*sl.plane,     n, MPI_FLOAT, down, 1, MPI_COMM_WORLD, &rq[3] );
	sl.t_halo += MPI_Wtime () - t0;
	sl.n_halo++;
	sl.bytes_halo += 2.0 * n * sizeof(float);
}

void halo_finish ( MPI_Request rq[4] )
{
	double t0 = MPI_Wtime ();
	MPI_Waitall ( 4, rq, MPI_STATUSES_IGNORE );
	sl.t_halo += MPI_Wtime () - t0;
}

/* Replica el campo completo (planos propios de cada proceso -> todos). */
void allgather_field ( float * x )
{
	if ( sl.nprocs == 1 ) return;
	double t0 = MPI_Wtime ();
	MPI_Allgatherv ( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, x, sl.counts, sl.displs, MPI_FLOAT, MPI_COMM_WORLD );
	sl.t_gather += MPI_Wtime () - t0;
	sl.n_gather++;
	sl.bytes_gather += (double)sl.counts[sl.rank] * sizeof(float);
}

/* ------------------------------------------------------------------ */
/* Solver                                                             */
/* ------------------------------------------------------------------ */

void add_source ( int M, int N, int O, float * x, float * s, float dt )
{
	(void)M; (void)N; (void)O;
	size_t i, lo = (size_t)sl.own_lo * sl.plane, hi = ((size_t)sl.own_hi + 1) * sl.plane;   /* [mod P3] solo mis planos */
	for ( i=lo ; i<hi ; i++ ) x[i] += dt*s[i];
}

/*
 * [mod P3] set_bnd sin comunicación: cada proceso rellena solo el contorno que le pertenece.
 *  - caras z (k=0, k=O+1) y bordes/esquinas en esos planos: rank 0 / último rank
 *  - caras x e y y bordes paralelos a k: los planos propios klo..khi
 * Todos los valores de contorno se calculan a partir de caras interiores de los mismos planos (o de
 * los planos k=1 / k=O, que pertenecen al rank 0 / al último), por lo que no hace falta comunicar.
 */
void  set_bnd ( int M, int N, int O, int b, float * x )
{
	int i, j;
	const int first = sl.rank == 0, last = sl.rank == sl.nprocs-1;

	/* caras z */
	if ( first ) for ( i=1 ; i<=M ; i++ ) for ( j=1 ; j<=N ; j++ )
		x[IX(i,j,0 )] = b==3 ? -x[IX(i,j,1)] : x[IX(i,j,1)];
	if ( last )  for ( i=1 ; i<=M ; i++ ) for ( j=1 ; j<=N ; j++ )
		x[IX(i,j,O+1)] = b==3 ? -x[IX(i,j,O)] : x[IX(i,j,O)];

	/* caras x e y, planos propios */
	for ( j=sl.klo ; j<=sl.khi ; j++ ) {
		for ( i=1 ; i<=N ; i++ ) {
			x[IX(0  ,i, j)] = b==1 ? -x[IX(1,i,j)] : x[IX(1,i,j)];
			x[IX(M+1,i, j)] = b==1 ? -x[IX(M,i,j)] : x[IX(M,i,j)];
		}
		for ( i=1 ; i<=M ; i++ ) {
			x[IX(i,0,j )] = b==2 ? -x[IX(i,1,j)] : x[IX(i,1,j)];
			x[IX(i,N+1,j)] = b==2 ? -x[IX(i,N,j)] : x[IX(i,N,j)];
		}
	}

	/* bordes paralelos a x (k=0 y k=O+1) */
	if ( first ) for ( i=1; i<=M; i++) {
		x[IX(i,  0,  0)] = 1.0/2.0*(x[IX(i,1,  0)]+x[IX(i,  0,  1)]);
		x[IX(i,N+1,  0)] = 1.0/2.0*(x[IX(i,N,  0)]+x[IX(i,N+1,  1)]);
	}
	if ( last ) for ( i=1; i<=M; i++) {
		x[IX(i,  0,O+1)] = 1.0/2.0*(x[IX(i,0,  O)]+x[IX(i,  1,O+1)]);
		x[IX(i,N+1,O+1)] = 1.0/2.0*(x[IX(i,N,O+1)]+x[IX(i,N+1,  O)]);
	}

	/* bordes paralelos a y (k=0 y k=O+1) */
	if ( first ) for ( i=1; i<=N; i++) {
		x[IX(0,  i,  0)] = 1.0/2.0*(x[IX(1,i,  0)]+x[IX(0,  i,  1)]);
		x[IX(M+1,i,  0)] = 1.0/2.0*(x[IX(M,i,  0)]+x[IX(M+1,i,  1)]);
	}
	if ( last ) for ( i=1; i<=N; i++) {
		x[IX(0,  i,O+1)] = 1.0/2.0*(x[IX(0,i,  O)]+x[IX(1,  i,O+1)]);
		x[IX(M+1,i,O+1)] = 1.0/2.0*(x[IX(M,i,O+1)]+x[IX(M+1,i,  O)]);
	}

	/* bordes paralelos a z, planos propios */
	for ( i=sl.klo; i<=sl.khi; i++) {
		x[IX(0,  0,  i)] = 1.0/2.0*(x[IX(0,  1,i)]+x[IX(1,  0,  i)]);
		x[IX(0,  N+1,i)] = 1.0/2.0*(x[IX(0,  N,i)]+x[IX(1,  N+1,i)]);
		x[IX(M+1,0,  i)] = 1.0/2.0*(x[IX(M,  0,i)]+x[IX(M+1,1,  i)]);
		x[IX(M+1,N+1,i)] = 1.0/2.0*(x[IX(M+1,N,i)]+x[IX(M,  N+1,i)]);
	}

	/* esquinas */
	if ( first ) {
		x[IX(0  ,0, 0  )] = 1.0/3.0*(x[IX(1,0,0  )]+x[IX(0  ,1,0)]+x[IX(0 ,0,1)]);
		x[IX(0  ,N+1, 0)] = 1.0/3.0*(x[IX(1,N+1, 0)]+x[IX(0  ,N, 0)] + x[IX(0  ,N+1, 1)]);
		x[IX(M+1,0, 0 )] = 1.0/3.0*(x[IX(M,0,0  )]+x[IX(M+1,1,0)] + x[IX(M+1,0,1)]) ;
		x[IX(M+1,N+1,0)] = 1.0/3.0*(x[IX(M,N+1,0)]+x[IX(M+1,N,0)]+x[IX(M+1,N+1,1)]);
	}
	if ( last ) {
		x[IX(0  ,0, O+1 )] = 1.0/3.0*(x[IX(1,0,O+1  )]+x[IX(0  ,1,O+1)]+x[IX(0 ,0,O)]);
		x[IX(0  ,N+1, O+1)] = 1.0/3.0*(x[IX(1,N+1, O+1)]+x[IX(0  ,N, O+1)] + x[IX(0  ,N+1, O)]);
		x[IX(M+1,0, O+1 )] = 1.0/3.0*(x[IX(M,0,O+1  )]+x[IX(M+1,1,O+1)] + x[IX(M+1,0,O)]) ;
		x[IX(M+1,N+1,O+1)] = 1.0/3.0*(x[IX(M,N+1,O+1)]+x[IX(M+1,N,O+1)]+x[IX(M+1,N+1,O)]);
	}
}

/*
 * [mod P3] Gauss-Seidel rojo-negro por slabs.
 * Por iteración: barrido de celdas con (i+j+k) par, intercambio de halos, barrido de impares,
 * intercambio de halos, set_bnd. Bucles en orden k-j-i (i contiguo en memoria, paso 2 por color).
 */
/* medio barrido de un color sobre los planos k0..k1 (orden k-j-i, i contiguo, paso 2) */
static void sweep ( int M, int N, float * x, float * x0, float a, float c, int color, int k0, int k1 )
{
	int i, j, k;
	for ( k=k0 ; k<=k1 ; k++ ) { for ( j=1 ; j<=N ; j++ ) {
		for ( i=1+((j+k+1+color)&1) ; i<=M ; i+=2 ) {
			x[IX(i,j,k)] = (x0[IX(i,j,k)] + a*(x[IX(i-1,j,k)]+x[IX(i+1,j,k)]+x[IX(i,j-1,k)]+x[IX(i,j+1,k)]+x[IX(i,j,k-1)]+x[IX(i,j,k+1)]))/c;
		}
	}}
}

/*
 * [mod P3] Gauss-Seidel rojo-negro por slabs.
 * Por iteración: medio barrido de celdas con (i+j+k) par, intercambio de halos, medio barrido de
 * impares, intercambio de halos, set_bnd.
 *  - Estrategia A (overlap=0): MPI_Sendrecv bloqueante tras cada medio barrido.
 *  - Estrategia B (overlap=1): se calculan primero los planos frontera (klo y khi), se lanzan los
 *    envíos/recepciones no bloqueantes y se calcula el interior (klo+1..khi-1) mientras viajan los mensajes.
 *    Es seguro: los planos interiores no leen los fantasma, y los planos enviados ya no se modifican.
 */
void lin_solve ( int M, int N, int O, int b, float * x, float * x0, float a, float c )
{
	int l, color;
	MPI_Request rq[4];

	halo_exchange ( x );   /* la estimación inicial de los planos fantasma debe ser la de los vecinos */
	for ( l=0 ; l<LINEARSOLVERTIMES ; l++ ) {
		for ( color=0 ; color<2 ; color++ ) {
			if ( sl.overlap && sl.nprocs > 1 ) {
				sweep ( M, N, x, x0, a, c, color, sl.klo, sl.klo );
				if ( sl.khi > sl.klo ) sweep ( M, N, x, x0, a, c, color, sl.khi, sl.khi );
				halo_start ( x, rq );
				if ( sl.khi - sl.klo >= 2 ) sweep ( M, N, x, x0, a, c, color, sl.klo+1, sl.khi-1 );
				halo_finish ( rq );
			} else {
				sweep ( M, N, x, x0, a, c, color, sl.klo, sl.khi );
				halo_exchange ( x );
			}
		}
		set_bnd ( M, N, O, b, x );
	}
}

void diffuse (  int M, int N, int O, int b, float * x, float * x0, float diff, float dt )
{
	int max = MAX(MAX(M, N), MAX(N, O));
	float a=dt*diff*max*max*max;
	lin_solve ( M, N, O, b, x, x0, a, 1+6*a );
}

/* [mod P3] requiere d0 replicado en todos los procesos (allgather previo); calcula solo mis planos */
void advect (  int M, int N, int O, int b, float * d, float * d0, float * u, float * v, float * w, float dt )
{
	int i, j, k, i0, j0, k0, i1, j1, k1;
	float x, y, z, s0, t0, s1, t1, u1, u0, dtx,dty,dtz;

	dtx=dty=dtz=dt*MAX(MAX(M, N), MAX(N, O));

	for ( k=sl.klo ; k<=sl.khi ; k++ ) { for ( j=1 ; j<=N ; j++ ) { for ( i=1 ; i<=M ; i++ ) {
		x = i-dtx*u[IX(i,j,k)]; y = j-dty*v[IX(i,j,k)]; z = k-dtz*w[IX(i,j,k)];
		if (x<0.5f) x=0.5f; if (x>M+0.5f) x=M+0.5f; i0=(int)x; i1=i0+1;
		if (y<0.5f) y=0.5f; if (y>N+0.5f) y=N+0.5f; j0=(int)y; j1=j0+1;
		if (z<0.5f) z=0.5f; if (z>O+0.5f) z=O+0.5f; k0=(int)z; k1=k0+1;

		s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1; u1 = z-k0; u0 = 1-u1;
		d[IX(i,j,k)] = s0*(t0*u0*d0[IX(i0,j0,k0)]+t1*u0*d0[IX(i0,j1,k0)]+t0*u1*d0[IX(i0,j0,k1)]+t1*u1*d0[IX(i0,j1,k1)])+
			s1*(t0*u0*d0[IX(i1,j0,k0)]+t1*u0*d0[IX(i1,j1,k0)]+t0*u1*d0[IX(i1,j0,k1)]+t1*u1*d0[IX(i1,j1,k1)]);
	}}}

	set_bnd (M, N, O, b, d );
}

/* [mod P3] requiere los planos fantasma de w actualizados (div usa w en k±1) */
void project ( int M, int N, int O, float * u, float * v, float * w, float * p, float * div )
{
	int i, j, k;

	for ( k=sl.klo ; k<=sl.khi ; k++ ) { for ( j=1 ; j<=N ; j++ ) { for ( i=1 ; i<=M ; i++ ) {
		div[IX(i,j,k)] = -1.0/3.0*((u[IX(i+1,j,k)]-u[IX(i-1,j,k)])/M+(v[IX(i,j+1,k)]-v[IX(i,j-1,k)])/M+(w[IX(i,j,k+1)]-w[IX(i,j,k-1)])/M);
		p[IX(i,j,k)] = 0;
	}}}

	set_bnd ( M, N, O, 0, div ); set_bnd (M, N, O, 0, p );

	lin_solve ( M, N, O, 0, p, div, 1, 6 );     /* termina con los fantasma de p actualizados */

	for ( k=sl.klo ; k<=sl.khi ; k++ ) { for ( j=1 ; j<=N ; j++ ) { for ( i=1 ; i<=M ; i++ ) {
		u[IX(i,j,k)] -= 0.5f*M*(p[IX(i+1,j,k)]-p[IX(i-1,j,k)]);
		v[IX(i,j,k)] -= 0.5f*M*(p[IX(i,j+1,k)]-p[IX(i,j-1,k)]);
		w[IX(i,j,k)] -= 0.5f*M*(p[IX(i,j,k+1)]-p[IX(i,j,k-1)]);
	}}}

	set_bnd (  M, N, O, 1, u ); set_bnd (  M, N, O, 2, v );set_bnd (  M, N, O, 3, w);
}

void dens_step ( int M, int N, int O, float * x, float * x0, float * u, float * v, float * w, float diff, float dt )
{
	add_source ( M, N, O, x, x0, dt );
	SWAP ( x0, x ); diffuse ( M, N, O, 0, x, x0, diff, dt );
	SWAP ( x0, x );
	allgather_field ( x0 );                                  /* [mod P3] advect lee x0 en cualquier posición */
	advect ( M, N, O, 0, x, x0, u, v, w, dt );
}

void vel_step ( int M, int N, int O, float * u, float * v,  float * w, float * u0, float * v0, float * w0, float visc, float dt )
{
	add_source ( M, N, O, u, u0, dt ); add_source ( M, N, O, v, v0, dt );add_source ( M, N, O, w, w0, dt );
	SWAP ( u0, u ); diffuse ( M, N, O, 1, u, u0, visc, dt );
	SWAP ( v0, v ); diffuse ( M, N, O, 2, v, v0, visc, dt );
	SWAP ( w0, w ); diffuse ( M, N, O, 3, w, w0, visc, dt );
	project ( M, N, O, u, v, w, u0, v0 );
	SWAP ( u0, u ); SWAP ( v0, v );SWAP ( w0, w );
	allgather_field ( u0 ); allgather_field ( v0 ); allgather_field ( w0 );   /* [mod P3] campo de velocidades replicado */
	advect ( M, N, O, 1, u, u0, u0, v0, w0, dt ); advect ( M, N, O, 2, v, v0, u0, v0, w0, dt );advect ( M, N, O, 3, w, w0, u0, v0, w0, dt );
	halo_exchange ( w );                                                       /* [mod P3] project necesita w en k±1 */
	project ( M, N, O, u, v, w, u0, v0 );
}
