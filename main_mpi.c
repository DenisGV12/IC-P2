/*
 * main_mpi.c — Práctica 3 (IC, UA 2026-27): solver de fluidos 3D con MPI
 *
 * Mismo escenario determinista que main.c (práctica 2). Cada proceso reserva la malla completa
 * (índices globales idénticos al secuencial) pero solo calcula sus planos k (ver slab.h).
 *
 * Uso: mpirun -np P ./fluid_mpi [N] [pasos] [dt] [modo]     modo 0: halos bloqueantes (por defecto), 1: solapados
 * Salida (solo rank 0): tiempos, checksums (idénticos para cualquier P) y una línea CSV:
 *   CSV,N,pasos,P,modo,total,vel_step,dens_step,comm_max,comm_media,halo_max,gather_max,n_halo,n_gather,MB_halo,MB_gather
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "slab.h"

#define IX(i,j,k) ((i)+(M+2)*(j) + (M+2)*(N+2)*(k))

extern void dens_step ( int M, int N, int O, float * x, float * x0, float * u, float * v, float * w, float diff, float dt );
extern void vel_step  ( int M, int N, int O, float * u, float * v, float * w, float * u0, float * v0, float * w0, float visc, float dt );

static int M, N, O;
static float dt     = 0.4f;
static float diff   = 0.0f;
static float visc   = 0.0f;
static float force  = 10.0f;
static float source = 200.0f;

static float *u, *v, *w, *u_prev, *v_prev, *w_prev, *dens, *dens_prev;

static int allocate_data ( void )
{
	size_t size = (size_t)(M+2)*(N+2)*(O+2);
	u = calloc ( size, sizeof(float) );  v = calloc ( size, sizeof(float) );
	w = calloc ( size, sizeof(float) );  u_prev = calloc ( size, sizeof(float) );
	v_prev = calloc ( size, sizeof(float) ); w_prev = calloc ( size, sizeof(float) );
	dens = calloc ( size, sizeof(float) ); dens_prev = calloc ( size, sizeof(float) );
	return u && v && w && u_prev && v_prev && w_prev && dens && dens_prev;
}

static void free_data ( void )
{
	free(u); free(v); free(w); free(u_prev); free(v_prev); free(w_prev);
	free(dens); free(dens_prev);
}

/* Igual que scenario_events() del secuencial, pero cada proceso solo pone a cero y rellena sus planos. */
static void scenario_events ( int step, float * d, float * uu, float * vv, float * ww )
{
	size_t i, lo = (size_t)sl.own_lo * sl.plane, hi = ((size_t)sl.own_hi + 1) * sl.plane;
	for ( i=lo ; i<hi ; i++ ) uu[i] = vv[i] = ww[i] = d[i] = 0.0f;

	#define MINE(k) ((k) >= sl.own_lo && (k) <= sl.own_hi)
	if ( step % 5 == 0 && MINE(O/2) )
		d[IX(M/2, N/2, O/2)] = source;

	if ( step % 15 == 0 ) {
		switch ( (step/15) % 3 ) {
		case 0: if ( MINE(O/2) ) uu[IX(2,   N/2, O/2)] = force*10; break;
		case 1: if ( MINE(O/2) ) vv[IX(M/2, 2,   O/2)] = force*10; break;
		case 2: if ( MINE(2)   ) ww[IX(M/2, N/2, 2  )] = force*10; break;
		}
	}
}

static double checksum ( const float * x )
{
	size_t i, size = (size_t)(M+2)*(N+2)*(O+2);
	double s = 0.0;
	for ( i=0 ; i<size ; i++ ) s += fabs ( (double)x[i] );
	return s;
}

int main ( int argc, char ** argv )
{
	MPI_Init ( &argc, &argv );
	int rank, np;
	MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
	MPI_Comm_size ( MPI_COMM_WORLD, &np );

	int size  = argc > 1 ? atoi ( argv[1] ) : 64;
	int steps = argc > 2 ? atoi ( argv[2] ) : 100;
	if ( argc > 3 ) dt = (float)atof ( argv[3] );
	int modo = argc > 4 ? atoi ( argv[4] ) : 0;
	if ( size < 4 || steps < 1 ) {
		if ( rank == 0 ) fprintf ( stderr, "Uso: mpirun -np P %s [N>=4] [pasos>=1] [dt] [modo 0|1]\n", argv[0] );
		MPI_Finalize (); return 1;
	}
	M = N = O = size;
	if ( !slab_init ( M, N, O ) ) {
		if ( rank == 0 ) fprintf ( stderr, "Error: hacen falta al menos N=%d planos para %d procesos\n", O, np );
		MPI_Finalize (); return 1;
	}
	sl.overlap = modo;
	if ( !allocate_data () ) { fprintf ( stderr, "[rank %d] cannot allocate data\n", rank ); MPI_Abort ( MPI_COMM_WORLD, 1 ); }

	double t_vel = 0, t_dens = 0, t0, t1;
	MPI_Barrier ( MPI_COMM_WORLD );
	double t_start = MPI_Wtime ();

	for ( int s = 0 ; s < steps ; s++ ) {
		scenario_events ( s, dens_prev, u_prev, v_prev, w_prev );
		t0 = MPI_Wtime ();
		vel_step ( M, N, O, u, v, w, u_prev, v_prev, w_prev, visc, dt );
		t1 = MPI_Wtime (); t_vel += t1 - t0; t0 = t1;
		dens_step ( M, N, O, dens, dens_prev, u, v, w, diff, dt );
		t1 = MPI_Wtime (); t_dens += t1 - t0;
	}

	MPI_Barrier ( MPI_COMM_WORLD );
	double local_total = MPI_Wtime () - t_start;   /* el tiempo paralelo es el del proceso más lento */

	/* --- estadísticas (fuera de la región medida) --- */
	double v_loc[6] = { local_total, t_vel, t_dens, sl.t_halo + sl.t_gather, sl.t_halo, sl.t_gather };
	double v_max[6], v_sum[6];
	MPI_Reduce ( v_loc, v_max, 6, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
	MPI_Reduce ( v_loc, v_sum, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	double mb_halo = sl.bytes_halo / 1048576.0, mb_gather = sl.bytes_gather / 1048576.0, mb_tot[2];
	double mb_loc[2] = { mb_halo, mb_gather };
	MPI_Reduce ( mb_loc, mb_tot, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

	/* checksums: se reúne el resultado completo para sumar en el mismo orden que el secuencial */
	allgather_field ( dens ); allgather_field ( u ); allgather_field ( v ); allgather_field ( w );

	if ( rank == 0 ) {
		printf ( "N=%d pasos=%d dt=%.3f procesos=%d modo=%d planos/proceso=%d..%d\n", size, steps, dt, np, modo,
		         O/np, O/np + (O%np ? 1 : 0) );
		printf ( "tiempo_total=%.4f s  (vel_step=%.4f  dens_step=%.4f)\n", v_max[0], v_max[1], v_max[2] );
		printf ( "comunicacion+sincronizacion: max=%.4f s  media=%.4f s  (%.1f%% del total)  halo=%.4f  allgather=%.4f\n",
		         v_max[3], v_sum[3]/np, 100.0*(v_sum[3]/np)/v_max[0], v_max[4], v_max[5] );
		printf ( "checksum dens=%.6e  vel=%.6e\n", checksum ( dens ), checksum ( u ) + checksum ( v ) + checksum ( w ) );
		printf ( "CSV,%d,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%ld,%ld,%.3f,%.3f\n", size, steps, np, modo,
		         v_max[0], v_max[1], v_max[2], v_max[3], v_sum[3]/np, v_max[4], v_max[5],
		         sl.n_halo, sl.n_gather, mb_tot[0], mb_tot[1] );
	}

	free_data ();
	slab_free ();
	MPI_Finalize ();
	return 0;
}
