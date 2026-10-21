/*
 * slab.h — descomposición en planos (slabs) a lo largo de k para la Práctica 3 (MPI)
 *
 * La malla es (M+2)x(N+2)x(O+2) con IX(i,j,k) = i + (M+2)*j + (M+2)*(N+2)*k, de modo que
 * cada plano k es un bloque contiguo de memoria. Cada proceso es responsable de un rango de
 * planos interiores [klo..khi] (el rango 0 también del plano de contorno 0 y el último
 * del plano O+1). Cada proceso reserva la malla completa (mismos índices globales que el
 * secuencial), pero solo calcula sus planos y recibe por mensajes lo que necesita del resto.
 */
#ifndef SLAB_H
#define SLAB_H

#include <stddef.h>
#include <mpi.h>

typedef struct {
	int rank, nprocs;
	int klo, khi;          /* planos interiores de cálculo propios: [klo..khi] dentro de [1..O] */
	int own_lo, own_hi;    /* planos propios incluyendo contorno: rank 0 -> desde 0, último -> hasta O+1 */
	size_t plane;          /* floats por plano: (M+2)*(N+2) */
	int overlap;           /* 0: halos bloqueantes (Sendrecv); 1: no bloqueantes solapados con el cálculo interior */
	int *counts, *displs;  /* para MPI_Allgatherv (en floats) */
	/* contadores de comunicación/sincronización */
	double t_halo, t_gather;
	long   n_halo, n_gather;
	double bytes_halo, bytes_gather;
} Slab;

extern Slab sl;

int  slab_init ( int M, int N, int O );   /* 0 si error */
void slab_free ( void );
void halo_exchange ( float * x );         /* actualiza planos fantasma klo-1 y khi+1 con los de los vecinos */
void halo_start ( float * x, MPI_Request rq[4] );   /* versión no bloqueante: inicia el intercambio */
void halo_finish ( MPI_Request rq[4] );              /* ...y espera a que termine */
void allgather_field ( float * x );       /* replica el campo completo en todos los procesos */

#endif
