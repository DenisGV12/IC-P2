/*
 * main.c — versión sin interfaz (headless) para la Práctica 2 de IC (UA, 2026-27)
 *
 * Basado en BlainMaguire/3dfluid (MIT, commit 1e58a59), a su vez basado en
 * Jos Stam, "Real-Time Fluid Dynamics for Games", GDC 2003.
 *
 * Cambios respecto al main.c original:
 *   - Eliminado todo GLUT/OpenGL (ventana, render, teclado, ratón).
 *   - Parámetros por línea de comandos: N (malla N^3), pasos, dt.
 *   - Las pulsaciones de teclas (X, W, D, S) se sustituyen por un guion
 *     determinista de fuentes y fuerzas (ver scenario_events()).
 *   - Medición con clock_gettime(CLOCK_MONOTONIC): tiempo total y por fase.
 *   - Checksums finales para comprobar que todas las compilaciones dan el mismo resultado.
 *
 * Uso: ./fluid [N] [pasos] [dt]
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define IX(i,j,k) ((i)+(M+2)*(j) + (M+2)*(N+2)*(k))

extern void dens_step ( int M, int N, int O, float * x, float * x0, float * u, float * v, float * w, float diff, float dt );
extern void vel_step  ( int M, int N, int O, float * u, float * v, float * w, float * u0, float * v0, float * w0, float visc, float dt );

static int M, N, O;
static float dt     = 0.4f;   /* igual que el original */
static float diff   = 0.0f;
static float visc   = 0.0f;
static float force  = 10.0f;  /* igual que el original */
static float source = 200.0f;

static float *u, *v, *w, *u_prev, *v_prev, *w_prev, *dens, *dens_prev;

static double now ( void )
{
	struct timespec t;
	clock_gettime ( CLOCK_MONOTONIC, &t );
	return t.tv_sec + t.tv_nsec * 1e-9;
}

static int allocate_data ( void )
{
	size_t size = (size_t)(M+2)*(N+2)*(O+2);
	u = calloc ( size, sizeof(float) );  v = calloc ( size, sizeof(float) );
	w = calloc ( size, sizeof(float) );  u_prev = calloc ( size, sizeof(float) );
	v_prev = calloc ( size, sizeof(float) ); w_prev = calloc ( size, sizeof(float) );
	dens = calloc ( size, sizeof(float) ); dens_prev = calloc ( size, sizeof(float) );
	if ( !u || !v || !w || !u_prev || !v_prev || !w_prev || !dens || !dens_prev ) {
		fprintf ( stderr, "cannot allocate data\n" );
		return 0;
	}
	return 1;
}

static void free_data ( void )
{
	free(u); free(v); free(w); free(u_prev); free(v_prev); free(w_prev);
	free(dens); free(dens_prev);
}

/*
 * Sustituye a get_force_source() y a las teclas del original. Guion fijo y reproducible:
 *   - cada 5 pasos:  fuente de densidad en el centro   (tecla X)
 *   - cada 15 pasos: fuerza alternando ejes x, y, z    (teclas D, W, S)
 */
static void scenario_events ( int M, int N, int O, int step, float * d, float * u, float * v, float * w )
{
	size_t i, size = (size_t)(M+2)*(N+2)*(O+2);
	for ( i=0 ; i<size ; i++ ) u[i] = v[i] = w[i] = d[i] = 0.0f;

	if ( step % 5 == 0 )
		d[IX(M/2, N/2, O/2)] = source;

	if ( step % 15 == 0 ) {
		switch ( (step/15) % 3 ) {
		case 0: u[IX(2,   N/2, O/2)] = force*10; break;
		case 1: v[IX(M/2, 2,   O/2)] = force*10; break;
		case 2: w[IX(M/2, N/2, 2  )] = force*10; break;
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
	int size  = argc > 1 ? atoi ( argv[1] ) : 64;
	int steps = argc > 2 ? atoi ( argv[2] ) : 100;
	if ( argc > 3 ) dt = (float)atof ( argv[3] );
	if ( size < 4 || steps < 1 ) {
		fprintf ( stderr, "Uso: %s [N>=4] [pasos>=1] [dt]\n", argv[0] );
		return 1;
	}
	M = N = O = size;
	if ( !allocate_data () ) return 1;

	double t_src = 0, t_vel = 0, t_dens = 0, t0, t1;
	double t_start = now ();

	for ( int s = 0 ; s < steps ; s++ ) {
		t0 = now ();
		scenario_events ( M, N, O, s, dens_prev, u_prev, v_prev, w_prev );
		t1 = now (); t_src += t1 - t0; t0 = t1;
		vel_step ( M, N, O, u, v, w, u_prev, v_prev, w_prev, visc, dt );
		t1 = now (); t_vel += t1 - t0; t0 = t1;
		dens_step ( M, N, O, dens, dens_prev, u, v, w, diff, dt );
		t1 = now (); t_dens += t1 - t0;
	}

	double total = now () - t_start;

	printf ( "N=%d pasos=%d dt=%.3f celdas=%lld mem=%.1f MB\n", size, steps, dt,
	         (long long)(M+2)*(N+2)*(O+2),
	         8.0 * (M+2)*(N+2)*(O+2) * sizeof(float) / (1024.0*1024.0) );
	printf ( "tiempo_total=%.4f s  (fuentes=%.4f  vel_step=%.4f  dens_step=%.4f)\n",
	         total, t_src, t_vel, t_dens );
	printf ( "checksum dens=%.6e  vel=%.6e\n",
	         checksum ( dens ), checksum ( u ) + checksum ( v ) + checksum ( w ) );
	/* línea CSV para los scripts de medición */
	printf ( "CSV,%d,%d,%.6f,%.6f,%.6f,%.6f\n", size, steps, total, t_src, t_vel, t_dens );

	free_data ();
	return 0;
}
