#include "../point.hpp"
#include <utility>
using namespace std;

const int WIDTH = 50;
const int HEIGHT = 50;

const double TAU = 10.0;
const point U0 = point ( 0.0042, 0 );
const point UF = point ( 0.0065, 0 );
const double RHO0 = 1;
const double RHOF = 1;

const double CS2 = 1.0/3.0;

const int DIRECTIONS_NUMBER = 9;

double rho[HEIGHT][WIDTH];
point u[HEIGHT][WIDTH];
bool is_obstacle[HEIGHT][WIDTH], is_source[HEIGHT][WIDTH];

double f[2][HEIGHT][WIDTH][DIRECTIONS_NUMBER], (*f_current)[WIDTH][DIRECTIONS_NUMBER] = f[1], (*f_previous)[WIDTH][DIRECTIONS_NUMBER] = f[0], f_eq[HEIGHT][WIDTH][DIRECTIONS_NUMBER];
const int E[DIRECTIONS_NUMBER][2] = {
	{0,0},
	{1,0},{0,1},{-1,0},{0,-1}, 
	{1,1},{-1,1},{-1,-1},{1,-1}
};
point e[DIRECTIONS_NUMBER] = {
	point(0,0),
	point(1,0),point(0,1),point(-1,0),point(0,-1), 
	point(1,1),point(-1,1),point(-1,-1),point(1,-1)
};
enum{
	CENTER = 0,
	EAST, NORTH, WEST, SOUTH,
	NORTHEAST, NORTHWEST, SOUTHWEST, SOUTHEAST
};
const double W[DIRECTIONS_NUMBER] = {
	4.0/9.0,
	1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,
	1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,
};

void save_velocity( const char *filename ){
	static int save_time = 0;
	++save_time;
	FILE *fh = fopen ( filename, "a" );
	if ( fh == NULL ) {
		fprintf ( stderr, "Cannot open %s for append\n", filename );
		return;
	}
	fprintf ( fh, "%d\n", save_time );
	for (int y = 0; y < HEIGHT; y++){
		for (int x = 0; x < WIDTH; x++){
			fprintf ( fh, "(%+.2e;%+.2e) ", u[y][x].x, u[y][x].y );
		}
		fputc ( '\n', fh );
	}
	fclose ( fh );
}
void save_velocity_abs( const char *filename ){
	static int save_time = 0;
	++save_time;
	FILE *fh = fopen ( filename, "a" );
	if ( fh == NULL ) {
		fprintf ( stderr, "Cannot open %s for append\n", filename );
		return;
	}
	fprintf ( fh, "%d\n", save_time );
	for (int y = 0; y < HEIGHT; y++){
		for (int x = 0; x < WIDTH; x++){
			fprintf ( fh, "%e\t", u[y][x].length() );
		}
		fputc ( '\n', fh );
	}
	fclose ( fh );
}
void save_velocity_x( const char *filename ){
	static int save_time = 0;
	++save_time;
	FILE *fh = fopen ( filename, "a" );
	fprintf ( fh, "%d\n", save_time );
	for (int y = 0; y < HEIGHT; y++){
		for (int x = 0; x < WIDTH; x++){
			fprintf ( fh, "%+.2e ", u[y][x].x, u[y][x].y );
		}
		fputc ( '\n', fh );
	}
	fclose ( fh );
}
void save_density( const char *filename ){
	static int save_time = 0;
	++save_time;
	FILE *fh = fopen ( filename, "a" );
	fprintf ( fh, "%d\n", save_time );
	for (int y = 0; y < HEIGHT; y++){
		for (int x = 0; x < WIDTH; x++){
			fprintf ( fh, "%2e ", rho[y][x] );
		}
		fputc ( '\n', fh );
	}
	fclose ( fh );
}
void save_particles( const char *filename ){
	static int save_time = 0;
	++save_time;
	FILE *fh = fopen ( filename, "a" );
	fprintf ( fh, "%d\n", save_time );
	for (int y = 0; y < HEIGHT; y++){
		for (int x = 0; x < WIDTH; x++){
			fprintf ( fh, "(%+.2e", f_current[y][x][0] );
			for (int a = 1; a < DIRECTIONS_NUMBER; a++){
				fprintf ( fh, ";%+.2e", f_current[y][x][a] );
				if ( f_current[y][x][a] > 10 ) {
					continue;
				}
			}
			fputc ( ')', fh );
		}
		fputc ( '\n', fh );
	}
	fclose ( fh );
}

void reflect( double before[DIRECTIONS_NUMBER], double after[DIRECTIONS_NUMBER] ){
/*	swap ( particles[WEST], particles[EAST] );
	swap ( particles[NORTH], particles[SOUTH] );
	swap ( particles[NORTHWEST], particles[SOUTHEAST] );
	swap ( particles[NORTHEAST], particles[SOUTHWEST] );
*/
	after[CENTER] = before[CENTER];
	after[WEST] = before[EAST];
	after[EAST] = before[WEST];
	after[NORTH] = before[SOUTH];
	after[SOUTH] = before[NORTH];
	after[NORTHWEST] = before[SOUTHEAST];
	after[SOUTHEAST] = before[NORTHWEST];
	after[NORTHEAST] = before[SOUTHWEST];
	after[SOUTHWEST] = before[NORTHEAST];
}

void calculate_equilibrium(){
#pragma omp parallel for
	for (int y = 0; y < HEIGHT; y++){
		for (int x = 0; x < WIDTH; x++){
			for (int a = 0; a < DIRECTIONS_NUMBER; a++){
				/*f_eq[y][x][a] = rho[y][x] * W[a] * (
					1+
					e[a] * u[y][x] / CS2 +
					sqr ( e[a] * u[y][x] ) / ( 2.0 * sqr ( CS2 ) ) -
					sqr ( u[y][x] ) / ( 2.0 * CS2 ) + 
					cube ( e[a] * u[y][x] ) / ( 2.0 * cube ( CS2 ) ) - 
					e[a] * u[y][x] * sqr ( u[y][x] ) / ( 2.0 * sqr ( CS2 ) )
				);*/
				f_eq[y][x][a] = rho[y][x] * W[a] * (
					1+
					u[y][x] * e[a] / CS2 +
					(sqr(u[y][x]*e[a])-CS2*sqr(u[y][x]))/(2*sqr(CS2))
				);
#ifdef _DEBUG
				if ( f_eq[y][x][a] < 0.0 ) {
					fprintf ( stderr, "equilibrium density %e in (%d;%d)[%d] less than 0\n", f_eq[y][x][a], x, y, a );
				}
#endif
			}
		}
	}
}

void collide(){
#pragma omp parallel for
	for (int y = 0; y < HEIGHT; y++){
		for (int x = 0; x < WIDTH; x++){
			if ( is_obstacle[y][x] ) {
				reflect ( f_current[y][x], f_previous[y][x] );
				continue;
			}
			for (int a = 0; a < DIRECTIONS_NUMBER; a++){
				f_previous[y][x][a] = f_current[y][x][a] + ( f_eq[y][x][a] - f_current[y][x][a] ) / TAU;
			}
		}
	}
}

void stream(){
#pragma omp parallel for
	for (int y = 1; y < HEIGHT-1; y++){
		for (int x = 1; x < WIDTH-1; x++){
			if ( is_source[y][x] ) {
				continue;
			}
			for (int a = 0; a < DIRECTIONS_NUMBER; a++){
				f_current[y][x][a] = f_previous[y-E[a][1]][x-E[a][0]][a];
			}
		}
	}
}
void boundary(){
	//East && West
	for (int y = 0; y < HEIGHT; y++){
		for (int a = 0; a < DIRECTIONS_NUMBER; a++){
			f_current[y][0][a] = f_previous[(y-E[a][1]+HEIGHT)%HEIGHT][(0-E[a][0]+WIDTH)%WIDTH][a];
			f_current[y][WIDTH-1][a] = f_previous[(y-E[a][1]+HEIGHT)%HEIGHT][(WIDTH-1-E[a][0]+WIDTH)%WIDTH][a];
		}
	}

	//North && Soth
	for (int x = 0; x < WIDTH; x++){
		for (int a = 0; a < DIRECTIONS_NUMBER; a++){
			if ( 0-E[a][1] >= 0  && 0-E[a][1] < HEIGHT && x-E[a][0] >= 0 && x-E[a][0] < WIDTH ){
				f_current[0][x][a] = f_previous[0-E[a][1]][x-E[a][0]][a];
			} else {
				f_current[0][x][a] = 0.0;
			}
			if ( (HEIGHT-1)-E[a][1] >= 0  && (HEIGHT-1)-E[a][1] < HEIGHT && x-E[a][0] >= 0 && x-E[a][0] < WIDTH ){
				f_current[HEIGHT-1][x][a] = f_previous[HEIGHT-1-E[a][1]][x-E[a][0]][a];
			} else {
				f_current[HEIGHT-1][x][a] = 0.0;
			}
		}
	}
}
void pressure(){
	const double VISC = CS2 * (TAU - 0.5);
	
	for (int y = 0; y < HEIGHT; y++){
		for (int x = 0; x < WIDTH; x++){
			double fpois = rho[y][x] * 8.0 * VISC * UF.length() / sqr ( HEIGHT ) / 6; // Poiseuille force
			f_current[y][x][EAST] += fpois;
			f_current[y][x][NORTHEAST] += fpois;
			f_current[y][x][SOUTHEAST] += fpois;

			f_current[y][x][WEST] -= fpois;
			f_current[y][x][NORTHWEST] -= fpois;
			f_current[y][x][SOUTHWEST] -= fpois;
		}
	}
}
void calculate_macro(){
	for (int y = 0; y < HEIGHT; y++){
		for (int x = 0; x < WIDTH; x++){
			if ( is_source[y][x] ) {
				continue;
			}
			if ( is_obstacle[y][x] ) {
				//continue;
			}
			rho[y][x] = 0;
			u[y][x] = point ( 0.0, 0.0 );
			for (int a = 0; a < DIRECTIONS_NUMBER; a++){
				rho[y][x] += f_current[y][x][a];
				u[y][x] += f_current[y][x][a] * e[a];
			}
			if ( rho[y][x] > EPS ) {
				u[y][x] /= rho[y][x];
			} else {
				rho[y][x] = 0.0;
				u[y][x] = point ( 0.0, 0.0 );
			}
		}
	}
}

void step(){
	calculate_equilibrium();
	collide();
	stream();
	pressure();
	calculate_macro();
#ifdef _DEBUG
	static int step = 0;
	printf ( "%d\n", ++step );
	save_velocity("velocity.log");
	save_particles("lattice.log");
	save_density("density.log");
	save_velocity_x ( "velocity_x.log" );
	save_velocity_abs ( "velocity_abs.log" );
#endif
}
void set_circle_obstacle( const int x0, const int y0, const int radius ){
	for ( int y = y0-radius; y <= y0+radius; ++y ) {
		int max_x = int ( x0 + sqrt ( .0 + sqr ( radius ) - sqr ( y - y0 ) ) );
		for (int x = (int)ceil ( x0 - sqrt ( .0 + sqr ( radius ) - sqr ( y - y0 ) ) ); x < max_x; x++){
			is_obstacle[y][x] = true;
		}
	}
}
void set_rectangle_obstacle( const int x0, const int y0, const int width, const int height ){
	for (int y = 0; y < height; y++){
		for (int x = 0; x < width; x++){
			is_obstacle[y+y0][x+x0] = true;
		}
	}
}
void set_rectangle_source( const int x0, const int y0, const int width, const int height, double rho0, point u0 ){
	for (int y = 0; y < height; y++){
		for (int x = 0; x < width; x++){
			is_source[y+y0][x+x0] = true;
			rho[y+y0][x+x0] = rho0;
			u[y+y0][x+x0] = u0;
		}
	}
}

void initialize(){
	for (int y = 0; y < HEIGHT; y++){
		for (int x = 0; x < WIDTH; x++){
			u[y][x] = U0;
			rho[y][x] = RHO0;
		}
	}
	set_rectangle_obstacle ( 0, 0, WIDTH, 1 );
	set_rectangle_obstacle ( 0, HEIGHT-1, WIDTH, 1 );
	/*set_rectangle_source ( 0, 1, 1, HEIGHT-2, RHO0, U0 );
	set_rectangle_source ( WIDTH-1, 1, 1, HEIGHT-2, RHOF, UF );	*/
	calculate_equilibrium();
	memcpy ( &f_current[0][0][0], &f_eq[0][0][0], sizeof ( double ) * HEIGHT * WIDTH * DIRECTIONS_NUMBER );
#ifdef _DEBUG
	remove ( "velocity.log" );
	remove ( "lattice.log" );
	remove ( "density.log" );
	remove ( "velocity_x.log" );
#endif
}