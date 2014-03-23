#include "lbe.hpp"
#define velocity u
#define density rho
#include "..\display.hpp"

int main(int argc, char *argv[]){

	initialize();

	//set_rectangle_obstacle ( WIDTH / 5,  HEIGHT/2-HEIGHT/8, 1, HEIGHT/4 );
	//set_circle_obstacle ( WIDTH / 5,  HEIGHT/2, HEIGHT/20 );

	/*save_particles();
	save_velocity();
	save_velocity_x();
	save_velocity_y();*/

	
	gl_main ( argc, argv );

	return 0;
}