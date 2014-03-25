/*
Magnetic pendulum model solving nonlinear second order differential equations
Code released under MIT LICENCE - more @ http://anze.mit-license.org/
Author: Alberto Anzellotti
co-Author: Giovanni Pederiva
citation would be appreciated :) thanks and good coding!
*/

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <chrono>
#include <array>
#include <cstring>
#include <algorithm>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

#define TIMING
 
#ifdef TIMING
#define INIT_TIMER auto start = std::chrono::high_resolution_clock::now();
#define START_TIMER  start = std::chrono::high_resolution_clock::now();
#define STOP_TIMER  std::cout << "It took " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now()-start).count() << " ms " << std::endl; 
#else
#define INIT_TIMER
#define START_TIMER
#define STOP_TIMER
#endif

int SEId(double array[], int size)
{
    int index = 0;

    for(int i = 1; i < size; i++)
    {
        if(array[i] < array[index])
            index = i;              
    }

    return index;
}



//pendulum_drop static vars
const double R=.01;		//friction
const double C=.01;		//gravity lin. approx. self restoring force
const double Q=30;		//magnetic approx. charge force
const double D=40;		//length of the pendulum
const double r=30;		//distance of the magnets from the centre. (!smaller than D)
double x_drop,y_drop; 		//empty x_drop and y_drop

//time static vars
const double s_time = 0;		//start time
const double e_time = 150;		//end time
const double d_time = 0.05;	//delta time  (time interval) !must return an integer when e_time/d_time is called

double out[(int)(e_time/d_time)][3];

struct attr{
double x;
double y;
double dist;
};

attr a,b,c;

char * a_color = "255   0   0   ";
char * b_color = "  0 255   0   ";
char * c_color = "  0   0 255   ";

void setup(){
c.x=0;
c.y=r;
b.x = c.x * (-0.49999999999999978) - ( c.y * 0.86602540378443871); //south-east point
b.y = c.x * 0.86602540378443871 + ( c.y * (-0.49999999999999978));
a.x = c.x *(-0.50000000000000044) - ( c.y * (-0.86602540378443837)); //south-west point
a.y = c.x * (-0.86602540378443837) + ( c.y * (-0.50000000000000044));
}

typedef boost::array< double , 4 > state_type;

state_type s;


void drop(double &x_drop, double &y_drop){
s[0] = 0;// x'(t_0)
s[1] = x_drop;// x (t_0)
s[2] = 0;// y'(t_0)
s[3] = y_drop;// y (t_0)
}

void mag_pend( const state_type &s_i, state_type &s_d , double /* t */ )
{
    s_d[0] = -R*s_i[0] - C*s_i[1] + Q* ((a.x-s_i[1])/pow(pow(a.x-s_i[1],2)+pow(a.y-s_i[3],2)+1e-12,1.5)+ (b.x-s_i[1])/pow(pow(b.x-s_i[1],2)+pow(b.y-s_i[3],2)+1e-12,1.5)+ (c.x-s_i[1])/pow(pow(c.x-s_i[1],2.0)+pow(c.y-s_i[3],2)+1e-12,1.5));
    s_d[1] = s_i[0];
    s_d[2] = -R*s_i[2] - C*s_i[3] + Q*( (a.y-s_i[3])/pow(pow(a.x-s_i[1],2)+pow(a.y-s_i[3],2)+1e-12,1.5)+ (b.y-s_i[3])/pow(pow(b.x-s_i[1],2)+pow(b.y-s_i[3],2)+1e-12,1.5)+ (c.y-s_i[3])/pow(pow(c.x-s_i[1],2)+pow(c.y-s_i[3],2)+1e-12,1.5) );
    s_d[3] = s_i[2];
}

//to display data while proceeding
void write_observer( const state_type &s_i , const double t )
{
    cout << t << '\t' << "| \t" << s_i[0] << '\t' << s_i[1] << '\t' << s_i[2] << '\t' << s_i[3] << endl;
}

//to save data to a n x 3 matrix
void save_observer( const state_type &s_i , const double t )
{
    out[(int)(t/d_time)][0]=t; out[(int)(t/d_time)][1]=s_i[1]; out[(int)(t/d_time)][2]=s_i[3];
}

char * asint(double (*data_matrix)[3]){
    double l_x=out[(int)(e_time/d_time)-1][1];
    double l_y=out[(int)(e_time/d_time)-1][2];
    //cout << a.x << " " << b.x << " " << c.x << endl;
    //cout << a.y << " " << b.y << " " << c.y << endl;
    //cout << l_x << " and " << l_y << endl;
    a.dist=pow(pow(l_x-a.x,2)+pow(l_y-a.y,2),0.5);
    b.dist=pow(pow(l_x-b.x,2)+pow(l_y-b.y,2),0.5);
    c.dist=pow(pow(l_x-c.x,2)+pow(l_y-c.y,2),0.5);
    char * colors[3];colors[0]=a_color;colors[1]=b_color;colors[2]=c_color;
    double dists[3]={a.dist, b.dist, c.dist};
    int s_index=SEId(dists, 3);
    //cout << a.dist << " " << b.dist << " " << c.dist << endl;
    if (dists[s_index]<D/3.0){
	return colors[s_index];
    }else{
	return "  0   0   0   ";
    }
}

int main(int argc, char **argv)
{
    setup();
    double x_=-15;
    double y_=-5;
    drop(x_,y_);
    INIT_TIMER;
    START_TIMER
    integrate_const(euler< state_type >(), mag_pend , s , s_time , e_time , d_time , save_observer ); // note 1
    STOP_TIMER
    int i=0;
    //printf("%d",(int)(e_time/d_time));
    /*
    while(i<(int)(e_time/d_time)){
	cout << out[i][0] << '\t' << "| \t" << out[i][1] << '\t' << out[i][2] << endl;
	i++;
    }
    */
    cout << asint(out) << endl;
    //printf("last X: %f",out[(int)(e_time/d_time)-1][1]);
    //printf("last Y: %f\n",out[(int)(e_time/d_time)-1][2]);


    //write out the picture

    cout << "press ENTER to map . . ." <<endl;
    cin.ignore(1);

    double img_d=3; //image density, 1 means 1 pixel for unit of the D (wich is the radius of the sphere of the pendulum)
    int img_s = (int)(img_d*2*D);

    FILE * map_file;
    map_file=fopen("map.ppm", "w+");
    fprintf(map_file, "P3\n%d %d\n255\n", img_s, img_s);
    fclose(map_file);
    map_file=fopen("map.ppm", "a");

    i = 0;
    while(i<img_s){
	int j=0;
	double y_pos=-D+(i/img_d);
	while(j<img_s){
		double x_pos=-D+(j/img_d);
		setup();
		drop(x_pos,y_pos);
		integrate_const(euler< state_type >(), mag_pend , s , s_time , e_time , d_time , save_observer );
		fprintf(map_file, "%s",asint(out));
		j++;
	}
	fprintf(map_file, "\n");
	cout << i << endl;
	i++;
    }
    fclose(map_file);

}

/*			################     NOTES     ###################

0)compile with -std=gnu++11 option!

1)euler integration, to use runge kutta substitute 'euler' with one of the following:

Runge-Kutta 4	runge_kutta4	The classical Runge Kutta scheme, good general scheme without error control
Cash-Karp	runge_kutta_cash_karp54	Good general scheme with error estimation
Dormand-Prince 5	runge_kutta_dopri5	Standard method with error control and dense output
Fehlberg 78	runge_kutta_fehlberg78	Good high order method with error estimation
Adams-Bashforth-Moulton	adams_bashforth_moulton	Multi-step method with high performance
Controlled Error Stepper	controlled_runge_kutta	Error control for the Runge-Kutta steppers
Dense Output Stepper	dense_output_runge_kutta	Dense output for the Runge-Kutta steppers
Bulirsch-Stoer	bulirsch_stoer	Stepper with step size, order control and dense output. Very good if high precision is required.
Implicit Euler	implicit_euler	Basic implicit routine
Rosenbrock 4	rosenbrock4	Solver for stiff systems with error control and dense output
Symplectic Euler	symplectic_euler	Basic symplectic solver for separable Hamiltonian system
Symplectic RKN McLachlan	symplectic_rkn_sb3a_mclachlan	Symplectic solver for separable Hamiltonian system with order 6

2)ppm picture declaration (Wikipedia):
P3
# The P3 means colors are in ASCII, then 3 columns and 2 rows,
# then 255 for max color, then RGB triplets
3 2
255
255   0   0     0 255   0     0   0 255


*/
