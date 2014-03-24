/*
Magnetic pendulum model solving nonlinear second order differential equations
Code released under MIT LICENCE - more @ http://anze.mit-license.org/
Author: Alberto Anzellotti
co-Author: Giovanni Pederiva
citation would be appreciated :) thanks and good coding!
*/



#include <iostream>
#include <math.h>
#include <chrono>
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
const double d_time = .05;	//delta time  (time interval)

struct attr{
int x;
int y;
int color[3];
};

attr a,b,c;

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
    s_d[0] = -R*s_d[0] - C*s_d[1] + Q* ((a.x-s_d[1])/pow(pow(a.x-s_d[1],2)+pow(a.y-s_d[3],2)+1e-12,3/2)+ (b.x-s_d[1])/pow(pow(b.x-s_d[1],2)+pow(b.y-s_d[3],2)+1e-12,3/2)+ (c.x-s_d[1])/pow(pow(c.x-s_d[1],2)+pow(c.y-s_d[3],2)+1e-12,3/2));
    s_d[1] = s_d[0];
    s_d[2] = -R*s_d[2] - C*s_d[3] + Q*( (a.y-s_d[3])/pow(pow(a.x-s_d[1],2)+pow(a.y-s_d[3],2)+1e-12,3/2)+ (b.y-s_d[3])/pow(pow(b.x-s_d[1],2)+pow(b.y-s_d[3],2)+1e-12,3/2)+ (c.y-s_d[3])/pow(pow(c.x-s_d[1],2)+pow(c.y-s_d[3],2)+1e-12,3/2) );
    s_d[3] = s_d[2];
}


void write_observer( const state_type &s_i , const double t )
{
    cout << t << '\t' << "| \t" << s_i[0] << '\t' << s_i[1] << '\t' << s_i[2] << '\t' << s_i[3] << endl;
}


int main(int argc, char **argv)
{
    setup();
    double x_=10;
    double y_=-10;
    drop(x_,y_);
    INIT_TIMER;
    START_TIMER
    integrate_const(euler< state_type >(), mag_pend , s , s_time , e_time , d_time , write_observer ); // note 1
    STOP_TIMER
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


*/

