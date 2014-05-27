/*
Magnetic pendulum model solving nonlinear second order differential equations
Code released under MIT LICENCE - more @ http://anze.mit-license.org/
Authors: Alberto Anzellotti Giovanni Pederiva
citation would be appreciated :) thanks and good coding!
*/

#include <iostream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <math.h>
#include <stdio.h>
#include <chrono>

using namespace std;
using namespace boost::numeric::odeint;

#define TIMING
#ifdef TIMING
#define INIT_TIMER auto start = std::chrono::high_resolution_clock::now();
#define START_TIMER start = std::chrono::high_resolution_clock::now();
#define STOP_TIMER std::cout << "It took " << std::chrono::duration_cast<std::chrono::seconds>( std::chrono::high_resolution_clock::now()-start).count() << " s " << std::endl;
#else
#define INIT_TIMER
#define START_TIMER
#define STOP_TIMER
#endif


const double R = 0.01;
const double C = 0.01;
const double Q = 30;
const double D = 40;
const double r = 30;


//time static vars
const double s_time = 0.0;	//start time
const double e_time = 1000.0;	//end time
const double d_time = 0.05;	//delta time (time interval)


//image static vars
const int img_q = 240;
const int radius = 2;

typedef vector< double > state_type;
typedef double time_type;

struct attr
{
    double x;
    double y;
    int color[3];
};

attr a, b, c;


void p_drop()
{
    a.color[0] = 255;     //red
    a.color[1] = 0; 
    a.color[2] = 0; 
    b.color[0] = 0;       //green
    b.color[1] = 255; 
    b.color[2] = 0; 
    c.color[0] = 0;       //blue
    c.color[1] = 0; 
    c.color[2] = 255; 
    
    c.x = 0;                                                                       //NORTH
    c.y = r;
    b.x = c.x * (-0.49999999999999978) - ( c.y * 0.86602540378443871);             //SOUTH EAST
    b.y = c.x * 0.86602540378443871 + ( c.y * (-0.49999999999999978));
    a.x = c.x *(-0.50000000000000044) - ( c.y * (-0.86602540378443837));           //SOUTH WEST
    a.y = c.x * (-0.86602540378443837) + ( c.y * (-0.50000000000000044));
}



void tractors( const state_type &x , state_type &dxdt , double  )
{
    dxdt[0] = -R*x[0] - C*x[1] + Q * ((a.x-x[1])/pow(pow(a.x-x[1],2)+pow(a.y-x[3],2)+0.00000000000001,3/2)+ (b.x-x[1])/pow(pow(b.x-x[1],2)+pow(b.y-x[3],2)+1e-12,3/2)+ (c.x-x[1])/pow(pow(c.x-x[1],2)+pow(c.y-x[3],2)+1e-12,3/2));
    dxdt[1] = x[0];
    dxdt[2] = -R*x[2] - C*x[3] + Q*( (a.y-x[3])/pow(pow(a.x-x[1],2)+pow(a.y-x[3],2)+0.00000000000001,3/2)+ (b.y-x[3])/pow(pow(b.x-x[1],2)+pow(b.y-x[3],2)+1e-12,3/2)+ (c.y-x[3])/pow(pow(c.x-x[1],2)+pow(c.y-x[3],2)+1e-12,3/2));
    dxdt[3] = x[2];
}

void mag_pend( const state_type &s_i, state_type &s_d , double  t  )
{
    s_d[0] = -R*s_i[0] - C*s_i[1] + Q* ((a.x-s_i[1])/pow(pow(a.x-s_i[1],2)+pow(a.y-s_i[3],2)+1e-12,3/2)+ (b.x-s_i[1])/pow(pow(b.x-s_i[1],2)+pow(b.y-s_i[3],2)+1e-12,3/2)+ (c.x-s_i[1])/pow(pow(c.x-s_i[1],2)+pow(c.y-s_i[3],2)+1e-12,3/2));
    s_d[1] = s_i[0];	
    s_d[2] = -R*s_i[2] - C*s_i[3] + Q*( (a.y-s_i[3])/pow(pow(a.x-s_i[1],2)+pow(a.y-s_i[3],2)+1e-12,3/2)+ (b.y-s_i[3])/pow(pow(b.x-s_i[1],2)+pow(b.y-s_i[3],2)+1e-12,3/2)+ (c.y-s_i[3])/pow(pow(c.x-s_i[1],2)+pow(c.y-s_i[3],2)+1e-12,3/2) );
    s_d[3] = s_i[2];
}


void write_tractors( const state_type &x , const double t )
{
    //cout << t << '\t' << '\t'<< x[0] << '\t' << '\t' << x[1] << '\t' << '\t' << x[2] << '\t' << '\t' << x[3] << endl;
}

int main(int argc, char **argv)
{
    p_drop();
    state_type x(4);
    
    const double dt = 0.1;
    
    int tol = 5; //tolerance
    
    FILE * output;
    output = fopen("output.ppm", "w+");
    fprintf(output, "P3\n%d %d \n255\n", img_q, img_q);
    for(int i_row = 0; i_row < img_q; i_row++)
    {
        INIT_TIMER;
        START_TIMER
        for(int i_col = 0; i_col < img_q; i_col++)
	{
	    if(sqrt(pow(radius*D*((-img_q/2)+i_col)/img_q,2)+pow(radius*D*((-img_q/2)+i_row)/img_q,2)) < D*radius)
	    {
	      x[0] = 0.0;
              x[1] = radius*D*((-img_q/2)+i_col)/img_q;
              x[2] = 0.0;
              x[3] = radius*D*((-img_q/2)+i_row)/img_q;
	      integrate_n_steps(controlled_runge_kutta< runge_kutta_dopri5< state_type > >(), tractors , x , 0.0 , dt , 2500 , write_tractors );
              //integrate( tractors , x , 0.0, 50.0 , dt , write_tractors );
	      //integrate_const(euler<state_type>(), mag_pend , x , 0.0, 200.0 , 0.05 , write_tractors);
	      
	      
	      
              if(x[1] < a.x + tol && x[1] > a.x - tol && x[3] < a.y + tol&& x[3] > a.y - tol)
                   fprintf(output, "%d %d %d  ", a.color[0], a.color[1], a.color[2]);
              else if(x[1] < b.x + tol && x[1] > b.x - tol && x[3] < b.y + tol && x[3] > b.y - tol)
                   fprintf(output, "%d %d %d  ", b.color[0], b.color[1], b.color[2]);
              else if(x[1] < c.x + tol && x[1] > c.x - tol && x[3] < c.y + tol && x[3] > c.y - tol)
                   fprintf(output, "%d %d %d  ", c.color[0], c.color[1], c.color[2]);
	      else
	      {
		   fprintf(output, "%d %d %d ", 255, 255, 255);
		   printf("%d %d \n", i_row, i_col);
	      }
	    }
	    else
	      fprintf(output, "%d %d %d ", 255, 255, 255);
        }
        printf("%d\n", i_row);
	fprintf(output, " \n");
	STOP_TIMER
    }
    fclose(output);
}
