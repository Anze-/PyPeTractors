#include <iostream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <math.h>
#include <stdio.h>
#include <chrono>

#define CHUNKSIZE 1

using namespace std;
using namespace boost::numeric::odeint;


#define R @R@
#define C @C@
#define Q @Q@
#define r @r@

#define img_d @img_d@
#define img_s (int)(img_d*2*D)
#define D @D@



#define N_ITER 35    //number of iterations for asintotic point

//time static vars
const double s_time = 0.0;	//start time
const double e_time = 1000.0;	//end time
const double d_time = 0.05;	//delta time (time interval)



typedef vector< double > state_type;

struct attr
{
	double x;
	double y;
	int color[3];
};

attr a, b, c;


void p_drop()
{
    a.color[0] = 1;     //red
    a.color[1] = 0; 
    a.color[2] = 0; 
    b.color[0] = 0;       //green
    b.color[1] = 1; 
    b.color[2] = 0; 
    c.color[0] = 0;       //blue
    c.color[1] = 0; 
    c.color[2] = 1; 
    
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

void asint_checker(const state_type &x , const double t ){}

void doTheMath(int row)
{
	p_drop();


	const double dt = 0.1;

	int i_col, i_iter,rep;

	char flag;
	char row_string[10];
	sprintf(row_string, "%d.pp", row);
	FILE * output;
	output = fopen(row_string, "w+");

	char mat[1000000];

	for(i_col = 0; i_col < img_s / 2; i_col++)
	{
		int tol = 5; //tolerance
		state_type x(4);
		x[0] = 0.0;
		x[1] = -D+((double)(i_col)/img_d);
		x[2] = 0.0;
		x[3] = -D+((double)(row)/img_d);

		if(x[1] < a.x + tol && x[1] > a.x - tol && x[3] < a.y + tol&& x[3] > a.y - tol)
		{
			mat[i_col] = 'r';
		 	mat[img_s - i_col - 1] = 'g';
		}
		else if(x[1] < b.x + tol && x[1] > b.x - tol && x[3] < b.y + tol && x[3] > b.y - tol)
		{
		  	mat[i_col] = 'g';
		  	mat[img_s - i_col - 1] = 'r';
		}
		else if(x[1] < c.x + tol && x[1] > c.x - tol && x[3] < c.y + tol && x[3] > c.y - tol)
		{
			mat[i_col] = 'b';
		  	mat[img_s - i_col - 1] = 'b';
		}
		else
		{
			flag = 'a';
			rep = 0;
			for(i_iter = 0; i_iter <= N_ITER; i_iter++)
			{
				integrate_n_steps(controlled_runge_kutta< runge_kutta_dopri5< state_type > >(), tractors , x , 0.0 , dt , 100 , asint_checker);
				if(x[1] < a.x + tol && x[1] > a.x - tol && x[3] < a.y + tol&& x[3] > a.y - tol)
				{
					if(flag == 'r' && rep == 4 || flag == 'r' && i_iter == N_ITER )
					{
						mat[i_col] = 'r';
						mat[img_s - i_col - 1] = 'g';
						break;
					}
					else if(flag == 'r') rep += 1;
					else 
					{
						flag ='r';
						rep = 0;
					}
				}
				else if(x[1] < b.x + tol && x[1] > b.x - tol && x[3] < b.y + tol && x[3] > b.y - tol)
				{
					if(flag == 'g' && rep == 4 || flag == 'r' && i_iter == N_ITER )
					{
						mat[i_col] = 'g';
						mat[img_s - i_col - 1] = 'r';
						break;
					}
					else if(flag == 'g') rep += 1;
					else 
					{
						flag ='g';
						rep = 0;
					}
				}
				else if(x[1] < c.x + tol && x[1] > c.x - tol && x[3] < c.y + tol && x[3] > c.y - tol)
				{
					if(flag == 'b' && rep == 4  || flag == 'r' && i_iter == N_ITER )
					{
						mat[i_col] = 'b';
						mat[img_s - i_col - 1] = 'b';
						break;
					}
					else if(flag == 'b') rep += 1;
					else 
					{
						flag ='b';
						rep = 0;
					}
				}
				else if (i_iter == N_ITER)
				{
					mat[i_col] = 'w';
					mat[img_s - i_col - 1] = 'w';
				}
			}
		}						   
	}

	for(i_col = 0; i_col < img_s; i_col++)
	{
		if(mat[i_col] == 'r')
			fprintf(output, "%d %d %d  ", a.color[0], a.color[1], a.color[2]);
		else if(mat[i_col] == 'g')
			fprintf(output, "%d %d %d  ", b.color[0], b.color[1], b.color[2]);
		else if(mat[i_col] == 'b')
			fprintf(output, "%d %d %d  ", c.color[0], c.color[1], c.color[2]);
		else if(mat[i_col] == 'w')
			fprintf(output, "%d %d %d  ", 1,1,1);
	}
	fprintf(output, "\n");
	fclose(output);
}
