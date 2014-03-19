#Magnetic pendulum model solving nonlinear second order differential equations
#Code released under MIT LICENCE - more @ http://anze.mit-license.org/
#Author: Alberto Anzellotti
#co-Author: Giovanni Pederiva
#citation would be appreciated :) thanks and good coding!

# initial conditions: x(t_0) y(t_0) t_0 attractors_x=[x1, x2, x3] attractors_y=[y1, y2, y3]

import numpy as np
from scipy.integrate import odeint
from matplotlib.pylab import *
from matplotlib import pylab
from math import *
from mpl_toolkits.mplot3d import Axes3D

 
def p_drop(x_drop, y_drop):
    #standard model constants
    global R, C, Q, D, r, s, a, b, c
    R=0.1		 #friction
    C=1		 #gravity lin. approx. self restoring force
    Q=10000000	 #magnetic approx. charge force
    D=25	 #length of the pendulum
    r=20	 #distance of the magnets from the centre. (!smaller than d)
    class a(object):
        pass
    class b(object):
        pass
    class c(object):
        pass
    c.x=0	 #north point
    c.y=r 
    b.x = c.x * 0.5 - ( c.y * 0.866025403784438) #south-east point
    b.y = c.x * 0.866025403784438 + ( c.y * 0.5)
    a.x = c.x * 0.5 - ( c.y * (-0.86602540378443)) #south-west point
    a.y = c.x * (-0.86602540378443) + ( c.y * 0.5)
    #initial conditions
    s = np.array([
    0,# x'(t_0)
    x_drop,# x(t_0)
    0,# y'(t_0)
    y_drop# y(t_0)
    ])


#fuuuuk il termine magnetico sommatoria blah e' costante e nulllo?
def d(s,t):
    return np.array([
    -R*s[0] - C*s[1] + Q* ((a.x-s[1])/np.power(np.power(a.x-s[1],2)+np.power(a.y-s[3],2)+D*D,3/2)+ (b.x-s[1])/np.power(np.power(b.x-s[1],2)+np.power(b.y-s[3],2)+D*D,3/2)+ (c.x-s[1])/np.power(np.power(c.x-s[1],2)+np.power(c.y-s[3],2)+D*D,3/2)),
    s[0],
    -R*s[2] - C*s[3] + Q*( (a.y-s[3])/np.power(np.power(a.x-s[1],2)+np.power(a.y-s[3],2)+D*D,3/2)+ (b.y-s[3])/np.power(np.power(b.x-s[1],2)+np.power(b.y-s[3],2)+D*D,3/2)+ (c.y-s[3])/np.power(np.power(c.x-s[1],2)+np.power(c.y-s[3],2)+D*D,3/2) ),
    s[2]]) 


def d1(sx,sy,t):
    return np.array([
    -R*s[0] - C*s[1] + (a.x-s[1])/np.power(np.power(a.x-s[1],2)+np.power(a.y-s[3],2)+D*D,3/2)+ (b.x-s[1])/np.power(np.power(b.x-s[1],2)+np.power(b.y-s[3],2)+D*D,3/2)+ (c.x-s[1])/np.power(np.power(c.x-s[1],2)+np.power(c.y-s[3],2)+D*D,3/2),
    s[0],
    -R*s[2] - C*s[3] + (a.y-s[3])/np.power(np.power(a.x-s[1],2)+np.power(a.y-s[3],2)+D*D,3/2)+ (b.y-s[3])/np.power(np.power(b.x-s[1],2)+np.power(b.y-s[3],2)+D*D,3/2)+ (c.y-s[3])/np.power(np.power(c.x-s[1],2)+np.power(c.y-s[3],2)+D*D,3/2) ,
    s[2]]) 


p_drop(16, -21)
d(s,0)

tmax = 100
dt = .0001
t = np.linspace(0, tmax, num=np.round(tmax/dt)+1)
'''
sol = odeint(d, s, t)
plt.plot(t, sol, label='lol')
plt.legend()
show()
'''
_,sol_x,_,sol_y = odeint(d, s, t).T
'''
plt.plot(sol_x, sol_y, label='X,Y dinamics')
plt.legend()

show()

plt.plot(t, sol_x, label='X,Y dinamics')
plt.legend()
show()

plt.plot(t, sol_y, label='X,Y dinamics')
plt.legend()
show()
'''
fig = plt.figure()
ax = Axes3D(fig)
# put 0s on the y-axis, and put the y axis on the z-axis
ax.plot(xs=sol_x, ys=sol_y, zs=t, zdir='z', label='xs=x, ys=y, zdir=time')
circle1=plt.Circle((a.x,a.y,0),.2,color='r')
fig.gca().add_artist(circle1)
xlim([-35,35])
fft_axes.set_autoscalex_on(False)
ylim([-35,35])

plt.show()

#    fig.savefig('coupled-ode-Python.png')
