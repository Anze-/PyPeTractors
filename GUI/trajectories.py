#Magnetic pendulum model solving nonlinear second order differential equations
#Code released under MIT LICENCE - more @ http://anze.mit-license.org/
#Authors: Alberto Anzellotti Giovanni Pederiva
#citation would be appreciated :) thanks and good coding!

import numpy as np
from scipy.integrate import odeint
from matplotlib.pylab import *
from math import *
from mpl_toolkits.mplot3d import Axes3D
from time import time
import sys
import getopt

#pendulum_drop static vars
XXX,YYY,R,C,Q,r,tmax,dt=sys.argv[1:]
XXX=float(XXX)
YYY=float(YYY)
R=float(R)
C=float(C)
Q=float(Q)
r=float(r)
tmax=float(tmax)
dt=float(dt)
#print Q
#sys.stdout.flush()
#print r
#sys.stdout.flush()
#print tmax
#sys.stdout.flush()
#print dt
#sys.stdout.flush()
class a(object):
    pass

class b(object):
    pass

class c(object):
    pass

a.name="a"
b.name="b"
c.name="c"
a.color="r"
b.color="g"
c.color="b"
c.x=0	 #north point
c.y=r
b.x = c.x * (-0.49999999999999978) - ( c.y * 0.86602540378443871) #south-east point
b.y = c.x * 0.86602540378443871 + ( c.y * (-0.49999999999999978))
a.x = c.x *(-0.50000000000000044) - ( c.y * (-0.86602540378443837)) #south-west point
a.y = c.x * (-0.86602540378443837) + ( c.y * (-0.50000000000000044))

def p_drop(x_drop, y_drop):
    #standard model constants
    global R, C, Q, D, r, s, a, b, c
    #initial conditions
    s = np.array([
    0,# x'(t_0)
    x_drop,# x(t_0)
    0,# y'(t_0)
    y_drop# y(t_0)
    ])


def d(s,t):
    return np.array([
    -R*s[0] - C*s[1] + Q* ((a.x-s[1])/np.power(np.power(a.x-s[1],2)+np.power(a.y-s[3],2)+1e-12,3/2)+ (b.x-s[1])/np.power(np.power(b.x-s[1],2)+np.power(b.y-s[3],2)+1e-12,3/2)+ (c.x-s[1])/np.power(np.power(c.x-s[1],2)+np.power(c.y-s[3],2)+1e-12,3/2)),
    s[0],
    -R*s[2] - C*s[3] + Q*( (a.y-s[3])/np.power(np.power(a.x-s[1],2)+np.power(a.y-s[3],2)+1e-12,3/2)+ (b.y-s[3])/np.power(np.power(b.x-s[1],2)+np.power(b.y-s[3],2)+1e-12,3/2)+ (c.y-s[3])/np.power(np.power(c.x-s[1],2)+np.power(c.y-s[3],2)+1e-12,3/2) ),
    s[2]])


#solve_oed static vars
def solve_oed():
    global sol_x,sol_y,d,s,t, dt, tmax
    t = np.linspace(0, tmax, num=np.round(tmax/dt)+1)
    _,sol_x,_,sol_y = odeint(d, s, t).T

def plot():
    global fig, ax, i, cStep, fft_axes, tmax
    fig = plt.figure()
    ax = Axes3D(fig)
    i=0
    cStep=int(0.4/(dt))#step piu' piccolo quando dt piu' grande
    while i < len(t.tolist()):
        ax.plot(xs=sol_x[i:i+cStep], ys=sol_y[i:i+cStep], zs=t[i:i+cStep],  color=plt.cm.jet(0.2+7*(np.sqrt(t.tolist()[i])/max(t.tolist()))))
        print(1+int(100*i/float(len(t.tolist()))))
        sys.stdout.flush()
        i=i+cStep-1
    ax.scatter(a.x,a.y, tmax, color=a.color, marker='.')
    ax.scatter(b.x,b.y, tmax, color=b.color, marker='.')
    ax.scatter(c.x,c.y, tmax, color=c.color, marker='.')
    ax.scatter(0, 0, tmax, marker='.')
    xlim([-35,35])
    ylim([-35,35])
    plt.show()

def asint():
    #define asint.
    ext_asint=np.array([np.mean(sol_x[-100:-1]),np.mean(sol_y[-100:-1])])
    #find closer
    a.dist=np.sqrt(np.power(a.x-ext_asint[0],2)+np.power(a.y-ext_asint[1],2))
    b.dist=np.sqrt(np.power(b.x-ext_asint[0],2)+np.power(b.y-ext_asint[1],2))
    c.dist=np.sqrt(np.power(c.x-ext_asint[0],2)+np.power(c.y-ext_asint[1],2))
    objs=np.array([a, b, c])
    dists=np.array([a.dist, b.dist, c.dist])
    closer = objs[np.where( dists == dists.min())[0][0]]
    #print "The reached point is "+closer.name
    #print "The associated color is "+closer.color
    if closer.color=='r':
        return (255,0,0)
    elif closer.color=='g':
        return (0,255,0)
    elif closer.color=='b':
        return (0,0,255)
p_drop(XXX,YYY)
solve_oed()
asint()
plot()

quit()
