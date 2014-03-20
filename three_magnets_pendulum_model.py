#Magnetic pendulum model solving nonlinear second order differential equations
#Code released under MIT LICENCE - more @ http://anze.mit-license.org/
#Author: Alberto Anzellotti
#co-Author: Giovanni Pederiva
#citation would be appreciated :) thanks and good coding!

# initial conditions: x(t_0) y(t_0) t_0 attractors_x=[x1, x2, x3] attractors_y=[y1, y2, y3]

import numpy as np
from scipy.integrate import odeint
from matplotlib.pylab import *
from math import *
from mpl_toolkits.mplot3d import Axes3D
from time import time


def p_drop(x_drop, y_drop):
    #standard model constants
    global R, C, Q, D, r, s, a, b, c
    R=.01		 #friction
    C=1*1e-2	 #gravity lin. approx. self restoring force
    Q=3*1e1	 #magnetic approx. charge force
    D=40	 #length of the pendulum
    r=30	 #distance of the magnets from the centre. (!smaller than D)
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
    #initial conditions
    s = np.array([
    0,# x'(t_0)
    x_drop,# x(t_0)
    0,# y'(t_0)
    y_drop# y(t_0)
    ])


def d_backup(s,t):
    return np.array([
    -R*s[0] - C*s[1] + Q* ((a.x-s[1])/np.power(np.power(a.x-s[1],2)+np.power(a.y-s[3],2)+D*D,3/2)+ (b.x-s[1])/np.power(np.power(b.x-s[1],2)+np.power(b.y-s[3],2)+D*D,3/2)+ (c.x-s[1])/np.power(np.power(c.x-s[1],2)+np.power(c.y-s[3],2)+D*D,3/2)),
    s[0],
    -R*s[2] - C*s[3] + Q*( (a.y-s[3])/np.power(np.power(a.x-s[1],2)+np.power(a.y-s[3],2)+D*D,3/2)+ (b.y-s[3])/np.power(np.power(b.x-s[1],2)+np.power(b.y-s[3],2)+D*D,3/2)+ (c.y-s[3])/np.power(np.power(c.x-s[1],2)+np.power(c.y-s[3],2)+D*D,3/2) ),
    s[2]])

def d(s,t):
    return np.array([
    -R*s[0] - C*s[1] + Q* ((a.x-s[1])/np.power(np.power(a.x-s[1],2)+np.power(a.y-s[3],2)+1e-12,3/2)+ (b.x-s[1])/np.power(np.power(b.x-s[1],2)+np.power(b.y-s[3],2)+1e-12,3/2)+ (c.x-s[1])/np.power(np.power(c.x-s[1],2)+np.power(c.y-s[3],2)+1e-12,3/2)),
    s[0],
    -R*s[2] - C*s[3] + Q*( (a.y-s[3])/np.power(np.power(a.x-s[1],2)+np.power(a.y-s[3],2)+1e-12,3/2)+ (b.y-s[3])/np.power(np.power(b.x-s[1],2)+np.power(b.y-s[3],2)+1e-12,3/2)+ (c.y-s[3])/np.power(np.power(c.x-s[1],2)+np.power(c.y-s[3],2)+1e-12,3/2) ),
    s[2]])

def solve_oed():
    global sol_x,sol_y,d,s,t, dt, tmax
    tmax = 150
    dt = .05
    t = np.linspace(0, tmax, num=np.round(tmax/dt)+1)
    _,sol_x,_,sol_y = odeint(d, s, t).T

def plot():
    global fig, ax, i, cStep, fft_axes, tmax
    fig = plt.figure()
    ax = Axes3D(fig)
    i=0
    cStep=int(0.4/(dt))#step piu' piccolo quando dt piu' grande
    while i < len(t.tolist()):
        i
        ax.plot(xs=sol_x[i:i+cStep], ys=sol_y[i:i+cStep], zs=t[i:i+cStep],  color=plt.cm.jet(0.2+7*(np.sqrt(t.tolist()[i])/max(t.tolist()))))
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

start = time()
p_drop(10,-10)
solve_oed()
asint()
elapsed = time() - start
plot()

import Image

img_q=40
img = Image.new( 'RGB', (img_q,img_q), "black") # create a new black image
pixels = img.load() # create the pixel map
t_ext = (elapsed*img_q*img_q*np.pi)/(4.0*60.0) #extimed total time in minutes
print "It will take about "+str(t_ext)+" minutes"
_=raw_input("press ENTER to proceed..")


i_row=0
while i_row < img_q:
    i_col=0
    while i_col < img_q:
        if np.sqrt(np.power(D*((-img_q/2)+i_row)/img_q,2)+np.power(D*((-img_q/2)+i_col)/img_q,2))<D/2 :
            p_drop(D*((-img_q/2)+i_col)/img_q,D*((-img_q/2)+i_row)/img_q)
            solve_oed()
            pixels[i_col, i_row]=asint()
        else:
            pixels[i_col, i_row]=(255,255,255)
	percent=i_row+i_col/float(100)
    	print str(percent)+" %"
        i_col=i_col+1
    i_row=i_row+1

img.show()
img.save("magnetic_pendulum_map.bmp")

"""


MORE PLOTTING
sol = odeint(d, s, t)
plt.plot(t, sol, label='lol')
plt.legend()
show()

plt.plot(sol_x, sol_y, label='X,Y dinamics')
plt.legend()

show()

plt.plot(t, sol_x, label='X,Y dinamics')
plt.legend()
show()

plt.plot(t, sol_y, label='X,Y dinamics')
plt.legend()
show()
"""
