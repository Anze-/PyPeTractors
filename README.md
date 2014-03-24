#PyPeTractors


Physics modeling of pendulum with three magnetic attractors in a gravity field in Python(mostly).

Denomination: Py[thon]Pe[ndulum][@]Tractors; read 'pie p ~~at~~ tractors'.

##NEWS
It looks like the odeint libraries for python are not the most efficent to solve a large number of differential equations (for a big set of initial conditions), moreover they are not develped to be vectorialized (and we would love to impoement GPU computing to create better-resolution pictures in a shorter time). This is why we have decided to use python with SciPy libraries to do the plots of the trajectories etc. while using C++ to create the map of the asintotic values given different initial condition, and we will probably script also a vectorialized CUDA script C++ script.

In a future we might wrap this whole package in a GUI software with related precompiled bins and/or exes.

####Well did you ask how is it faster?
We have computed the mean value of over 100 times measure of the same integration with the same resolution, the same pendulum model and the same  initial conditions using the same laptop CPU forcing the max speed to 80% so that the performance was the same. Obtaining a very little variance we measured that python took about 3 sec per integration, while C++ with the odeint library (still with no CUDA) took about 350 ms, which is almost 10 times faster.

![cpp looks faster](/plots/cpp_looks_faster.png)


##Plots

Up to now have plotted some trajecories:

![trajectory side view](/plots/magnetic_pendulum_1a.png)

x,y on the floor, time on the z, 3 attractors -> (r,g,b) dots, pendulum pin -> black dot
![trajectory top view](/plots/magnetic_pendulum_1b.png)


... and some maps:

![160x160 map](/plots/magnetic_pendulum_map_160x160.bmp)
to each initial position (x,y) of the picture is associated the color of the attractor reached by the pendulum
![40x40 map](/plots/magnetic_pendulum_RGB_map_40x40.bmp)


... using Python.

We will hopefully achive better results in the maps using C++, CUDA and more processing time.

##Sources

######Look at that crazy attractor:

http://en.wikipedia.org/wiki/Attractor

######Here are some infos:

http://bugman123.com/Fractals/index.html

######Now, are they just approximating a magnetic force using a linear equation? Seriously?... yep...:

http://www.codeproject.com/Articles/16166/The-magnetic-pendulum-fractal#idOverview

##Getting started:


######Python differential equations:

http://www.space-kerala.org/freeelectron/wp-content/uploads/2013/11/PythonScientific-simple.pdf

http://stackoverflow.com/questions/19779217/need-help-solving-a-second-order-non-linear-ode-in-python

http://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html

3d plot example:



    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = Axes3D(fig)
    x = [6,3,6,9,12,24]
    y = [3,5,78,12,23,56]
    # put 0s on the y-axis, and put the y axis on the z-axis
    ax.plot(xs=x, ys=[0]*len(x), zs=y, zdir='z', label='ys=0, zdir=z')
    plt.show()

Fixed axes example fro stackowerflow:

    fig = pylab.figure(figsize=(12,9))
    signal_axes = fig.add_subplot(211)
    signal_axes.plot(xs,rawsignal)
    fft_axes = fig.add_subplot(212)
    fft_axes.set_title("FFT")
    fft_axes.set_autoscaley_on(False)
    fft_axes.set_ylim([0,1000])
    fft = scipy.fft(rawsignal)
    fft_axes.plot(abs(fft))
    pylab.show()
    
    
    
##CREDITS & LICENSING

Author: Alberto Anzellotti

Co-Author: Giovanni Pederiva

Citations would be appreciated! : - ) Thanks.

ALL CONTENT RELEASED UNDER MIT LICENCE -  MORE @ http://anze.mit-license.org
