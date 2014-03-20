PyPeTractors
============

Phisics modeling of pendulum with three magnetic attractors in a gravity field in Python(mostly).

Sources
--------------
Look at that crazy atrtactor:

http://en.wikipedia.org/wiki/Attractor

Here are some infos:

http://bugman123.com/Fractals/index.html

Now, are they just approximating a magnetic force using a linear equation? Seriously?... yep...:

http://www.codeproject.com/Articles/16166/The-magnetic-pendulum-fractal#idOverview

Getting started:
--------------

Python differential equations:

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
