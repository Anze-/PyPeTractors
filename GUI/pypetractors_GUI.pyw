#!/usr/bin/python
# -*- coding: utf-8 -*-

"Sample gui attractors"

from __future__ import with_statement

__author__ = "Alberto Anzellotti and Giovanni Pederiva (alberto.anzellotti@gmail.com)"
__copyright__ = "Copyright (C) 2014- Anze"
__license__ = "MIT Licence"
__version__="0.0.1"

# images were taken from Pythoncard's proof and widgets demos

import datetime     # base imports, used by some controls and event handlers
import decimal
import time, math
import os, sys, subprocess, platform, signal

base_folder=os.getcwd()
import gui          # import gui2py package (shortcuts)

# set default locale to handle correctly numeric format (maskedit):
import wx, locale
#locale.setlocale(locale.LC_ALL, u'es_ES.UTF-8')
loc = wx.Locale(wx.LANGUAGE_DEFAULT, wx.LOCALE_LOAD_DEFAULT)


def open_folder(evt):
    global base_folder
    path=base_folder
    if platform.system() == "Windows":
        os.startfile(path)
    elif platform.system() == "Darwin":
        subprocess.Popen(["open", path])
    else:
        subprocess.Popen(["gnome-open", path])
progress=0

def plot_t(evt):
    global all
    try:
       p1.kill()
    except:
        pass
    p1=subprocess.Popen(["python","trajectories.py",str(mywin['notebook']['tab0']['XXX'].value),str(mywin['notebook']['tab0']['YYY'].value),str(mywin['notebook']['tab0']['R_box'].value),str(mywin['notebook']['tab0']['C_box'].value),str(mywin['notebook']['tab0']['Q_box'].value),str(mag_rad),str(max_time),str(t_delta)], stdout=subprocess.PIPE,bufsize=0)
    while p1.poll() is None:
        global progress
        progress = p1.stdout.readline()
        mywin['gauge'].value = int(progress)
        sys.stdout.write(progress)
        sys.stdout.flush()

def plot_m(evt):
    print("ok!")

# --- here go your event handlers ---

def button_press(evt):
    gui.alert("button clicked! %s" % evt.target.name, "Alert!!!")
    mywin['txtTest'].value = "hello world!!!!!"


def expand_item(event):
    "tree lazy evaluation example: virtually add children at runtime"
    tv = event.target
    if not event.detail.get_children_count():
        for i in range(5):
            it = tv.items.add(parent=event.detail, text="lazy child %s" % i)
            it.set_has_children()  # allow to lazy expand this child too

def t0_slider0_click(evt):
    # move the progress bar according the slider ("scroll bar")
    global max_time
    print("Slider value: %s" % (evt.target.value*5))
    max_time=5*mywin['notebook']['tab0']['slider'].value
    #mywin['gauge'].value = slider_val
    mywin['notebook']['tab0']['slider_value'].text=str(max_time)

def t0_slider1_click(evt):
    # move the progress bar according the slider ("scroll bar")
    print("Slider value: %s" % evt.target.value)
    global mag_rad
    mag_rad=mywin['notebook']['tab0']['slider1'].value/2
    #mywin['gauge'].value = slider_val
    mywin['notebook']['tab0']['slider1_value'].text=str(mag_rad)

def t0_slider2_click(evt):
    # move the progress bar according the slider ("scroll bar")
    global t_delta
    print("Slider value: %s" % evt.target.value)
    t_delta=math.exp(-2.2+(mywin['notebook']['tab0']['slider2'].value)/15.0)/300.0
    #mywin['gauge'].value = slider_val
    mywin['notebook']['tab0']['slider2_value'].text=str(t_delta)

def t1_slider0_click(evt):
    # move the progress bar according the slider ("scroll bar")
    global i_rad
    print("Slider value: %s" % evt.target.value)
    i_rad=mywin['notebook']['tab1']['slider'].value
    #mywin['gauge'].value = slider_val
    mywin['notebook']['tab1']['slider_value'].text=str(i_rad)

def t1_slider1_click(evt):
    # move the progress bar according the slider ("scroll bar")
    global i_dens
    i_dens=mywin['notebook']['tab1']['slider11'].value/2.0
    mywin['notebook']['tab1']['slider11_value'].text=str(i_dens)

def help(evt):
    msg = []
    msg.append("""Trajectories:
    N mags - choose number of magnets
    Color line - the line changes color accordingly to the time
    3D - plot the time on the z-axis
    Pin dot - plot blak dot which represents the pendulum pin spot
    R - air friction
    C - gravity force intensity
    Q - magnetic charge intensity
    End Time, Time delta - bigger means better but slower
    Mag radius - Distance from each magnet to the pin
Map:
    Image radius - distance from the pin to the side of the map
    Image density - square root of calculations per map square unit.
Buttons:
->Plot Trajectory - given x,y initial conditions will plot the x,y  positions of the pendulum over time
->Plot Map - will plot a map each pixel of which represents with a colour the magnet reached with the initial conditions given by the pixel coordinates in the map
    """)
    gui.alert('\n'.join(msg), "Help", scrolled=True)

def credits(evt):
    msg = []
    msg.append("""                                                    Thanks to:
Alberto Anzellotti and Giovani Pederiva, University of Trento,         Degree in Physics.

All the code is licensed under MIT License, see more at http://anze.mit-license.org, citations would be appreciated!

Featuring GitHub Education, Python, gui2py, SciPy, C++, Odeint, OpenMP, Linux and Gnome.


""")
    gui.alert('\n'.join(msg), "Credits and more", scrolled=True)

def image_buton_pressed(evt):
    if platform.system() == "Windows":
        os.startfile(map_path)
    elif platform.system() == "Darwin":
        subprocess.Popen(["open", map_path])
    else:
        subprocess.Popen(['gnome-open',map_path])


t0 = time.time()    # user for basic timing

i_rad=80

max_time=150

t_delta=0.05

mag_rad=30

i_dens=2.0

map_path="map.ppm"

#helps

h_3D=u'Plot time on the z axis'


# --- gui2py designer generated code starts ---

with gui.Window(name='mywin', title=u'gui2py sample app', resizable=False,
                height='565px', left='180', top='24', width='396px',
                bgcolor=u'#E0E0E0', image='tile.bmp', tiled=True, ):
#    gui.ListBox(name='lstTest', height='96', left='277', top='28',
#                width='103', data_selection='two',
#                items=[u'one', u'two', u'tree'], selection=1,
#                string_selection=u'two', )
    gui.Label(name='lblTest', alignment='right', transparent=True, left='38',
              top='30', width='350', height="50", text=u'Welcome to PyPeTractors, with this app you will be able to run a few magnetic pendulum simulations ', )
#    gui.TextBox(name='txtTest', left='100', top='31', width='152',
#               value=u'mariano', )
    gui.Line(name='line_25_556', height='3', left='25', top='80',
             width='349', )
    with gui.Notebook(name='notebook', height='320', left='21', top='85',
                      width='360', selection=0, ):
        with gui.TabPanel(id=133, name='tab0', selected=True, text=u'Trajetories', ):
            with gui.Panel(label=u'Options: ', name='panel', height='145',
                           left='15', top='10', width='110', image='', ):
                gui.RadioButton(id=274, label=u'3 mags', name='mag3',
                                left='14', top='23', width='85', value=True, )
                gui.CheckBox(label=u'Colr line', name='color', left='14',
                             top='50', value=True, )
                gui.CheckBox(label=u'3 D', name='threeD', left='14',
                             top='85', value=True, )
                gui.CheckBox(label=u'Pin dot', name='pin', left='14',
                             top='120', value=True,)

            gui.Label(name='label_159_27', left='147', top='23',
                      text=u'init X:', )
            gui.TextBox(mask='****.***', name='XXX',
                        alignment='right', left='220', top='19', width='110',
                        value="-010.000", )
            gui.Label(name='label_153_56', left='147', top='60',
                      text=u'init Y:', )
            gui.TextBox(mask='****.***',allowNegative=True, name='YYY', alignment='right',
                        left='220', top='55', width='110', value=u'0010.000', )
            gui.Label(name='R_label', left='155', top='92',
                      text=u'R:', )
            gui.TextBox(mask='##.####', name='R_box', left='132', top='115',
                        width='70', value=0.01, )
            gui.Label(name='C_label', left='225', top='92',
                      text=u'C:', )
            gui.TextBox(mask='##.####', name='C_box', left='202', top='115',
                        width='70', value=0.01, )
            gui.Label(name='Q_label', left='295', top='92',
                      text=u'Q:', )
            gui.TextBox(mask='###.###', name='Q_box', left='272', top='115',
                        width='70', value=30,)
            gui.Label(name='slider_label', left='10', top='175',text='End Time', )
            gui.Slider(name='slider', left='20', top='190', width='240', freq=0, value=max_time*0.2)
            gui.Label(name='slider_value', left='300', top='189',text=str(max_time), )
            gui.Label(name='slider2_label', left='10', top='210',text='Time delta', )
            gui.Slider(name='slider2', left='20', top='225', width='240', freq=0, value=74,)
            gui.Label(name='slider2_value', left='300', top='225',text=str(t_delta), )
            gui.Label(name='slider1_label', left='10', top='240',text='Mag radius', )
            gui.Slider(name='slider1', left='20', top='255', width='240', freq=0, value=mag_rad*2,)
            gui.Label(name='slider1_value', left='300', top='255',text=str(mag_rad), )


        with gui.TabPanel(id=163, name='tab1', selected=False, text=u'Map', ):
            gui.Button(name='image', height='150', left='100', top='10',width='150', filename="sample.bmp",)
            gui.Label(name='slider_label', left='10', top='182',text='Image radius', )
            gui.Slider(name='slider', left='20', top='200', width='240', freq=0, value=i_rad)
            gui.Label(name='slider_value', left='300', top='200',text=str(i_rad), )

            gui.Label(name='slider11_label', left='10', top='232',text='Image density', )
            gui.Slider(name='slider11', left='20', top='250', width='250', freq=0, value=i_dens*2)
            gui.Label(name='slider11_value', left='300', top='250',text=str(i_dens), )
            #            with gui.GridView(name='gridview', height='100%', left='0',
            #                              top='0', width='100%', ):
            #                gui.GridColumn(name='col1', text='Col A', type='text',
            #                               width=75, )
            #                gui.GridColumn(name='col2', text='Col 2', type='long',
            #                              width=75, )
            #                gui.GridColumn(name='col3', text='Col B', type='float',
            #                               width=75, )
            pass
    with gui.MenuBar(name='menubar', fgcolor=u'#000000', ):
        with gui.Menu(name='menu', fgcolor=u'#000000', ):
            gui.MenuItem(label=u'Show Folder', id=142, name='show',help=base_folder,)
            gui.MenuItemSeparator(id=130, name='menuitemseparator_130', )
            gui.MenuItem(label=u'Quit', id=140, name='quit', onclick='exit()',help=u'Quit the program',)
        with gui.Menu(label=u'About', name='about', ):
            gui.MenuItem(help=u'Help', id=225,
                         label=u'Help', name='help', )
            gui.MenuItem(help=u'Infos about the program',
                         id=235, label=u'Credits and more', name='credits', )
    gui.Button(label=u'Plot Trajectory', name='plt_t', left='40', top='420',
               width='220', height='50',)
    gui.Button(label=u'Plot Map', name='plt_m', left='280', top='420',
               width='85', height='50', )
    gui.Gauge(name='gauge', height='18', left='13', top='480', width='367',
              value=0, )

    gui.StatusBar(name='statusbar', )
    #gui.TreeView(name='treeview', default_style=True, has_buttons=True,
    #             height='98', left='223', top='312', width='154',
    #             onitemselected="print('selected TreeItem: %s' % event.detail.text)", )

# --- gui2py designer generated code ends ---


t1 = time.time()
print("basic creation timing: t1 - t0: %s" % (t1 - t0))

# get a reference to the Top Level Window (used by designer / events handlers):
mywin = gui.get("mywin")

# assign (bind) some event handlers
# (move to "__main__" block if you don't want to be executed in design mode):

#mywin.onload =
#mywin['btnTest'].onclick = button_press
mywin['notebook']['tab0']['slider'].onclick = t0_slider0_click
mywin['notebook']['tab0']['slider1'].onclick = t0_slider1_click
mywin['notebook']['tab0']['slider2'].onclick = t0_slider2_click
mywin['notebook']['tab1']['slider'].onclick = t1_slider0_click
mywin['notebook']['tab1']['slider11'].onclick = t1_slider1_click
mywin['menubar']['menu']['show'].onclick=open_folder
mywin['menubar']['about']['help'].onclick=help
mywin['menubar']['about']['credits'].onclick=credits
mywin['plt_t'].onitemexpanding = expand_item
mywin['plt_m'].onitemexpanding = expand_item
#mywin['treeview'].onitemexpanding = expand_item
#mywin['menubar']['list']['del'].onclick = del_an_item
#mywin['menubar']['list']['update'].onclick = update_items
#mywin['menubar']['grid']['add'].onclick = add_a_row
#mywin['menubar']['grid']['del'].onclick = del_sel_rows
#mywin['menubar']['grid']['clear'].onclick = clear_rows
#mywin['menubar']['grid']['update'].onclick = update_rows
mywin['plt_t'].onclick = plot_t
mywin['plt_m'].onclick = plot_m
mywin['notebook']['tab1']['image'].onmousedown = image_buton_pressed

if __name__ == "__main__":
    mywin.show()
    mywin.title = "PyPeTractors "+__version__
    # update the status bar text with python and wxpython version info:
    import wx, sys, platform
    mywin['statusbar'].text = "Developed by Alberto Anzellotti and Giovanni Pederiva"
    gui.main_loop()
