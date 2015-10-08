import numpy as np
import matplotlib.pyplot as plt
import os

#PLOTTING PARAMETERS
savedir = 'results/'
#Save filetype
savetype = '.png'
#Plot Colors
colors = ['green', 'red', 'blue', 'cyan']
#Plot Markers
markers = ['.', 'x', '*', 'o']
#Line Styles
mark = 5
minimark = 0.75
line = 1.5
#Font Styles
font_ttl = {'family' : 'serif',
            'color'  : 'black',
            'weight' : 'normal',
            'size'   : 18,
            }
font_lbl = {'family' : 'serif',
            'color'  : 'black',
            'weight' : 'normal',
            'size'   : 18,
            }
font_box = {'family' : 'arial',
            'color'  : 'black',
            'weight' : 'normal',
            'size'   : 12,
            }
font_tick = 16

def PlotStart(title, xlbl, ylbl, horzy='vertical'):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.title(title, fontdict=font_ttl)
    plt.xlabel(xlbl, fontdict=font_lbl)
    plt.xticks(fontsize=font_tick)
    plt.ylabel(ylbl, fontdict=font_lbl, rotation=horzy)
    plt.yticks(fontsize=font_tick)
    #increase title spacing
    ttl = ax.title
    ttl.set_position([.5, 1.025])
    return fig, ax

def MakeOutputDir(savedir):
    """make results output directory if it does not already exist.
    instring --> directory path from script containing folder
    """
    #split individual directories
    splitstring = savedir.split('/')
    prestring = ''
    for string in splitstring:
        prestring += string + '/'
        try:
            os.mkdir(prestring)
        except Exception:
            pass

def GetParentDir(savename):
    """Get parent directory from path of file"""
    #split individual directories
    splitstring = savename.split('/')
    parent = ''
    #concatenate all dirs except bottommost
    for string in splitstring[:-1]:
        parent += string + '/'
    return parent

def SavePlot(savename, overwrite=1):
    """Save file given save path.  Do not save if file exists
    or if variable overwrite is 1"""
    if os.path.isfile(savename):
        if overwrite == 0:
            print('     Overwrite is off')
            return
        else: os.remove(savename)
    MakeOutputDir(GetParentDir(savename))
    plt.savefig(savename, bbox_inches='tight')

def TextBox(ax, boxtext, x=0.005, y=0.95, fontsize=font_box['size'],
                                                    alpha=0.5, props=None):
    if props == None:
        props = dict(boxstyle='round', facecolor='white', alpha=alpha)
    ax.text(x, y, boxtext, transform=ax.transAxes, fontsize=fontsize,
            verticalalignment='top', bbox=props)

########################################################################
### LOAD DATA ##########################################################
########################################################################
x = np.loadtxt('x.dat', delimiter='\n')
y = np.loadtxt('y.dat', delimiter='\n')

#FIT DATA
m = 1.84419537
b = 0.190936923
r = 0.948222518
def lineeqn(x, m, b):
    return m * x + b
xfit = np.linspace(min(x), max(x))
yfit = []
for xx in xfit:
    yfit.append(lineeqn(xx, m, b))
print(xfit, yfit)

#PlOT DATA
title = 'Least-Squares Fit'
xlbl = 'x'
ylbl = 'y'
_, ax = PlotStart(title, xlbl, ylbl, 'horizontal')
ax.plot(x, y, label='Original Data', color='blue', linestyle='', linewidth=line,
                                        marker='.', markersize=mark*2)
ax.plot(xfit, yfit, label='Least-Squares Fit', color='red',
                                        linestyle='-', linewidth=line)
boxtext = 'Least-Squares Fit:\n\n$y={}x+{}$\nR={}'.format(m, b, r)
TextBox(ax, boxtext,)
plt.legend(loc='lower right', fancybox=True, framealpha=0.5)
plt.show()
# plt.savefig('PlotFit.png', bbox_inches='tight')
SavePlot('PlotFit.png')

