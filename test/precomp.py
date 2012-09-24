import numpy as np
import matplotlib as mp
import sympy as sp
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.animation as animation

# Principal lattice on the reference triangle
def princlat2(d):
    B = np.zeros(((d+2)*(d+1)/2,2))
    m = 0
    for l in range(d+1):
        for k in range(d-l+1):
            B[m,0] = float(k)/d
            B[m,1] = float(l)/d
            m+=1
    return B

# Lagrange basis on the principal lattice
def lagr2(d):
    B = princlat2(d)
    x, y = sp.symbols('x,y')
    bary = [sp.Poly(1-x-y, x, y), sp.Poly(x, x, y), sp.Poly(y, x, y)]
    lagr = [1 for k in range(B.shape[0])]
    for k in range(B.shape[0]):
        for p in range(3):
            barypk = float(bary[p].eval({x: B[k,0], y:B[k,1]}))
            if barypk > 0:
                for q in range(int(d*barypk)):
                    lagr[k] *= (bary[p] - float(q)/d) / (barypk - float(q)/d)
    return lagr   

def makesurf(ax, poly):
    X = np.arange(0, 1, 0.04)
    Y = np.arange(0, 1, 0.04)
    #X, Y = np.meshgrid(X, Y)

    Z = np.zeros((len(X),len(Y)))
    for i in range(len(X)):
        for j in range(len(Y)):
            Z[i,j] = float(poly.eval({x: X[i], y:Y[j]}))


    X, Y = np.meshgrid(X, Y)

    Zm = np.ma.masked_where(X+Y>1, Z)
    Z[np.where(np.ma.getmask(Zm)==True)] = np.nan
    
    lev = np.arange(0.0,1.0,0.01)
    norml = mp.colors.BoundaryNorm(lev, 200)

#    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet, norm = norml,
#                           linewidth=0, antialiased=False)
    surf = ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)


azim = 0 
def updatefig(num, fig, axes):
    global azim
    for ax in axes:
        ax.view_init(45, azim)
    azim += 5
    if azim <= 360:
        fig.savefig('pic/lagr23/' + str(azim/5) + '.png', dpi = (150))

x, y = sp.symbols('x,y')

d = 1

B = princlat2(d)
lagr = lagr2(d)

fig = plt.figure()

axes = []
#ffmpeg -sameq -i %d.png lagr1.mp4
for k in range(len(lagr)):
    ax = fig.add_subplot(1, 3, k+1, projection='3d', aspect='equal')
    
    ax.plot(B[:,0], B[:,1], np.zeros_like(B[:,0]), 'ro')
    makesurf(ax, sp.Poly(lagr[k],x,y))
    ax.set_zlim(0, 1.2)
    #ax.zaxis.set_major_locator(LinearLocator(10))
    #ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    ax.zaxis.set_ticks([])
    ax._axis3don = False

    axes.append(ax)
    plt.title(k)

fig.suptitle("Lagrange basis polynomials on the reference triangle (d=1)", fontsize=14)


line_ani = animation.FuncAnimation(fig, updatefig, 25, fargs=(fig, axes), interval=50, blit=False)


plt.show()
