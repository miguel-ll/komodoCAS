import cplot
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from numpy import sqrt,meshgrid,linspace,arange,pi,cos

def complexplot(eq, tup=(-3.0,3.0,400)):
    f = lambda x: eval(eq)
    plt = cplot.plot(f, tup, tup)
    plt.show()

def polarplot(eq,r=2):
    f = lambda x: eval(eq)
    plt.axes(projection='polar')
    rads = arange(0,(2*pi),0.01)
    for rad in rads:
        r = f(rad)
        plt.polar(rad,r,'g.')
    plt.show()

def funcplot(eq,num=100,s=3):
    f = lambda x: eval(eq)
    x = linspace(-s,s,num)
    fx = []
    for i in range(len(x)):
        fx.append(f(x[i]))
    plt.grid()
    plt.axvline()
    plt.axhline()
    plt.plot(x,fx)
    plt.show()

def vfplot2d(a="x",b="y"):
    x,y = meshgrid(linspace(-5,5,10),linspace(-5,5,10))
    u = eval(a)
    v = eval(b)
    plt.quiver(x,y,u,v)
    plt.show()
    
def plotvec2d(vecs,size=10,point=0):
    ax = plt.axes()
    for v in vecs:
        ax.arrow(0,0, v[0],v[1] , head_width=0.2, head_length=0.2)
    plt.ylim(point,size)
    plt.xlim(point,size)
    plt.show()

def plotvec3d(vecs, size=10):
    fig = plt.figure()
    ax = plt.axes(projection = '3d')
    ax.set_xlim([-size,size])
    ax.set_ylim([-size,size])
    ax.set_zlim([0,size])

    start = [0,0,0]
    for v in vecs:
        ax.quiver(start[0],start[1],start[2],v[0],v[1],v[2])
    plt.show()
    
def vfplot3d(a="x",b="y",c="z"):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x, y, z = meshgrid(arange(-0.8, 1, 0.2), arange(-0.8, 1, 0.2), arange(-0.8, 1, 0.8))
    u = eval(a)
    v = eval(b)
    w = eval(c)
    ax.quiver(x, y, z, u, v, w, length=0.1, color = 'black')
    plt.show()

#vec=[2,3]
#vec2=[6,2]
#vec3=[8,5]
#vecs = [vec,vec2,vec3]
#plotvec2d(vecs,size=12,point=-12)
#vecs = [u,v,q]
#plotvec3d(vecs)

#funcplot("x**5-5*x+3")
#complexplot("x**3+1")
#polarplot("cos(4*x)")
#vfplot2d(a="-y/(sqrt(x*x + y*y))", b="x/(sqrt(x*x + y*y))")
#vfplot3d(a="y",b="z",c="x")
