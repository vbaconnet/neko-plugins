import numpy as np

def S(x):
    if (x <= 0.0):
        return 0.0
    elif (x >= 1.0):
        return 1.0
    else:
        return 1.0 / ( 1.0 + np.exp( 1.0/(x-1.0) + 1.0/x  ))

def fringe(x, xs, xe, xmin, xmax, dr, df):
    if (x < xmin): return 0.0
    elif (x > xmax): return 0.0
    else:
        return S((x - xs)/dr) - S((x-xe)/df + 1.0)

xmin =-0.1
xmax = 0.1
alpha= 0.15

#xs =-0.075
xe = xmax
xs =xmin
#xe = 0.075
dr = (xmax-xmin)*alpha
df = dr

x = np.linspace(xmin, xmax, 501)
y = np.zeros_like(x)

for i in range(len(x)):
    y[i] = fringe(x[i], xs, xe, xmin, xmax, dr, df)

import matplotlib.pyplot as plt

def line(x):
    plt.plot([x,x], [0.0, 1.2], "k--", alpha=0.5)

plt.plot(x,y)
plt.xlim(xmin, xmax)
plt.ylim(0.0,1.1)
plt.xlabel("u")
plt.ylabel(r"$\lambda_u$")
plt.xticks([xmin, xs, xs+dr, xe-df, xe, xmax], [r"$u_{min}$", r"$u_{start}$", r"$u_{start}+\delta_{rise}$", r"$u_{end} - \delta_{fall}$", r"$u_{end}$", r"$u_{max}$"])
line(xmin)
line(xmax)
line(xs)
line(xs+dr)
line(xe)
line(xe-df)
line(xmax)
plt.tight_layout()
plt.savefig("fringe.png")
plt.show()
