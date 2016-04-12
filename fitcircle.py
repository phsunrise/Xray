import numpy as np
from scipy.optimize import leastsq

def dist_func(var, x, y, w):
    A, D, theta = var
    
    P = (A * (x**2+y**2) 
       + np.sqrt(1+4*A*D)*(x*np.cos(theta)+y*np.sin(theta))
       + D)
    d = 2. * P / (1+np.sqrt(1+4*A*P))

    return d * np.sqrt(w) # weighted distance


def fitcircle(X, Y, weights): # parameters: xarray, yarray, weights
    if len(X) <= 3:
        raise ValueError("Array too short for circle fitting")

    # first calculate the initial guess x0, y0, r0 from three points
    # method for choosing the three points:
    #   if range in x is larger, then find points corresponding to
    #   xmin, xmax, and closest to 0.5*(xmin+xmax); else use y
    if np.max(X)-np.min(X) > np.max(Y)-np.min(Y):
        ind1 = np.argmax(X)
        ind2 = np.argmin(X)
        ind3 = np.argmin(np.abs(X-0.5*(np.min(X)+np.max(X))))
    else: 
        ind1 = np.argmax(Y)
        ind2 = np.argmin(Y)
        ind3 = np.argmin(np.abs(Y-0.5*(np.min(Y)+np.max(Y))))
    x1 = X[ind1]*1.
    x2 = X[ind2]*1.
    x3 = X[ind3]*1.
    y1 = Y[ind1]*1.
    y2 = Y[ind2]*1.
    y3 = Y[ind3]*1.
    M11 = np.linalg.det([[x1,y1,1],[x2,y2,1],[x3,y3,1]])
    M12 = np.linalg.det([[x1**2+y1**2,y1,1],
                         [x2**2+y2**2,y2,1],
                         [x3**2+y3**2,y3,1]])
    M13 = np.linalg.det([[x1**2+y1**2,x1,1],
                         [x2**2+y2**2,x2,1],
                         [x3**2+y3**2,x3,1]])
    M14 = np.linalg.det([[x1**2+y1**2,x1,y1],
                         [x2**2+y2**2,x2,y2],
                         [x3**2+y3**2,x3,y3]])
    x0 = 0.5*M12/M11
    y0 = -0.5*M13/M11
    r0 = np.sqrt(x0**2+y0**2+M14/M11)
    print "Initial guess: x0, y0, r0 =", x0, y0, r0

    A0 = 1. / (2.*r0)
    D0 = (x0**2 + y0**2 - r0**2) / (2.*r0)
    theta0 = np.arctan(y0*1./x0)
    if np.cos(theta0) * (-x0) < 0: # possible shift by pi
        theta0 += np.pi
    print A0, D0, theta0

    res, success = leastsq(dist_func, \
                [A0, D0, theta0], args=(X, Y, weights))
    print res

    r0 = 1./(2.*np.abs(res[0]))
    x0 = -np.sqrt(1+4*res[0]*res[1])*np.cos(res[2]) / (2*res[0])
    y0 = -np.sqrt(1+4*res[0]*res[1])*np.sin(res[2]) / (2*res[0])

    if success:
        return x0, y0, r0 
    else:
        raise RuntimeError("Fit not successful...")
