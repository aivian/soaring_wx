import pdb

import numpy

import matplotlib.pyplot as plt

n = 40

z = numpy.array([0.0, 100.0e3, numpy.inf])
N = numpy.array([0.02047, 0.004188])
U = 10.0

z = numpy.array([0.0, 1.5e3, 3.0e3, numpy.inf])
N = numpy.array([0.022, 0.02, 0.0140])
U = numpy.array([5.0, 15.0, 23.0])

aij = lambda z, g: numpy.exp(z * g)
bij = lambda z, g: numpy.exp(-z * g)

def topo(x):
    """ The topography. We'll use the top from their example
    http://www.atmos.uw.edu/2010Q1/536/2003AP_lee_waves.pdf

    Arguments:
        x: horizontal coordinate (m)

    Returns:
        z: topographic height (m)
    """
    h0 = 250.0
    a = 1000.0
    xx = x - 25.0e3
    z = h0 * numpy.power(a, 2.0) / (numpy.power(xx, 2.0) + numpy.power(a, 2.0))
    #z = h0 * numpy.cos(xx * 2.0 * numpy.pi / a)
    return z

def dtopo_dx(x):
    """gradient of the topography)
    """
    h0 = 250.0
    a = 1000.0
    xx = x - 25.0e3
    dzdx = - (2.0 * numpy.power(a, 2.0) * h0 * xx) / numpy.power(
        numpy.power(xx, 2.0) + numpy.power(a, 2.0), 2.0)
    #dzdx = - h0 * 2.0 * numpy.pi / a * numpy.sin(xx * 2.0 * numpy.pi / a)
    u = numpy.interp(topo(x), z[:-1], U)
    return dzdx * u

def compute_fourier(f, T, N):
    f_sample = 2 * N
    t, dt = numpy.linspace(0, T, f_sample + 2, endpoint=False, retstep=True)
    y = numpy.fft.rfft(f(t)) / t.size
    y *= 2
    return y[0].real, y[1:-1].real, -y[1:-1].imag

coeffs_bc = compute_fourier(dtopo_dx, 50.0e3, n)

xtest = numpy.linspace(0, 50.0e3, 1000)
yy = numpy.zeros(xtest.shape)

for n_i, a, b in zip(range(n), coeffs_bc[1], coeffs_bc[2]):
    n_i += 1
    yy = yy + 0.0 * a * numpy.cos(2.0 * numpy.pi * xtest * n_i / 50.0e3) + b * numpy.sin(2.0 * numpy.pi * xtest * n_i / 50.0e3)

plt.plot(xtest, yy, xtest, dtopo_dx(xtest))
plt.show()

scorer = numpy.power(N/U, 2.0)
K = (numpy.arange(n) + 1) * 2.0 * numpy.pi / 50.0e3

X = []
XX = []
for k, a in zip(K, coeffs_bc[2]):
    M = numpy.array([a, 0.0, 0.0, 0.0], ndmin=2).T

    A = [[1, 1, 0, 0],]
    for idx in range(N.size-1):
        gamma_i = numpy.sqrt(numpy.power(k, 2.0) - scorer[idx] + 0j)
        gamma_ip1 = numpy.sqrt(numpy.power(k, 2.0) - scorer[idx+1] + 0j)
        zi = z[idx+1]
        aii = aij(zi, gamma_i)
        bii = bij(zi, gamma_i)
        aiip1 = aij(zi, gamma_ip1)
        biip1 = bij(zi, gamma_ip1)
        A.append([aii, bii, aiip1, biip1])
        A.append([gamma_i * aii, -gamma_i * bii, -gamma_ip1 * aiip1, gamma_ip1 * biip1])

    A.append([0, 0, numpy.inf, 0.0])

    A = numpy.array(A)
    X.append(numpy.linalg.solve(A, M))

Y = []
for zi in numpy.array([000.0]):
    y = numpy.zeros(xtest.shape)
    for x, k in zip(X, K):
        gamma = numpy.sqrt(numpy.power(k, 2.0) - scorer[0] + 0j)
        a = x[0] * aij(zi, gamma)
        b = x[1] * bij(zi, gamma)
        y = y + a * numpy.cos(xtest * k) + b * numpy.sin(xtest * k)
    Y.append(y)

plt.plot(xtest, y, xtest, dtopo_dx(xtest), xtest, yy)
plt.show()

xstart = xtest[0]
z_start = numpy.linspace(0.0, 4000.0, 10)
dt = 50.0
Zp = []
Xp = []
for zs in z_start:
    z_p = [zs]
    x_p = [xstart]
    while x_p[-1] < xtest[-1]:
        w = 0.0
        for x, k in zip(X, K):
            gamma = numpy.sqrt(numpy.power(k, 2.0) - scorer[0] + 0j)
            a = x[0] * aij(z_p[-1], gamma)
            b = x[1] * bij(z_p[-1], gamma)
            w = w + b * numpy.cos(x_p[-1]* k) + a * numpy.sin(x_p[-1] * k)
        z_p.append(z_p[-1] + dt * numpy.real(w))
        u = numpy.interp(z_p[-1], z[:-1], U)
        x_p.append(x_p[-1] + dt * u)
    plt.plot(x_p, z_p)
    Zp.append(z_p)
    Xp.append(x_p)

plt.plot(xtest, topo(xtest))
plt.show()
