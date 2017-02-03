import numpy

import matplotlib.pyplot as plt

aij = lambda z, g: numpy.exp(z * g)
bij = lambda z, g: numpy.exp(-z * g)

n = 10

def topo(x):
    """ The topography. We'll use the top from their example
    http://www.atmos.uw.edu/2010Q1/536/2003AP_lee_waves.pdf

    Arguments:
        x: horizontal coordinate (m)

    Returns:
        z: topographic height (m)
    """
    h0 = 573.0
    a = 10000.0
    x -= 25.0e3
    z = h0 * numpy.power(a, 2.0) / (numpy.power(x, 2.0) + numpy.power(a, 2.0))
    return z

def dtopo_dx(x):
    """gradient of the topography)
    """
    h0 = 573.0
    a = 10000.0
    x -= 25.0e3
    dzdx = - (2.0 * numpy.power(a, 2.0) * h0 * x) / numpy.power(
        numpy.power(x, 2.0) + numpy.power(a, 2.0), 2.0)
    return dzdx * 10.0

def compute_fourier(x, y, N):
    a = []
    b = []
    P = x[-1] - x[0]
    for n in range(N+1):
        a.append(2.0 / P *
            numpy.trapz(y * numpy.cos(2.0 * numpy.pi * n * x / P), x))
        b.append(2.0 / P *
            numpy.trapz(y * numpy.sin(2.0 * numpy.pi * n * x / P), x))

xtest = numpy.linspace(0, 50.0e3, 1000)

ybc = topo(xtest)
coeffs_bc = compute_fourier(xtest, ybc, n)

y = numpy.zeros(xtest.shape)

for n_i, a, b in zip(range(n), coeffs_bc[1], coeffs_bc[2]):
    n_i += 1
    y = y + a * numpy.cos(2.0 * numpy.pi * xtest * n_i / 50.0e3) + b * numpy.cos(2.0 * numpy.pi * xtest * n_i / 50.0e3)
y += coeffs_bc[0]

plt.plot(xtest, -y * 1000.0, xtest, topo(xtest))
plt.show()
