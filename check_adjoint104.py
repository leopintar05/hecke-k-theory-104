"""Rebuild the degree-104, level-three adjoint trace tests with Python 3.

In weight 53 a cusp basis is C*A^(47-3j)*B^j for j=0,...,15.  A is
the hexagonal theta series, B = eta(3z)^9/eta(z)^3, and
C = eta(z)^6*eta(3z)^6.  The first powers are q through q^16.
The weight-53 valence bound for Gamma_0(3) is 53/3, and a cusp form
also vanishes at the other cusp.  This bounds the dimension by 16.
"""

from math import gcd, isqrt, prod


def series_product(left, right):
    """Multiply truncated q-series of the same precision."""
    return [sum(left[i] * right[n - i] for i in range(n + 1))
            for n in range(len(left))]


def series_powers(series, highest):
    powers = [[1] + [0] * (len(series) - 1)]
    for _ in range(highest):
        powers.append(series_product(powers[-1], series))
    return powers


def euler_factor(series, step, exponent):
    """Multiply series in place by (1-q^step)^exponent."""
    for _ in range(abs(exponent)):
        if exponent > 0:
            for n in range(len(series) - 1, step - 1, -1):
                series[n] -= series[n - step]
        else:
            for n in range(step, len(series)):
                series[n] += series[n - step]


def cusp_basis():
    order = 7 * 16  # T_7 needs coefficients up to q^(7*16).
    theta = [0] * (order + 1)
    limit = 2 * isqrt(order) + 2  # x^2+xy+y^2 >= 3*max(x^2,y^2)/4
    for x in range(-limit, limit + 1):
        for y in range(-limit, limit + 1):
            norm = x*x + x*y + y*y
            if norm <= order:
                theta[norm] += 1
    b = [0, 1] + [0] * (order - 1)
    c = b.copy()
    for n in range(1, order + 1):
        euler_factor(b, n, -3)
        euler_factor(b, 3*n, 9)
        euler_factor(c, n, 6)
        euler_factor(c, 3*n, 6)
    theta_powers = series_powers(theta, 47)
    b_powers = series_powers(b, 15)
    return [series_product(c, series_product(theta_powers[47-3*j], b_powers[j]))
            for j in range(16)]


def hecke_matrix(forms, prime, character):
    """Compute T_p using a_pn + eta(p)*p^52*a_(n/p) in the triangular basis."""
    dimension = len(forms)
    for n, form in enumerate(forms, start=1):
        if form[:n] != [0] * n or form[n] != 1:
            raise ArithmeticError('Unexpected leading Fourier coefficient')
    columns = []
    for form in forms:
        coordinates = []
        for n in range(1, dimension + 1):
            value = form[prime*n]
            if n % prime == 0:
                value += character * prime**52 * form[n//prime]
            value -= sum(c * forms[j][n] for j, c in enumerate(coordinates))
            coordinates.append(value)
        columns.append(coordinates)
    return list(zip(*columns))


def square(matrix):
    columns = list(zip(*matrix))
    return [[sum(a*b for a, b in zip(row, col)) for col in columns]
            for row in matrix]


def add_scalar(matrix, scalar):
    return [[entry + scalar if i == j else entry for j, entry in enumerate(row)]
            for i, row in enumerate(matrix)]


def determinant(matrix):
    """Fraction-free Bareiss elimination with exact division and row pivoting."""
    a = [list(row) for row in matrix]
    size = len(a)
    previous, sign = 1, 1
    for k in range(size - 1):
        pivot_row = next((i for i in range(k, size) if a[i][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, size):
            for j in range(k + 1, size):
                value = pivot * a[i][j] - a[i][k] * a[k][j]
                a[i][j], remainder = divmod(value, previous)
                if remainder:
                    raise ArithmeticError('Inexact Bareiss division')
        previous = pivot
    return sign * a[-1][-1]


def main():
    forms = cusp_basis()
    t2 = hecke_matrix(forms, 2, -1)
    t7 = hecke_matrix(forms, 7, 1)
    # At inert 2 the field Frobenius is rational Frobenius squared.
    trace2 = add_scalar(square(t2), 2**53)
    d2 = determinant(add_scalar(square(trace2), -(1 + 2**104)**2))
    d7 = determinant(add_scalar(square(t7), -(1 + 7**52)**2))
    primes = (43279339111, 490440141803, 63421406191963401438831238997)
    common = gcd(d2, d7)
    expected = 3**76 * 5**12 * 11**2 * 19**2 * prod(primes)**2
    if common != expected:
        raise ArithmeticError(f'Unexpected determinant gcd: {common}')
    print('gcd =', common)
    print('gcd = 3^76 * 5^12 * 11^2 * 19^2 * (p1*p2*p3)^2')
    print('Large prime factors:', *primes)
    print('PASS degree-104 adjoint trace tests at 2 and 7 (dimension 16).')


if __name__ == '__main__':
    main()
