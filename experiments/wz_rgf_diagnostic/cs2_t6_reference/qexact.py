"""Exact complex rationals RC = re + i*im with re, im Fractions. Rejects bool, float, complex, numpy and mpmath values."""
from fractions import Fraction


def exact_rational(x, what="value"):
    """Return x as a Fraction if it is an int (not bool) or a Fraction; raise TypeError otherwise."""
    if isinstance(x, bool) or not isinstance(x, (int, Fraction)):
        raise TypeError(f"{what} must be int or Fraction (not bool/float/complex/numpy/mpmath), got {type(x).__name__}")
    return Fraction(x)


class RC:
    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        object.__setattr__(self, "re", exact_rational(re, "real part"))
        object.__setattr__(self, "im", exact_rational(im, "imaginary part"))

    def __setattr__(self, *_):
        raise AttributeError("RC is immutable")

    @staticmethod
    def of(x):
        """Coerce an RC, an exact rational, or a (re, im) pair of exact rationals."""
        if isinstance(x, RC):
            return x
        if isinstance(x, tuple) and len(x) == 2:
            return RC(x[0], x[1])
        return RC(x, 0)

    def __add__(self, o):
        o = RC.of(o)
        return RC(self.re + o.re, self.im + o.im)

    __radd__ = __add__

    def __sub__(self, o):
        o = RC.of(o)
        return RC(self.re - o.re, self.im - o.im)

    def __rsub__(self, o):
        return RC.of(o) - self

    def __neg__(self):
        return RC(-self.re, -self.im)

    def __mul__(self, o):
        o = RC.of(o)
        return RC(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = RC.of(o)
        n = o.re * o.re + o.im * o.im
        if n == 0:
            raise ZeroDivisionError("division by exact zero")
        return RC((self.re * o.re + self.im * o.im) / n, (self.im * o.re - self.re * o.im) / n)

    def conj(self):
        return RC(self.re, -self.im)

    def abs2(self):
        return self.re * self.re + self.im * self.im

    def is_zero(self):
        return self.re == 0 and self.im == 0

    def __eq__(self, o):
        if isinstance(o, RC):
            return self.re == o.re and self.im == o.im
        return NotImplemented

    def __hash__(self):
        return hash((self.re, self.im))

    def __repr__(self):
        return f"RC({self.re}, {self.im})"


ZERO = RC(0)
ONE = RC(1)
