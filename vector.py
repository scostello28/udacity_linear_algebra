import numpy as np
from math import pi, acos, sqrt, cos
from decimal import Decimal, getcontext

getcontext().prec = 30

class Vector(object):

    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(self.coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')
            
    def __iter__(self):
            return iter(self.coordinates)

    def __getitem__(self,index):
            return self.coordinates[index]

    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)

    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def plus(self, v):
        new_coordinates = [x+y for x, y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)

    def minus(self, v):
        new_coordinates = [x-y for x, y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)

    def times_scalar(self, c):
        new_coordinates = [Decimal(c)*x for x in self.coordinates]
        return Vector(new_coordinates)

    def magnitude(self):
        return np.sqrt(np.sum([x**2 for x in self.coordinates]))

    def norm(self):
        try:
            magnitude = self.magnitude()
            return self.times_scalar(Decimal('1.0') / magnitude)
        except ZeroDivisionError:
            raise Exception(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)

    def dot(self, v):
        return np.sum([x*y for x,y in zip(self.coordinates, v.coordinates)])

    def angle(self, v, degrees=False):
        # # Mine works but his is more pythonic
        # rad = acos((self.dot(v)) / (self.magnitude() * v.magnitude()))
        # degrees = rad * 180. / pi
        # return rad, degrees
        try:
            u1 = self.norm()
            u2 = v.norm()
            angle_in_rad = acos(u1.dot(u2))

            if degrees:
                degrees_per_rad = 180. / pi
                return angle_in_rad * degrees_per_rad
            else:
                return angle_in_rad

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('Cannot compute an angle with the zero vector')
            else:
                raise e

    def parallel(self, v):
        return (self.zero() or
                v.zero() or
                self.angle(v) == 0 or
                self.angle(v) == pi)

    def orthogonal(self, v, tolerance=1e-10):
        return abs(self.dot(v)) < tolerance

    def zero(self, tolerance=1e-10):
        return self.magnitude() < tolerance

    def v_projection(self, b):
        try:
            u = b.norm()
            weight = self.dot(u)
            return u.times_scalar(weight)

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            else:
                raise e

    def v_orthogonal(self, b):
        try:
            return self.minus(self.v_projection(b))

        except Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(self.NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG)
            else:
                raise e

    def cross_product(self, w):
        if self.dimension == 3:
            v = self.coordinates
            w = w.coordinates
            one =    v[1] * w[2] - w[1] * v[2]
            two = -( v[0] * w[2] - w[0] * v[2] )
            three =  v[0] * w[1] - w[0] * v[1]
            return Vector([one, two, three])
        else:
            return 'Only works for 3 dimensional vectors'

    def area_of_parallelogram(self, w):
        cross_prod = self.cross_product(w)
        return cross_prod.magnitude()

    def area_of_triangle(self, w):
        return Decimal('0.5') * self.area_of_parallelogram(w)

if __name__ == '__main__':
    # print(my_vector)

    '''Plus, Minus, Scalar Multiply'''
    # my_vector = Vector([8.218, -9.341])
    # my_vector2 = Vector([-1.129, 2.111])
    # my_vector3 = Vector([7.119, 8.215])
    # my_vector4 = Vector([-8.223, 0.878])
    # my_vector5 = Vector([1.671, -1.012, -0.318])
    #
    # print(my_vector.plus(my_vector2))
    # print(my_vector3.minus(my_vector4))
    # print(my_vector5.times_scalar(7.41))

    '''Magnitude and Direction'''
    # my_vector6 = Vector([-0.221, 7.437])
    # my_vector7 = Vector([8.813, -1.331, -6.247])
    # my_vector8 = Vector([5.581, -2.136])
    # my_vector9 = Vector([1.996, 3.108, -4.554])
    #
    # print(my_vector6.magnitude())
    # print(my_vector7.magnitude())
    # print(my_vector8.norm())
    # print(my_vector9.norm())

    '''Dot Product and Angle'''
    # v1 = Vector([7.887, 4.138])
    # w1 = Vector([-8.802, 6.776])
    # v2 = Vector([-5.955, -4.904, -1.874])
    # w2 = Vector([-4.496, -8.755, 7.103])
    #
    # print(v1.dot(w1))
    # print(v2.dot(w2))
    #
    # v3 = Vector([3.183, -7.627])
    # w3 = Vector([-2.668, 5.319])
    # v4 = Vector([7.35, 0.221, 5.188])
    # w4 = Vector([2.751, 8.259, 3.985])
    #
    # print(v3.angle(w3))
    # print(v4.angle(w4, degrees=True))

    '''Parallel and Orthogonal Vectors'''
    # v1 = Vector([-7.579, -7.88])
    # w1 = Vector([22.737, 23.64])
    # v2 = Vector([-2.029, 9.97, 4.172])
    # w2 = Vector([-9.231, -6.639, -7.245])
    # v3 = Vector([-2.328, -7.284, -1.214])
    # w3 = Vector([-1.821, 1.072, -2.94])
    # v4 = Vector([2.118, 4.827])
    # w4 = Vector([0, 0])
    #
    # print('Parallel {} | Orthoganal {}'.format(v1.parallel(w1), v1.orthogonal(w1)))
    # print('Parallel {} | Orthoganal {}'.format(v2.parallel(w2), v2.orthogonal(w2)))
    # print('Parallel {} | Orthoganal {}'.format(v3.parallel(w3), v3.orthogonal(w3)))
    # print('Parallel {} | Orthoganal {}'.format(v4.parallel(w4), v4.orthogonal(w4)))

    '''Projecting Vectors'''
    # v1 = Vector([3.039, 1.879])
    # b1 = Vector([0.825, 2.036])
    # v2 = Vector([-9.88, -3.264, -8.159])
    # b2 = Vector([-2.155, -9.353, -9.473])
    # v3 = Vector([3.009, -6.172, 3.692, -2.51])
    # b3 = Vector([6.404, -9.144, 2.759, 8.718])
    #
    # print('Projection Component: {}'.format(v1.v_projection(b1)))
    # print('\n')
    # print('Orthogonal Component: {}'.format(v2.v_orthogonal(b2)))
    # print('\n')
    # print('Projection Component: {} | Orthogonal Component: {}'.format(v3.v_projection(b3), v3.v_orthogonal(b3)))

    '''Cross Product'''
    v1 = Vector([8.462, 7.893, -8.187])
    w1 = Vector([6.984, -5.975, 4.778])
    v2 = Vector([-8.987, -9.838, 5.031])
    w2 = Vector([-4.268, -1.861, -8.866])
    v3 = Vector([1.5, 9.547, 3.691])
    w3 = Vector([-6.007, 0.124, 5.772])

    print('Cross Product {}'.format(v1.cross_product(w1)))
    print('Area of parallelogram {}'.format(v2.area_of_parallelogram(w2)))
    print('Area of triangle {}'.format(v3.area_of_triangle(w3)))
