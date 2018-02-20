# -*- coding: utf-8 -*-

from sympy import symbols
from core import minimisation, maximisation

if __name__ == '__main__':
    X1, X2, X3, X4, X5, X6 = symbols('X1 X2 X3 X4 X5 X6')

    arguments_2 = {
        'z': 3 * X1 + 2 * X2,
        'sc': [
            {'expr': 2 * X1 + X2 + X3, 'result': 18},
            {'expr': 2 * X1 + 3 * X2 + X4, 'result': 42},
            {'expr': 3 * X1 + X2 + X5, 'result': 24}
        ]
    }

    arguments = {
        'z': -8 * X1 + -6 * X2,
        'sc': [
            {'expr': 5 * X1 + 3 * X2 + X3, 'result': 30},
            {'expr': 2 * X1 + 3 * X2 + X4, 'result': 24},
            {'expr': 1 * X1 + 3 * X2 + X5, 'result': 18}
        ]
    }

    arguments_3 = {
        'z': X1 + 2 * X2,
        'sc': [
            {'expr': X1 + 3 * X2 + X3, 'result': 21},
            {'expr': -1 * X1 + 3 * X2 + X4, 'result': 18},
            {'expr': X1 + -1 * X2 + X5, 'result': 5},
        ]
    }

    arguments_4 = {
        'z': 6 * X1 + 7 * X2 + 8 * X3,
        'sc': [
            {'expr': X1 + 2 * X2 + X3 + X4, 'result': 100},
            {'expr': 3 * X1 + 4 * X2 + 2 * X3 + X5, 'result': 120},
            {'expr': 2 * X1 + 6 * X2 + 4 * X3 + X6, 'result': 200},
        ]
    }

    result = maximisation(arguments_2)
    max = result['z'].as_coefficients_dict()
    print(max)
    print('Max de Z = {0}'.format(max[1]))

    print('\n\n\n\n')
    minimisation(arguments)

    print('\n\n\n\n')
    maximisation(arguments_3)

    print('\n\n\n\n')
    result = maximisation(arguments_4)
    max = result['z'].as_coefficients_dict()
    print('Max de Z = {0}'.format(max[1]))




