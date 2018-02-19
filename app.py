# -*- coding: utf-8 -*-

from sympy import symbols
from core import minimisation, maximisation

if __name__ == '__main__':
    X1, X2, X3, X4, X5 = symbols('X1 X2 X3 X4 X5')

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

    maximisation(arguments_2)

    print('\n\n\n\n')
    minimisation(arguments)


