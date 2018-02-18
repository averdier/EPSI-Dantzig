# -*- coding: utf-8 -*-

from sympy.solvers import solve


def eq_minus_var(eq, var):
    """
    Return new equation without variable

    :param eq:
    :param var:
    :return:
    """
    if var in eq.as_coefficients_dict():
        return eq - var * eq.coeff(var)
    return eq


def eq_minus_vars(eq, var_list):
    """
    Return new equation without variables

    :param eq:
    :param var_list:
    :return:
    """
    result = eq
    for var in var_list:
        result = eq_minus_var(result, var)
    return result


def subs_expr(expr, var, sub):
    """
    Substitute var in expression and return new expression

    :param expr:
    :param var:
    :param sub:
    :return:
    """
    result = {
        'expr': expr['expr'].subs(var, sub),
        'result': expr['result']
    }

    coeffs = result['expr'].as_coefficients_dict()
    # coeffs[1] is a number without symbols
    if coeffs.get(1, None):
        result['result'] -= coeffs[1]
        result['expr'] -= coeffs[1]

    return result


def find_ve(eq, mode='max'):
    """
    Find Ve in equation for minimisation

    If maximisation
    Ve = Var with Max { Coeff in expr }

    If minimisation
    Ve = Var with Min { Coeff in expr }

    :param eq:
    :param mode:
    :return:
    """
    s = list(eq.as_coefficients_dict().keys())
    ve = s[0] if s[0] != 1 else s[1]

    for var in s:
        if var != 1:
            if mode == 'max':
                if eq.coeff(var) > eq.coeff(ve):
                    ve = var

            if mode == 'min':
                if eq.coeff(var) < eq.coeff(ve):
                    ve = var

    return ve


def max_find_ve(eq):
    """
    Find Ve in equation for maximisation
    Ve = Var with Max { Coeff in expr }

    :param eq:
    :return:
    """
    s = list(eq.as_coefficients_dict().keys())
    ve = s[0] if s[0] != 1 else s[1]

    for var in s:
        if var != 1:
            if eq.coeff(var) > eq.coeff(ve):
                ve = var

    return ve


def min_find_ve(eq):
    """
    Find Ve in equation for minimisation
    Ve = Var with Max { Coeff in expr }

    :param eq:
    :return:
    """
    s = list(eq.as_coefficients_dict().keys())
    ve = s[0] if s[0] != 1 else s[1]

    for var in s:
        if var != 1:
            if eq.coeff(var) < eq.coeff(ve):
                ve = var

    return ve


def find_ve_value(exchange_expr, ve):
    """
    Find Ve value with exchange expression

    :param exchange_expr:
    :param ve:
    :return:
    """
    return solve(exchange_expr['expr'] - exchange_expr['result'], ve)[0]


def find_vs(exchange_expr, out_vars):
    """
    Find Vs in exchange expression
    <!> Call before remove Ve in out_var <!>

    :param exchange_expr:
    :param out_vars:
    :return:
    """
    return list(eq_minus_vars(exchange_expr['expr'], list(out_vars.keys())).as_coefficients_dict().keys())[0]


def find_exchange_expr(expr_list, ve):
    """
    Find the exchange expression
    The exchange expression have ve and Min { result / coeff(ve)}

    :param expr_list:
    :param ve:
    :return:
    """

    # Find first candidate
    result = None
    for expr in expr_list:
        if expr['result'] / expr['expr'].coeff(ve) > 0:
            result = expr
            break

    if result is not None:

        # Search minimal ratio
        for expr in expr_list:
            coeff = expr['expr'].coeff(ve)

            if expr['result'] / coeff > 0 and expr['result'] / coeff < result['result'] / result['expr'].coeff(ve):
                result = expr
    else:
        raise Exception('No exchange expression found')

    return result


def find_min_exchange_expr(expr_list, ve):
    """
    Find the exchange expression for minimisation
    The exchange expression have ve and Min { result / coeff(ve)}

    :param expr_list:
    :param ve:
    :return:
    """

    # Find first candidate
    result = None
    for expr in expr_list:
        if expr['result'] / expr['expr'].coeff(ve) > 0:
            result = expr
            break

    if result is not None:

        # Search minimal ratio
        for expr in expr_list:
            coeff = expr['expr'].coeff(ve)

            if expr['result'] / coeff > 0 and expr['result'] / coeff > result['result'] / result['expr'].coeff(ve):
                result = expr
    else:
        raise Exception('No exchange expression found')

    return result


def init_dantzig(args):
    """
    Do Dantzig initialization

    :param args:
    :return:
    """
    result = {
        'out_vars': {},
        'in_vars': {},
        'z': args['z'],
        'sc': args['sc']

    }

    # Z = 0 ==> Symbols Z = 0
    for var in args['z'].as_coefficients_dict():
        result['out_vars'][var] = 0

    # Symbols Z = 0 ==> Last symbol = result
    for eq in args['sc']:
        temp_eq = eq_minus_vars(eq['expr'], list(result['out_vars'].keys()))
        result['in_vars'][list(temp_eq.as_coefficients_dict().keys())[0]] = eq['result']

    return result


def dantzig_iteration(state, mode='max'):
    """
    Do dantzig iteration

    :param state:
    :param mode:
    :return:
    """
    ve = find_ve(state['z'], mode)
    exchange_expr = find_exchange_expr(state['sc'], ve)
    vs = find_vs(exchange_expr, state['out_vars'])

    coeff_ve = exchange_expr['expr'].coeff(ve)
    ve_value = find_ve_value(exchange_expr, ve)
    constraint_ve = {
        'expr': exchange_expr['expr'] / coeff_ve,
        'result': exchange_expr['result'] / coeff_ve
    }

    # Trick
    state['in_vars'][ve] = ve_value.as_coefficients_dict()[1]
    del state['out_vars'][ve]

    state['out_vars'][vs] = 0
    del state['in_vars'][vs]

    new_constraints = [constraint_ve]
    for expr in state['sc']:
        if expr != exchange_expr:
            new_expression = subs_expr(expr, ve, ve_value)
            new_constraints.append(new_expression)

            temp_eq = eq_minus_vars(new_expression['expr'], list(state['out_vars'].keys()))
            state['in_vars'][list(temp_eq.as_coefficients_dict().keys())[0]] = new_expression['result']

    new_z = state['z'].subs(ve, ve_value)

    print('#\tVE = {0}'.format(ve))
    print('#\tVS = {0}'.format(vs))
    print('#\tExchange expr = {0}'.format(exchange_expr))

    return {
        'in_vars': state['in_vars'],
        'out_vars': state['out_vars'],
        'z': new_z,
        'sc': new_constraints
    }


def vars_upper_zero(eq):
    s = list(eq.as_coefficients_dict().keys())
    for sym in s:
        if sym != 1 and eq.coeff(sym) < 0:
            return False
    return True


def vars_under_zero(eq):
    s = list(eq.as_coefficients_dict().keys())
    for sym in s:
        if sym != 1 and eq.coeff(sym) > 0:
            return False
    return True


def minimisation(problem):
    print("======================== Dantzig Minimisation =====================")
    print("# Inputs:")
    print("# \tZ  : {0}".format(problem['z']))
    print("# \tSC :")
    for e in problem['sc']:
        print('#\t\t{0}'.format(e))

    print('#\n# Initialization:')
    step = init_dantzig(problem)
    print('#\tOut vars: {0}'.format(step['out_vars']))
    print('#\tIn vars : {0}'.format(step['in_vars']))

    cpt = 0
    while not vars_upper_zero(step['z']):
        print('#\n# Iteration {0}'.format(cpt))
        step = dantzig_iteration(step, 'min')

        print('#\n#\tOut vars : {0}'.format(step['out_vars']))
        print('#\tIn vars  : {0}'.format(step['in_vars']))
        print('#\tZ : {0}'.format(step['z']))
        print('#\tSC :')
        for e in step['sc']:
            print('#\t\t{0}'.format(e))

        cpt += 1

    print("#\n====================================================================")
    print('#\n#\tResolved after : {0} iteration'.format(cpt))
    print('#\t Results: ')

    result = {}
    for s in list(problem['z'].as_coefficients_dict().keys()):
        if s != 1:
            if s in step['out_vars']:
                result[s] = step['out_vars'][s]
            else:
                result[s] = step['in_vars'][s]

            print('#\t\t{0} = {1}'.format(s, result[s]))

    print("#\n====================================================================")
    return result


def maximisation(problem):
    print("======================== Dantzig Maximisation =====================")
    print("# Inputs:")
    print("# \tZ  : {0}".format(problem['z']))
    print("# \tSC :")
    for e in problem['sc']:
        print('#\t\t{0}'.format(e))

    print('#\n# Initialization:')
    step = init_dantzig(problem)
    print('#\tOut vars: {0}'.format(step['out_vars']))
    print('#\tIn vars : {0}'.format(step['in_vars']))

    cpt = 0
    while not vars_under_zero(step['z']):
        print('#\n# Iteration {0}'.format(cpt))
        step = dantzig_iteration(step, 'max')

        print('#\n#\tOut vars : {0}'.format(step['out_vars']))
        print('#\tIn vars  : {0}'.format(step['in_vars']))
        print('#\tZ : {0}'.format(step['z']))
        print('#\tSC :')
        for e in step['sc']:
            print('#\t\t{0}'.format(e))

        cpt += 1

    print("#\n====================================================================")
    print('#\n#\tResolved after : {0} iteration'.format(cpt))
    print('#\t Results: ')

    result = {}
    for s in list(problem['z'].as_coefficients_dict().keys()):
        if s != 1:
            if s in step['out_vars']:
                result[s] = step['out_vars'][s]
            else:
                result[s] = step['in_vars'][s]

            print('#\t\t{0} = {1}'.format(s, result[s]))

    print("#\n====================================================================")
    return result
