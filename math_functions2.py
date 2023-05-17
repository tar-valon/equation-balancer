from math import lcm
import fractions


def normalize(decimals):
    frac = [fractions.Fraction(number).limit_denominator() for number in decimals]
    frac_denominators = [number.denominator for number in frac]
    new_lcm = lcm(*frac_denominators)
    new_numbers = [int(number * new_lcm) for number in frac]
    return new_numbers


def row_of_zeros(row):
    is_zero = True
    for number in row:
        if number != 0:
            is_zero = False
            break
    return is_zero


if __name__ == '__main__':
    a = [0, 0, 0, 0, 0, 0]
    is_row_zero = row_of_zeros(a)
    print(is_row_zero)
