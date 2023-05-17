import chemparse as cp
import re
import gaussian_elimination as rref
from math_functions2 import normalize, row_of_zeros


def reaction_type(equation):
    reactant_compounds, product_compounds, _ = equation_parser(equation)
    if "O2" in reactant_compounds and ("C" in product_compounds or "CO" in product_compounds):
        return "This is an incomplete combustion reaction"
    elif "O2" in reactant_compounds and ("CO2" in product_compounds or "H2O" in product_compounds):
        return "This is a complete combustion reaction"
    else:
        return "This is a non combustion reaction"


def equation_parser(equation):
    # Split equation into 2 parts
    reactants, products = re.split(r'=|->', equation)

    # Split each part into multiple reactants and products
    reactant_compounds = re.split(r'[+]', reactants)
    product_compounds = re.split(r'[+]', products)

    # Remove whitespace from compunds
    reactant_compounds = [compound.strip() for compound in reactant_compounds]
    product_compounds = [compound.strip() for compound in product_compounds]

    # Create list with Parsed elements
    compounds_dict = []
    for compound in reactant_compounds:
        parsed_compound = cp.parse_formula(compound)
        compounds_dict.append(parsed_compound)

    # Create list with Parsed elements, negative
    # value because the product coefficients must
    # move to the other side of the equation
    for compound in product_compounds:
        parsed_compound = cp.parse_formula(compound)
        parsed_compound = {key: -value for key, value in parsed_compound.items()}
        compounds_dict.append(parsed_compound)

    return reactant_compounds, product_compounds, compounds_dict


# Function to check if equation is balanced
# Gets a total for each element for both reactants and products
# Then compares the two dictionaries to see if they are equal
def equation_checker(compounds_dict, unique_elements):
    product_dict = {}
    reactant_dict = {}
    for element in unique_elements:
        if element not in product_dict:
            reactant_dict[element] = 0
            product_dict[element] = 0
        for compound in compounds_dict:
            if element not in compound:
                continue
            elif int(compound[element]) > 0:
                reactant_dict[element] = int(reactant_dict[element]) + int(compound[element])
            else:
                product_dict[element] = int(product_dict[element]) + (-1 * int(compound[element]))

    if product_dict == reactant_dict:
        return True
    else:
        return False


def balancer(equation):
    reactant_compounds, product_compounds, compounds_dict = equation_parser(equation)

    # Find unique elements
    unique_elements = []
    for compound in compounds_dict:
        for key in compound:
            if key not in unique_elements:
                unique_elements.append(key)

    # Check if equation is balanced, if so then return to main
    is_balanced = equation_checker(compounds_dict, unique_elements)
    if is_balanced:
        return "Equation is already balanced"

    # Initialize rows and column for
    # the augmented matrix to solve.
    # +1 is for the column with 0's
    matrix_columns = len(compounds_dict) + 1
    matrix_rows = len(unique_elements)

    # Create an empty matrix
    elements_matrix = []
    for i in range(matrix_rows):
        # Append an empty sublist inside the list
        elements_matrix.append([])
        for j in range(matrix_columns):
            elements_matrix[i].append(0)

    # Fill matrix with coefficients
    # len is -1 because the last column should
    # contain only 0's
    # zip() used for parallel iteration
    for (column, compound) in zip(range(len(elements_matrix[0]) - 1), compounds_dict):
        for (row, key) in zip(range(len(elements_matrix)), unique_elements):
            if key in compound:
                elements_matrix[row][column] = compound[key]
            else:
                elements_matrix[row][column] = 0

    print("\nAugmented matrix of the chemical equation")
    rref.print_matrix(elements_matrix)

    # Solve the matrix
    rref.gauss_jordan(elements_matrix)

    # Check if the solved matrix contains
    # a row of zeros, if so remove items
    for row_index in range(len(elements_matrix)):
        contains_only_zeros = row_of_zeros(elements_matrix[row_index])
        if contains_only_zeros:
            del elements_matrix[row_index]

    final_coefficients = []

    # last column does not count the augmented part
    if len(elements_matrix[0]) == (len(elements_matrix) + 1):
        last_column = len(elements_matrix[0]) - 1
    # If there is a free variable, the free variable
    # column and augmented column are not counted
    elif len(elements_matrix[0]) == (len(elements_matrix) + 2):
        last_column = len(elements_matrix[0]) - 2

    # get the solutions from solved matrix and store it into
    # final_coefficients after multipltying by -1 and converting
    # decimals to fractions, -1 becuase the coefficients move to the
    # other side of the equation. Also checks if the coefficient is an integer

    for row_number in range(len(elements_matrix)):
        elements_matrix[row_number][last_column] *= -1
        final_coefficients.append(elements_matrix[row_number][last_column])

    # Add the free variable to final_coefficients
    final_coefficients.append(1)

    new_coefficients = normalize(final_coefficients)
    print(new_coefficients)

    # Join coefficients with compounds
    balanced_compound = []
    for (i, j) in zip(new_coefficients, reactant_compounds+product_compounds):
        balanced_compound.append(str(i)+j)

    # Join products and reactants
    balanced_equation = ""
    balanced_equation += ' + '.join(balanced_compound[:len(reactant_compounds)])
    balanced_equation += ' -> '
    balanced_equation += ' + '.join(balanced_compound[len(reactant_compounds):])

    return balanced_equation
