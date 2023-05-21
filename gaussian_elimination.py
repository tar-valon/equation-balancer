# Perform gaussian elimination on an augmented matrix
# matrix to solve using gaussian elimination

def multiply_row(row, scalar):
    for number in range(len(row)):
        row[number] *= scalar
    return row


def divide_row(row, scalar):
    for number in range(len(row)):
        row[number] /= scalar
    return row


def subtract_row(row1, row2):
    for number in range(len(row2)):
        row2[number] = row1[number] - row2[number]
    return row2


def print_matrix(matrix):
    for row in range(len(matrix)):
        for column in range(len(matrix[0])):
            print(matrix[row][column], end='\t')
        print("\n")


def swap_row(aug_mat, r1, r2):
    tmp = aug_mat[r1].copy()
    aug_mat[r1] = aug_mat[r2][:]
    aug_mat[r2] = tmp[:]


def new_subtract_row(row_to_transform, row_to_subtract_from):
    for number in range(len(row_to_transform)):
        row_to_transform[number] = row_to_transform[number] - row_to_subtract_from[number]
    return row_to_transform


def gauss_jordan(aug_mat):

    row_start_position = 0
    column_count = len(aug_mat[0])
    row_count = len(aug_mat)
    stop_at_column = column_count
    stop_at_row = row_count

    # When rows is greater than columns
    if row_count > column_count:
        stop_at_column = column_count - 2
        stop_at_row = stop_at_column
    # When columns is greater than rows
    elif column_count > row_count:
        extra_columns = column_count - row_count
        stop_at_column = column_count - extra_columns
    # When both are same
    else:
        pass

    # Get the number of columns, does not include the final column which is the augmented part
    for column_number in range(stop_at_column):
        # Iterate through each row of the previous column
        for row_number in range(row_start_position, stop_at_row):
            # Check for diagonal element
            if row_number == column_number:
                # if diagonal element is 1, swap row with the next row
                # containing a non zero element from same column
                if aug_mat[row_number][column_number] == 0:
                    for i in range(row_number, stop_at_row):
                        if aug_mat[i][column_number] != 0:
                            swap_row(aug_mat, i, row_number)
                            scalar = aug_mat[row_number][column_number]
                            aug_mat[row_number] = divide_row(aug_mat[row_number], scalar)
                            break
                # Check if diagonal element is 1, if not then divide row by element to turn into 1
                elif aug_mat[row_number][column_number] != 1:
                    scalar = aug_mat[row_number][column_number]
                    aug_mat[row_number] = divide_row(aug_mat[row_number], scalar)
                else:
                    continue

            # Check if element is 0
            elif aug_mat[row_number][column_number] == 0:
                continue
            # Convert element to 0
            else:
                scalar = aug_mat[row_number][column_number]
                temp_row = aug_mat[row_start_position].copy()
                temp_row = multiply_row(temp_row, scalar)
                aug_mat[row_number] = subtract_row(temp_row, aug_mat[row_number])

        # The next row position should start from 1 row below(it should start from the next diagonal element)
        row_start_position += 1

   #print("\nMatrix in row echelon form")
   #print_matrix(aug_mat)

    # Continue to transform matrix into reduced row echelon form

    # Start 1 row above from the last diagonal element
    new_row_start_position = stop_at_row - 2
    # Start at the last row for temporary row, -1 because index starts at 0
    temp_row_position = stop_at_row - 1
    # Start from second last column to first column
    # Iterate from index of second last column to first column
    # -2 because of reverse iteration, len returns from 1 to max, we need index from 0 to max - 1
    # another - 1 because we want to exclude the augmented column
    for new_column_num in range(stop_at_column - 1, 0, -1):
        # stops at -1 to stop at row index 0, since stop number is excluded
        for new_row_num in range(new_row_start_position, -1, -1):
            if new_row_num == new_column_num:
                continue
            # Check if element is 0
            elif aug_mat[new_row_num][new_column_num] == 0:
                continue
            # Convert element to 0
            else:
                scalar = aug_mat[new_row_num][new_column_num]
                temp_row = aug_mat[temp_row_position].copy()
                temp_row = multiply_row(temp_row, scalar)
                aug_mat[new_row_num] = new_subtract_row(aug_mat[new_row_num], temp_row)

        # The next row position should start from 1 row below(it should start from the diagonal element)
        new_row_start_position -= 1
        # The next temporary row position should be 1 row above
        temp_row_position -= 1

    print("\nMatrix in reduced row echelon form")
    print_matrix(aug_mat)


def main(aug_mat):
    gauss_jordan(aug_mat)


if __name__ == "__main__":
    aug_mat = [[2, 4, 8, 2, 0],
               [1, 0, 3, 7, 7],
               [2, 1, 6, 8, 5],
               [1, 9, 7, 0, 6]]
    main(aug_mat)
