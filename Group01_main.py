from chemistry import balancer, reaction_type

if __name__ == '__main__':

    try:
        equation = input("Enter a chemical equation: ")
        balanced_equation = balancer(equation)

        print(f"\nBalanced equation: {balanced_equation}")
        print(reaction_type(equation))

    except:
        print("This equation cannot be balanced")
