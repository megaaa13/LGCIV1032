"""
mu_to_omega.py

Author : Luca Hafhouf

Permet de retrouver aisément le omega à partir de la valeur de mu.

"""

value_ref = [(0.01, 0.009), (0.02, 0.019), (0.03, 0.029), (0.04, 0.038), (0.05, 0.048), (0.06, 0.058), (0.07, 0.069), (0.08, 0.080), (0.09, 0.091), (0.1, 0.102), (0.11, 0.114), (0.12, 0.125), (0.13, 0.137), (0.14, 0.149), (0.15, 0.161), (0.16, 0.173), (0.17, 0.185), (0.18, 0.198), (0.19, 0.211), (0.2, 0.224), (0.21, 0.237), (0.22, 0.251), (0.23, 0.264), (0.24, 0.278), (0.25, 0.293), (0.26, 0.307), (0.27, 0.322), (0.28, 0.338), (0.29, 0.353), (0.3, 0.369), (0.31, 0.381), (0.32, 0.403), (0.33, 0.420), (0.34, 0.439), (0.35, 0.457), (0.36, 0.477), (0.37, 0.497)]

def interpolateur_2000(value, value_ref, iterator = 0):
    if iterator+1 == len(value_ref):
        raise ValueError("La valeur entrée n'est pas dans l'intervalle")
    else:
        if value_ref[iterator][0] <= value <= value_ref[iterator+1][0]:
            return (value_ref[iterator][1] + ((value-value_ref[iterator][0])*((value_ref[iterator+1][1]-value_ref[iterator][1])/(value_ref[iterator+1][0]-value_ref[iterator][0])))) 
        else:
            return interpolateur_2000(value, value_ref, iterator+1)


user_input = float(input("Enter a value between 0.01 and 0.37: "))
result = interpolateur_2000(user_input, value_ref)
print("Voici le oméga:", round(result,3))