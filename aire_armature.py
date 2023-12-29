from math import pi

armatures = [6, 8, 10, 12, 14, 16, 20, 25, 28, 32, 40]

def choix_armatures(aire_min):
    possibilities = []
    for i in armatures:
        aire_unit = i**2 / 4 * pi * 1e-2
        nb_armatures = aire_min // aire_unit + 1
        if nb_armatures <= 2:
            continue
        possibilities.append((f"{int(nb_armatures)} "+ u"\u00D8" + f" {i}", nb_armatures*aire_unit, aire_unit, nb_armatures))
        possibilities.sort(key=lambda x: x[3]/x[2]+(x[1]-aire_min))
    print(f"Les possibilités d'armatures pour une aire de {aire_min} cm\u00B2 sont :")
    for i in possibilities:
        print(f"\t- {i[0]} \t({i[1]:3.2f} cm\u00B2)")
    print("/!\ N'hésitez pas à vérifier à la main !!!")



choix_armatures(47.46)