from collections import Counter

def calcular_ENC(frecuencias_codones):
    
    total_codones = sum(frecuencias_codones.values())

    sumatoria = sum((fi / total_codones) ** 2 for fi in frecuencias_codones.values())
    print(type(sumatoria))

    ENC = 2 / (1 + sumatoria)
    print(type(ENC))
    
    return ENC

# Ejemplo de uso
frecuencias_codones_ejemplo = {'TTT': 10, 'TTC': 15, 'TTA': 5, 'TTG': 8}
enc_resultado = calcular_ENC(frecuencias_codones_ejemplo)
print(f"El valor de ENC es: {enc_resultado}")
