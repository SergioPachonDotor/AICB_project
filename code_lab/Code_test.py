class GeneCode:
    def __init__(self) -> None:
        pass
    
# count_of_codon_usage:dict[str, int]  = dict(Counter(filtered_genetic_code.values()))
# claves_a_eliminar:list[str]          = [codon for codon, aa in count_of_codon_usage.items() if aa == 1]
# for clave in claves_a_eliminar:
#     del filtered_genetic_code[clave]