from chemicalconverters import NamesConverter

converter = NamesConverter(model_name="knowledgator/IUPAC2SMILES-canonical-base")
print(converter.iupac_to_smiles('hexane'))
