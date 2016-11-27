from Bio.Affy import CelFile
with open("Affy/406335477_B.CEL", "rb") as f:
    record = CelFile.read(f, strict=True)
