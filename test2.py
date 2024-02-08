sentence = "Ikar łapał raki"

# Podziel zdanie przez słowo "łapał"
parts = sentence.split("łapał", 1)

# Usuń ewentualne puste znaki na początku i końcu każdej części
parts = [part.strip() for part in parts]


print(parts)