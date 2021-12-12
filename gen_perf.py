
vertices = [5, 10, 100, 1000, 10000]
weights = [1, 100]

for v in vertices:
    e = round(v**(3/2))
    for w in weights:
        print(f"{v} {e} {w} scratch/data-{v}-{e}-{w}.txt")
