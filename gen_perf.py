
vertices = [5, 10, 100]
weights = [1, 100]

for v in vertices:
    e = round(v**(3/2))
    for w in weights:
        print(f"{v} {e} {w} perf/data-{v}-{e}-{w}.txt")
