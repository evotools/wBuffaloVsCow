import sys, os


def main():
    infld = sys.argv[1]
    minS = int(sys.argv[2])
    naut = int(sys.argv[3])
    savedBreeds = []
    for f in os.listdir(infld):
        nrows = len([l for l in open(os.path.join(infld, f))])
        if nrows >= minS:
            savedBreeds.append(f.strip(".txt"))
    sys.stderr.write("Considered {} breeds.\n".format(len(savedBreeds)) )

    pairs = [[savedBreeds[b1], savedBreeds[b2]] for b1 in range(0, len(savedBreeds)-1) for b2 in range(b1, len(savedBreeds)) if b1 != b2]
    sys.stderr.write("Created {} pairwise comparisons.\n".format(len(pairs)) )
    for (breed1,breed2) in pairs:
        for nchr in range(0, naut):
            print("{}\t{}\t{}".format(nchr + 1, breed1, breed2))

if __name__ == "__main__":
    main()
