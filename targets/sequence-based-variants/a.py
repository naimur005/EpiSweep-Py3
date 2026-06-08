import sys
f = open("lst_cwbd_MSA.fasta").readlines()

for n, i in enumerate(f):
    if i.count(sys.argv[1]):
        break

print(" ".join(list(f[n+1][:-1])))
