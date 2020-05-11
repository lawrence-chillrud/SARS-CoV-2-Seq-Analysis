f = 'mult.txt'
nf = '2msa.clustal'
tfile = open(f, 'r')
nfile = open(nf, 'w')
lines = []
for line in tfile.readlines():
    if line[0:8] != 'MT044257' and line[0:8] != 'MT072688':
        lines.append(line)

nfile.writelines(lines)

tfile.close()
nfile.close()
