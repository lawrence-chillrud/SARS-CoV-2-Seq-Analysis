f = 'testrun.out'
nf = 'new.out'
tfile = open(f, 'r')
nfile = open(nf, 'w')
lines = []
for i, line in enumerate(tfile.readlines()):
    #print("Line #%d: %s" % (i, line))
    if line[0] == '>':
        new_line = ''
        found = False
        for char in line:
            if char == '.':
                found = True
            if not found:
                new_line += char
        #nfile.write(new_line)
        #nfile.write('\n')
        lines.append(new_line+'\n')
    else:
        lines.append(line)

nfile.writelines(lines)

tfile.close()
nfile.close()
