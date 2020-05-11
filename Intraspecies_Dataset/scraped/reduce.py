
f = open('percent_identity_mtx.pim.txt', 'r')
care_s = open('care-S.gb', 'w')
lines = f.readlines()
seqs = []
sim = []
for i in range(6, len(lines)):
    arr = lines[i].split()
    seqs.append(arr[1])
    count = 0
    for i in range(2, len(arr)):
        if arr[i] == "100.00":
            count += 1
    sim.append(count)
seen = []
care = []
i = 0
care.append('NC_045512')
while i < len(seqs):
    if sim[i] != 220:
        care.append(seqs[i])

    i += sim[i]

for fi in care:
    temp = open(fi + '.txt', 'r')
    temp_lines = temp.readlines()
    care_s.writelines(temp_lines)
    temp.close()
f.close()
care_s.close()
