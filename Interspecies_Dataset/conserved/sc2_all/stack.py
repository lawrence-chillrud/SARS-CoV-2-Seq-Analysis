import numpy as np
import os
import glob

FILE = 'msa_play.txt'
SC2 = 'NC_045512.2_21563-25384'
#EXCLUDE = ['NC_019843.3_21456-25517.mfe', 'NC_004718.3_21492-25259.mfe', 'NC_014470.1_21391-25170.mfe']
EXCLUDE = []
LIST_OF_FILES = glob.glob("*.mfe")
SEQS = len(LIST_OF_FILES)
MOD = SEQS + 2
MTX_PATH = "mtx_%s.npy" % (SEQS)

def load_file():
    f = open(FILE, 'r')
    return f.readlines()

def make_dict(lines):
    seqs_dict = {}
    arr = [[] for i in range(SEQS)]
    for i in range(3, len(lines)):
        if (i - 3) % MOD < SEQS:
            cline = lines[i].split()
            if cline[0] not in seqs_dict:
                seqs_dict[cline[0]] = (i-3) % MOD
            arr[(i-3) % MOD] += list(cline[1])

    for key in seqs_dict.keys():
        idx = seqs_dict[key]
        sequence = ""
        for char in arr[idx]:
            sequence += char
        seqs_dict[key] = sequence.replace('T', 'U')

    np.save(MTX_PATH, arr)
    return np.array(arr), seqs_dict

def make_dist(mtx):
    dist = {}
    mtx = np.transpose(mtx)
    for i in range(mtx.shape[0]):
        key = ""
        for char in mtx[i]:
            key += char
        if key in dist:
            dist[key] += 1
        else:
            dist[key] = 1

    return dist

def insert_gaps(seqs_dict):
    list_of_files = LIST_OF_FILES
    dbs_w_gaps = {}

    for f in list_of_files:
        mfe_file = open(f, 'r')
        lines = mfe_file.readlines()
        dot_bracket = lines[len(lines)-1].split()[0]
        seq_w_gaps = seqs_dict[f[0:len(f)-4]]
        db_w_gaps = ""
        db_counter = 0

        for i in range(len(seq_w_gaps)):
            if seq_w_gaps[i] != "-":
                db_w_gaps += dot_bracket[db_counter]
                db_counter += 1
            else:
                db_w_gaps += "-"
        dbs_w_gaps[f[0:len(f)-4]] = db_w_gaps

    return dbs_w_gaps

'''
MN996532.1_21545-25354   RatG13
NC_004718.3_21492-25259  SARS
NC_014470.1_21391-25170  Bat BM48-31
NC_019843.3_21456-25517  MERS
NC_045512.2_21563-25384  SC2
'''
def make_db_array(dbs_w_gaps):
    l = [[] for i in range(SEQS)]
    
    for i, species in enumerate(dbs_w_gaps.keys()):
        print(species)
        l[i] = list(dbs_w_gaps[species])
        #print('\n\n')
        #print(l[i])

    arr = np.array(l)
    return arr

'''
MN996532.1_21545-25354   RatG13
NC_004718.3_21492-25259  SARS
NC_014470.1_21391-25170  Bat BM48-31
NC_019843.3_21456-25517  MERS
NC_045512.2_21563-25384  SC2
'''
def compare2(seqs_dict, dbs_w_gaps, threshold=20):
    input_mtx = make_db_array(dbs_w_gaps)
    mtx = np.transpose(input_mtx)
    streaks = []
    cstreak = 0
    sites = []
    start_site = 0
    end_site = 0
    for i in range(mtx.shape[0]):
        site = mtx[i]
        if site[0] == site[SEQS-1]:
            if cstreak == 0:
                start_site = i
            cstreak += 1
        elif cstreak != 0:
            if cstreak >= threshold:
                streaks.append(cstreak)
                end_site = i
                sites.append((start_site, end_site))
            cstreak = 0
    print(streaks)
    print(sites)
    print("max: ", max(streaks))
    print("total more than %d: %d" % (threshold, len(streaks)))
    print_conserved(mtx, dbs_w_gaps, seqs_dict, sites, [0, SEQS-1])

def print_conserved(db, dbs_w_gaps, seqs_dict, sites, species_list):
    keys = list(dbs_w_gaps.keys())
    for j in species_list:
        f = open('%s_cons.txt' % (keys[j]), 'w')
        for loc in sites:
            f.write('> ' + str(loc) + ' ' + keys[j] + '\n')
            line1 = seqs_dict[keys[j]][loc[0]:loc[1]]
            line2 = dbs_w_gaps[keys[j]][loc[0]:loc[1]]
            #line2 = ''
            #for i in range(loc[0], loc[1]):
                #line2 += db[i][j]
            f.write(line1 + '\n')
            f.write(line2 + '\n')
        f.close()

def count_all_same(db_dist):
    counts = np.zeros(SEQS)
    for entry in db_dist.keys():
        count = db_dist[entry]
        counts[len(set(entry))-1] += count
    
    print("[all same, 4 same, 3 same, 2 same, all diff]")
    print(counts)

def translate(entry, sc2_idx):
    translation = ['x', 'q', 'w', 'y', 'z']
    entrylist = list(entry)
    x = entrylist[sc2_idx]
    new_entry = []
    seen = [0]
    for i in range(SEQS):
        if entrylist[i] == x:
            new_entry.append(0)
        elif entrylist[i] in seen:
            new_entry.append(seen.index(entrylist[i]))
        else:
            new_entry.append(len(seen))
            seen.append(entrylist[i])
    ne = ''
    for j in range(SEQS):    
        ne += translation[new_entry[j]]

    return ne

def collapse(distribution, sc2_idx):
    new_d = {}
    total = 0
    for entry in distribution:
        ne = translate(entry, sc2_idx)  
        if ne in new_d:
            new_d[ne] += distribution[entry]
            total += distribution[entry]
        else:
            new_d[ne] = distribution[entry]
            total += distribution[entry]

    return new_d

def standardize(seqs_dict, dbs_w_gaps, threshold=0):
    files = LIST_OF_FILES
    seqs_mtx = [[] for i in range(SEQS)]
    dbs_mtx = [[] for i in range(SEQS)]
    keys = []
    for i, f in enumerate(files):
        key = f[0:len(f)-4]
        keys.append(key)
        seqs_mtx[i] = list(seqs_dict[key])
        dbs_mtx[i] = list(dbs_w_gaps[key])

    seqs = np.array(seqs_mtx)
    dbs = np.array(dbs_mtx)
    dist = make_dist(seqs)
    sites = np.transpose(seqs)
    # now onto the stack:
    covid = keys.index(SC2)
    # FLAG
    collapsed_dist = collapse(dist, covid)
    #collapsed_dist = dist
    covid_seq = seqs[covid]
    covid_db = dbs[covid]
    #print(dbs_w_gaps[SC2][617:691])
    #print('\n')
    covid_db_nogaps = dbs_w_gaps[SC2].replace('-', '')
    #print(covid_db_nogaps[471:516])
    #print(covid_db_nogaps[517:666])
    #print(covid_db_nogaps[666:701])
    dim = seqs.shape[1]
    stack = []
    streaks = []
    locs = []
    cstreak = 0
    start = 0
    end = 0
    i = 0
    while i < dim:
        dot_brac = covid_db[i]
        site = sites[i]
        entry = ''
        for char in site:
            entry += char
        te = translate(entry, covid)
        if dot_brac == '(':
            stack.append((te, i))
        elif dot_brac == ')':
            #print('stack on loc %d:' % (i))
            #print(stack)
            check = stack.pop()
            if te == check[0]:
                if covid_db[i+1] != ')':
                    end = i
                    cstreak = end - check[1]
                    if cstreak >= threshold:
                        streaks.append(cstreak)
                        locs.append((check[1], end))
                        #cstreak = 0
                    #start = i
            else:
                prev_stack = check[1]
                j = 0
                while j < len(stack):
                    if stack[-1][1] - prev_stack == 1:
                        prev_stack = stack[-1][1]
                        stack.pop()
                        i += 1
                    else:
                        j = len(stack)
        i += 1
    #cons(locs, seqs_dict[SC2], dbs_w_gaps[SC2])
    return streaks, locs, collapsed_dist, sites, covid, seqs_dict[SC2].replace('-', '')

def cons(locs, seq, db):
    f = open('covid_cons_no_gaps.txt', 'w')
    for loc in locs:
        f.write('> ' + str(loc[0]) + '-' + str(loc[1]) + ' log likelihood: ' + str(loc[2]) + 'random seq of same size log likelihood: ' + str(loc[3]) + '\n')
        line1 = seq[loc[0]:loc[1]+1]
        line2 = db[loc[0]:loc[1]+1]
        f.write(line1 + '\n')
        f.write(line2 + '\n\n')
    f.close()

def find(seqs, locs):
    no_gaps = seqs[SC2].replace('-', '')
    seq = seqs[SC2]
    ans = []
    for loc in locs:
        find = seq[loc[0]:loc[1]+1].replace('-', '')
        start = no_gaps.index(find)
        end = start + len(find)
        ans.append((start, end))
    
    starts = []
    ends = []
    for loc in ans:
        starts.append(loc[0])
        ends.append(loc[1])
    fans = []
    for i in range(len(ans)):
        cloc = ans[i]
        if cloc[0] - 1 in starts and cloc[1] + 1 in ends and starts.index(cloc[0] - 1) == ends.index(cloc[1] + 1):
            pass
        else:
            fans.append(cloc)

    return fans

def getnum(blist):
    current_num = 0
    for i in range(len(blist)):
        if blist[i] == ')' and blist[i-1] != ')':
            current_num += 1
    return current_num

def nest(newlocs, cov_db_nogaps):
    #print(cov_db_nogaps[7:3220])
    final_locs = []
    current_stretch = newlocs.pop()
    blist = list(cov_db_nogaps[current_stretch[0]:current_stretch[1]+1])
    current_num = getnum(blist)
    track = 1
    i = 0
    size = len(newlocs)
    while i < size:
        check = newlocs.pop()
        blist = list(cov_db_nogaps[check[0]: check[1]+1])
        #print(current_stretch[0])
        #print(check[0])
        if current_stretch[0] < check[0]:
            if current_num - getnum(blist) == track:
                track += 1
            else:
                current_stretch = check
                current_num = getnum(blist)
                track = 1
            #print("nested!")
            #current_num -= 1
        else:
            #print(current_num)
            if current_num - track == 0:
                final_locs.append(current_stretch)
            current_stretch = check
            current_num = getnum(blist)
            track = 1
        if i == size - 1 and current_num - track == 0:
            final_locs.append(current_stretch)
        i += 1
    #print(current_num)
    return final_locs

def tack_on(loc, cov_db_nogaps, sites):
    start = loc[0]
    end = loc[1]
    before = start - 1
    entry = ""
    for char in sites[before]:
        entry += char
    e = translate(entry, covid)
    addone = False
    while before >= 0 and cov_db_nogaps[before] == '.':
        centry = ""
        for char in sites[before]:
            centry += char
        ce = translate(centry, covid)
        if ce == e and cov_db_nogaps[before-1] == '.':
            before -= 1
        else:
            break
    if before == start - 1:
        before += 1
    
    after = end
    #print(cov_db_nogaps[after])
    entry = ""
    for char in sites[after]:
        entry += char
    e = translate(entry, covid)
    while after < sites.shape[0] and cov_db_nogaps[after] == '.':
        centry = ""
        for char in sites[after]:
            centry += char
        ce = translate(centry, covid)
        if ce == e:
            after += 1
        else:
            #after -= 1
            break
    
    laced_loc = (before, after-1)

    return laced_loc

def lace(flocs, cov_db_nogaps, sites, cdist, covid):
    intsarr = []
    for i in range(sites.shape[0]):
        if sites[i][covid] == '-':
             intsarr.append(i)
    sites = np.delete(sites, intsarr, axis=0)
    final_dist = make_dist(np.transpose(sites))
    tot = 0
    final_dist = collapse(final_dist, covid)
    for item in final_dist:
        tot += final_dist[item]

    laced_locs = []
    for i in range(len(flocs)):
        loc = flocs.pop()
        laced_locs.append(tack_on(loc, cov_db_nogaps, sites))
    
    fin = []
    for i in range(len(laced_locs)-1):
        c = laced_locs[i]
        cn = laced_locs[i+1]
        if cn[0] - c[1] <= 1:
            fin.append((c[0], cn[1]))
        else:
            if len(fin) == 0:
                fin.append(c)
            elif fin[-1][1] != c[1]:
                fin.append(c)
    if fin[-1][1] != laced_locs[-1][1]:
        fin.append(laced_locs[-1])
    fin2 = []
    for i in range(len(fin)-1):
        c = fin[i]
        cn = fin[i+1]
        if cn[0] <= c[1]:
            fin2.append((c[0], cn[1]))
        else:
            if len(fin2) == 0:
                fin2.append(c)
            elif fin2[-1][1] != c[1]:
                fin2.append(c)
    if fin2[-1][1] != fin[-1][1]:
        fin2.append(fin[-1])
    fin3 = []
    for i in range(len(fin2)-1):
        c = fin2[i]
        cn = fin2[i+1]
        if cn[0] <= c[1]:
            fin3.append((c[0], cn[1]))
        else:
            if len(fin3) == 0:
                fin3.append(c)
            elif fin3[-1][1] != c[1]:
                fin3.append(c)
    if fin3[-1][1] != fin2[-1][1]:
        fin3.append(fin2[-1])
    return fin3, final_dist, tot, sites


def probs(laced_flocs, fd, tot, final_sites, covid):
    s_gene_probs = []
    str_sites = []
    for site in final_sites:
        entry = ""
        for char in site:
            entry += char
        s_gene_probs.append(fd[translate(entry, covid)] / tot)
        str_sites.append(translate(entry, covid))
    s_gene_probs = np.log(s_gene_probs)
    sigma = np.var(s_gene_probs)
    mu = np.mean(s_gene_probs) 
    # np.random.normal(mu, sigma, length of seq)
    fprobs = []
    samplesarray = []
    for loc in laced_flocs:
        samples = (loc[1]+1) - loc[0]
        samplesarray.append(samples)
        prob = 0
        for i in range(loc[0], loc[1]+1):
            prob += s_gene_probs[i]
        draw = np.random.normal(mu, sigma, samples)
        fprobs.append((loc[0], loc[1], prob, np.sum(draw)))
    print(max(samplesarray))
    return fprobs

if __name__ == '__main__':
    lines = load_file()
    _, seqs_dict = make_dict(lines)
    for item in EXCLUDE:
        LIST_OF_FILES.remove(item)
        del seqs_dict[item[0:len(item)-4]]
    SEQS = len(LIST_OF_FILES)
    dbs_w_gaps = insert_gaps(seqs_dict)
    streaks, locs, cdist, sites, covid, sd = standardize(seqs_dict, dbs_w_gaps)
    #print("num of locs from stand: ", len(locs))
    newlocs = find(seqs_dict, locs)
    #print("num of newlocs from find: ", len(newlocs))
    cov_db_nogaps = dbs_w_gaps[SC2].replace('-', '')
    #print("NEW LOCS: ", newlocs)
    #print("look at new locs:\n")
    '''
    for loc in newlocs:
        print(str(loc))
        print(cov_db_nogaps[loc[0]:loc[1]])
        print('\n')
    '''
    flocs = nest(newlocs, cov_db_nogaps)
    #print("FINAL LOCS: ", flocs)
    #print("look at final locs:\n")
    '''
    for loc in flocs:
        print(str(loc))
        print(cov_db_nogaps[loc[0]:loc[1]])
        print('\n')
    '''
    laced_flocs, fd, tot, final_sites = lace(flocs, cov_db_nogaps, sites, cdist, covid)
    attached_probs = probs(laced_flocs, fd, tot, final_sites, covid)
    #print('-' * 80)
    #print("LACED FINAL LOCS: ", laced_flocs)
    #print("557-581")
    #print(cov_db_nogaps[557:582])
    cons(attached_probs, sd, cov_db_nogaps)
    print(len(attached_probs))

    #for loc in attached_probs:
        #print(str(loc))
        #print(cov_db_nogaps[loc[0]:loc[1]+1])
        #print("This sequence had probability: %f. A random one of similar size had prob: %f" % (loc[2], loc[3]))
        #print('\n')
    #nest(newlocs)
    # when SEQS=5:
    # mtx[0] = NC_019843, mtx[1] = MN996532, mtx[2] = NC_045512, mtx[3] = NC_004718, mtx[4] = NC_014470
    # 0 = MERS, 1 = RaTG13, 2 = SC2, 3 = SARS, 4 = Bat BM48-31
    #compare2(seqs_dict, dbs_w_gaps)
    #print(db_mtx)
    #print(np.transpose(db_mtx))
    #db_dist = make_dist(db_mtx)
    #count_all_same(db_dist)
    #make_dist(mtx)
