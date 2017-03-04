#!/fs/sz-user-supported/Linux-x86_64/bin/python

import sys
import math

delta = open(sys.argv[1])
matches = [ ]
contigs = { }
refs = { }

# all the alignments between one ref and one contig
class NucMatch:
    def __init__(self, ref, qry, rlen, qlen):
        self.ref = ref
        self.qry = qry
        self.rlen = int(rlen)
        self.qlen = int(qlen)
        # list of alignments
        self.align = []
        self.type = "u"

# a single alignment
class Align:
    def __init__(self, rstr, rend, qstr, qend, errors):
        self.rstr = int(rstr)
        self.rend = int(rend)
        self.qstr = int(qstr)
        self.qend = int(qend)
        self.errors = int(errors)
        self.reverse = 0

# a single contig containing all the nucmer matchs to it
class Contig:
    def __init__(self, id, len):
        self.id = id
        self.len = len
        self.nucs = []
        self.type = "u"

# a single ref containing all the nucmer matchs to it
class Ref:
    def __init__(self, id, len):
        self.id = id
        self.len = len
        self.nucs = []



#parse delta file
for line in delta:
    line = line.strip()
    fields = line.split()

    if len(fields) == 4:
        nmatch = NucMatch(fields[0], fields[1], fields[2], fields[3])
        matches.append(nmatch)
        contig = contigs.get(nmatch.qry, Contig(nmatch.qry, nmatch.qlen))
        contig.nucs.append(nmatch)
        contigs[nmatch.qry] = contig

        ref = refs.get(nmatch.ref, Contig(nmatch.ref, nmatch.rlen))
        ref.nucs.append(nmatch)
        refs[nmatch.ref] = ref
        
    elif len(fields) == 7:
        a = Align(fields[0], fields[1], fields[2], fields[3], fields[4])
        if a.qstr > a.qend:
            a.qstr, a.qend = a.qend, a.qstr
            a.reverse = 1

        nmatch.align.append(a)



        
# step 1 check for good matches (contig that completing match reference, 5 bp on ends)
# and check for slight mis-assemblies (90%)
good = { }
slight = { }
for id, c in contigs.items():
    for m in c.nucs:

        if c.type != "g":

            for a in m.align:
                # look for good match
                if (a.qstr < 6) and (a.qend > (m.qlen - 6)):
                    good[id] = m
                    m.type = "g"
                    c.type = "g"
                    if slight.has_key(id):
                        del slight[id]
                    break

                # look for slight mis-assembly
                elif math.floor(m.qlen * .9) <= (a.qend - a.qstr):
                    slight[id] = m
                    c.type = "s"
                    m.type = "s"



# Step 2 look for mis-assemblies
# heavy1, heavy2, chimera mis-assemblies as defined in Mihai's validation document
heavy1 = { }
heavy2 = { }

# amount of overlap allowed between to alignments
WINDOW = 25

for id, c in contigs.items():

    if (c.type != "g") and (c.type != "s"):

        for m in c.nucs:
            # check heavy1
            # only one alignment so it is a heavy1
            if (len(m.align) == 1) and (len(c.nucs) == 1):
                c.type = "h1"
                m.type = "h1"
                heavy1[id] = c

            else:
                # check for heavy2
                # need more than one align to single ref to be heavy2
                if len(m.align) > 1:
                    m.align.sort(key=lambda align: align.qstr)
                    cur_str = 0
                    cur_end = 0
                    for a in m.align:
                        if ((a.qstr + WINDOW) > cur_end) and (a.qend > cur_end):
                            c.type = "h2"
                            m.type = "h2"
                            heavy2[id] = c
                            break

                        cur_str = a.qstr
                        cur_end = a.end

# Step 3 check for chimera
chimera = { }
for id, c in contigs.items():
    if (c.type != "g") and (c.type != "s") and (c.type != "h1") and (c.type != "h2"):
        #need multiple alignments to two different refs to be chimeric
        if c.nucs > 1:
            all_aligns = [ ]
            for m in c.nucs:
                all_aligns.extend( m.align)
            all_aligns.sort(key=lambda align: align.qstr)
            cur_str = 0
            cur_end = 0
            for a in all_aligns:
                if (a.qstr + WINDOW) > cur_end:
                    c.type = "c"
                    chimera[id] = c
                    break
                
                cur_str = a.qstr
                cur_end = a.end
            
for id, c in contigs.items():
    if (c.type != "g") and (c.type != "s") and (c.type != "h1") and (c.type != "h2") and (c.type != "c"):
        c.type = "h1"
        heavy1[id] = c
        




# report stats
#print "number of places contigs align " + str(len(matches))

print "description\tTotal number of Contigs\tcontigs <300\t300< contigs <1000\t1000< contigs <10000\t10000< contigs"

c0 = 0
c300 = 0
c1000 = 0
c10000 = 0
for id, c in contigs.items():
    if c.len <= 300:
        c0 += 1
    elif c.len > 300 and c.len <= 1000:
        c300 += 1
    elif c.len > 1000 and c.len <= 10000: 
        c1000 += 1
    elif c.len > 10000:
        c10000 += 1
    else:
        print "***error " + str(c.len)

print "All contigs mapped to reference\t%s\t%d\t%d\t%d\t%d" %  (str(len(contigs)), c0, c300, c1000, c10000)


c0 = 0
c300 = 0
c1000 = 0
c10000 = 0
for id, c in good.items():
    if c.qlen <= 300:
        c0 += 1
    elif c.qlen > 300 and c.qlen <= 1000:
        c300 += 1
    elif c.qlen > 1000 and c.qlen <= 10000: 
        c1000 += 1
    elif c.qlen > 10000:
        c10000 += 1
    else:
        print "***error " + str(c.qlen)

print "Good Contigs\t%s\t%d\t%d\t%d\t%d" %  (str(len(good)), c0, c300, c1000, c10000)

c0 = 0
c300 = 0
c1000 = 0
c10000 = 0
for id, c in slight.items():
    if c.qlen <= 300:
        c0 += 1
    elif c.qlen > 300 and c.qlen <= 1000:
        c300 += 1
    elif c.qlen > 1000 and c.qlen <= 10000: 
        c1000 += 1
    elif c.qlen > 10000:
        c10000 += 1
    else:
        print "***error " + str(c.qlen)

print "Slight Mis-assembly (90 percent identity)\t%s\t%d\t%d\t%d\t%d" %  (str(len(slight)), c0, c300, c1000, c10000)


c0 = 0
c300 = 0
c1000 = 0
c10000 = 0
for id, c in heavy1.items():
    if c.len <= 300:
        c0 += 1
    elif c.len > 300 and c.len <= 1000:
        c300 += 1
    elif c.len > 1000 and c.len <= 10000: 
        c1000 += 1
    elif c.len > 10000:
        c10000 += 1
    else:
        print "***error " + str(c.len)

print "Heavy1 Mis-assembly \t%s\t%d\t%d\t%d\t%d" %  (str(len(heavy1)), c0, c300, c1000, c10000)


c0 = 0
c300 = 0
c1000 = 0
c10000 = 0
for id, c in heavy2.items():
    if c.len <= 300:
        c0 += 1
    elif c.len > 300 and c.len <= 1000:
        c300 += 1
    elif c.len > 1000 and c.len <= 10000: 
        c1000 += 1
    elif c.len > 10000:
        c10000 += 1
    else:
        print "***error " + str(c.len)

print "Heavy2 Mis-assembly \t%s\t%d\t%d\t%d\t%d" %  (str(len(heavy2)), c0, c300, c1000, c10000)

c0 = 0
c300 = 0
c1000 = 0
c10000 = 0
for id, c in chimera.items():
    if c.len <= 300:
        c0 += 1
    elif c.len > 300 and c.len <= 1000:
        c300 += 1
    elif c.len > 1000 and c.len <= 10000: 
        c1000 += 1
    elif c.len > 10000:
        c10000 += 1
    else:
        print "***error " + str(c.len)

print "Chimera Mis-assembly \t%s\t%d\t%d\t%d\t%d" %  (str(len(chimera)), c0, c300, c1000, c10000)



# stats for the ref
# cov is length of all contigs divided by length of reference

print ""
print "Reference id\tReference length\tLength of all aligned contigs\tCoverage"
for id, r in refs.items():
    ctg_total = 0
    cur_m = 0
    for m in r.nucs:
        if m.type == "g" or m.type == "s":
            if m.qlen > cur_m or m.type == "g":
                cur_m = m.qlen


    ctg_total += cur_m

    if ctg_total != 0:
        cov = ctg_total/float(r.len)
        print "%s\t%d\t%d\t%g" % (r.id, r.len, ctg_total, cov)

        
