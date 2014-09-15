lines = open('patran.rpt', 'r').readlines()[14:]

nheader = 2
headers = ''
for n in xrange(nheader):
    header = lines[n].strip('\n\r\t')
    headers += header
headers = headers.split('-')

headers2 = []
for header in headers:
    if header:
        headers2.append(header)
headers = headers2
#headers2 = [header if header.strip() for header in headers]
print headers

i = 2
j = 0
res = 0
results = {0:[], 1:[]}
while i < len(lines):
    data = []

    for n in xrange(nheader):
        data += lines[i+n].strip().split()
    if "MSC.Patran" in data:
        res += 1
        i += 11
        continue
    if len(data) == 0:
        break
    results[res].append(data)
    #print "i=%s j=%s data =%s" %(i, j, data)
    i += nheader
    j += 1
    if int(float(data[0])) == 0:
        asdf

key_map = {}
for ikey, key in enumerate(headers):
    key_map[key] = ikey

for key, rows in results.iteritems():
    data2 = {}
    iz = key_map['Z Location']
    iozz = key_map['Z Component']
    for row in rows:
        z = float(row[iz])
        ozz = float(row[iozz])
        print "z=%s ozz=%s" % (z, ozz)
    break
