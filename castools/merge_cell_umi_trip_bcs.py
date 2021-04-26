import sys

cellumitrip_count = {}
for f1 in sys.argv[1:]:
    with open(f1) as cell_umi_trip_count_f1:
        for line in cell_umi_trip_count_f1:
            cell, umi, trip, count = line.rstrip("\n").split(" ")
            count = int(count)
            uid = cell + "\t" + umi + "\t" + trip
            if uid not in cellumitrip_count:
                cellumitrip_count[uid] = 0
            cellumitrip_count[uid] += count

for uid in cellumitrip_count:
    print(uid + "\t" + str(cellumitrip_count[uid]))
