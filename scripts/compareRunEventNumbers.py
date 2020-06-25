import re
import json

cern_files = [
    ''
]

tallinn_files = [
    ''
]

json_files = [
    ''
]

# Good run ranges
good_run_ranges = {}
for json_file in json_files:
    good_run = open(json_file, 'r')
    good_run_dict = json.loads(good_run.read())
    for key, val in good_run_dict.iteritems():
        good_run_ranges[int(key)] = val
    good_run.close()

def event_in_json(event_info):
    run, lumi, evt = event_info
    if run in good_run_ranges:
        lumiOK = False
        for lumi_range in good_run_ranges[run]:
            if lumi >= lumi_range[0] and lumi <= lumi_range[1]:
                lumiOK = True
                break
        return lumiOK
    else:
        return False

llr_matcher = re.compile('(?P<run>\d*):(?P<lumi>\d*):(?P<evt>\d*)')
#tifr_matcher = re.compile('\*(?P<garbage>[0-9 ]+)\*(?P<run>[0-9 ]+)\*(?P<lumi>[0-9 ]+)\*(?P<evt>[0-9 ]+)\*')
tifr_matcher = llr_matcher
#llr_matcher = tifr_matcher 

llr_events = set([])
tifr_events = set([])

def fill(files, matcher, event_set):
    for filename in files:
        file = open(filename, 'r')
        for line in file.readlines():
            match = matcher.match(line)
            if match:
                event_set.add(tuple(
                    map(int, map(match.group, ['run', 'lumi', 'evt']))
                ))

fill(llr_files, llr_matcher, llr_events)
fill(tifr_files, tifr_matcher, tifr_events)

print "Found %i LLR events" % len(llr_events)
print "Found %i TIFR events" % len(tifr_events)

print "There are %i common events" % len(llr_events.intersection(tifr_events))

llr_only = llr_events - tifr_events
print "======= LLR only events (%i) ========" % len(llr_only)
for event in llr_only:
    print ":".join(map(str, event))

tifr_only = tifr_events - llr_events
print "======= TIFR only events (%i) ========" % len(tifr_only)
for event in tifr_only:
    print ":".join(map(str, event))

tifr_output_file = open('tifr_events.txt', 'w')
for event in tifr_events:
    tifr_output_file.write(":".join(map(str, event)) + "\n")

llr_only_file = open('llr_only_events.txt', 'w')
for event in llr_only:
    llr_only_file.write(":".join(map(str, event)) + "\n")

tifr_only_file = open('tifr_only_events.txt', 'w')
for event in tifr_only:
    tifr_only_file.write(":".join(map(str, event)) + "\n")

print "There are %i LLR  events that fail JSON cut" % len([
    event for event in llr_events if not event_in_json(event)])
print "There are %i TIFR events that fail JSON cut" % len([
    event for event in tifr_events if not event_in_json(event)])

