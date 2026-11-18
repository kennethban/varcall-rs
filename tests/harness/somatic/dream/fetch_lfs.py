#!/usr/bin/env python3
"""fetch_lfs.py OUT_DIR POINTER...: download Git LFS objects of the GATK repository for the given pointer files."""
import json, sys, os, hashlib, urllib.request
out = sys.argv[1]
objs = []
for p in sys.argv[2:]:
    t = open(p).read()
    if not t.startswith('version https://git-lfs'):
        # not a pointer: copy as is
        open(os.path.join(out, os.path.basename(p)), 'w').write(t); continue
    oid = [l.split('sha256:')[1] for l in t.splitlines() if l.startswith('oid')][0]
    size = int([l.split()[1] for l in t.splitlines() if l.startswith('size')][0])
    objs.append((os.path.basename(p), oid, size))
req = urllib.request.Request('https://github.com/broadinstitute/gatk.git/info/lfs/objects/batch',
    data=json.dumps({'operation': 'download', 'transfers': ['basic'], 'objects': [{'oid': o, 'size': s} for _, o, s in objs]}).encode(),
    headers={'Accept': 'application/vnd.git-lfs+json', 'Content-Type': 'application/vnd.git-lfs+json'})
resp = json.load(urllib.request.urlopen(req))
href = {o['oid']: o['actions']['download']['href'] for o in resp['objects']}
for name, oid, size in objs:
    data = urllib.request.urlopen(href[oid]).read()
    assert len(data) == size and hashlib.sha256(data).hexdigest() == oid, name
    open(os.path.join(out, name), 'wb').write(data)
    print(f'{name} {size} ok')
