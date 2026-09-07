"""Verify this extracted route-browser companion with the Python standard library."""
import csv
import hashlib
import json
from pathlib import Path


def main():
    root = Path(__file__).resolve().parent
    with (root / 'CHECKSUMS.csv').open(encoding='utf-8', newline='') as handle:
        manifest = list(csv.DictReader(handle))
    for row in manifest:
        file = (root / row['path']).resolve()
        assert root in file.parents, row['path']
        raw = file.read_bytes()
        assert len(raw) == int(row['bytes']), row['path']
        assert hashlib.sha256(raw).hexdigest() == row['sha256'], row['path']
    count = 0
    for dataset, size, methods in [('PaRoutes120',120,9),('RF25',25,4)]:
        data = json.loads((root / f'data/{dataset}.json').read_text(encoding='utf-8'))
        assert len(data['targets']) == size and len(data['methods']) == methods
        assert set(data['routes']) == {t['id'] for t in data['targets']}
        for tid, entries in data['routes'].items():
            assert set(entries) == {m['id'] for m in data['methods']}
            for method, route in entries.items():
                count += 1
                assert route['target_id'] == tid and route['method'] == method
                assert route['status'] != 'missing_source'
                if method == 'reference':
                    assert dataset == 'PaRoutes120' and route['strict_closed'] is None
                    assert route['status'] == 'reference' and route['graph']['reactions']
                    for source in route['sources']:
                        assert hashlib.sha256((root/source['package_path']).read_bytes()).hexdigest() == source['sha256']
                if route['strict_closed']:
                    assert route['graph'] and route['graph']['reactions']
                graphs = [route['graph']]
                if route.get('reviewed_variant'):
                    graphs.append(route['reviewed_variant']['graph'])
                for graph in filter(None, graphs):
                    ids = {n['id'] for n in graph['nodes']}
                    assert len(ids) == len(graph['nodes']) and graph['target_matches_input']
                    assert set(graph['roots']) <= ids
                    for reaction in graph['reactions']:
                        assert reaction['product'] in ids and set(reaction['reactants']) <= ids
    assert count == 1180
    print(f'PASS: {len(manifest)} file checksums; {count} target-method positions; all graph references.')


if __name__ == '__main__':
    main()
