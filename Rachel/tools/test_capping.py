#!/usr/bin/env python
"""smart_cap / custom_cap 全面测试。"""
import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
sys.stdout.reconfigure(encoding='utf-8')
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from Rachel.chem_tools.smart_cap import custom_cap, suggest_capping
from Rachel.chem_tools._rdkit_utils import parse_mol, canonical

fails = []
total = 0

def test(name, smi, bond, ci, cj, expect_ok=True, expect_ring=None,
         expect_frag_count=None, expect_contains=None):
    global total, fails
    total += 1
    r = custom_cap(smi, bond, ci, cj, name)
    ok = r.get('ok', False)
    frags = r.get('fragments', [])
    ring = r.get('ring_opening', False)
    err = r.get('error', '')

    issues = []
    if ok != expect_ok:
        issues.append(f'ok={ok} expected {expect_ok}')
    if ok:
        if expect_ring is not None and ring != expect_ring:
            issues.append(f'ring={ring} expected {expect_ring}')
        if expect_frag_count is not None and len(frags) != expect_frag_count:
            issues.append(f'frag_count={len(frags)} expected {expect_frag_count}')
        for f in frags:
            if parse_mol(f) is None:
                issues.append(f'invalid SMILES: {f}')
        for f in frags:
            m = parse_mol(f)
            if m:
                for a in m.GetAtoms():
                    if a.GetAtomicNum() == 0:
                        issues.append(f'dummy atom in fragment: {f}')
                        break
        if expect_contains:
            for substr in expect_contains:
                found = any(substr in f for f in frags)
                if not found:
                    issues.append(f'expected substr "{substr}" not in frags')

    status = 'PASS' if not issues else 'FAIL'
    if issues:
        fails.append((name, issues, frags, err))
    print(f'{status:4s} | {name}' + (f' | {issues}' if issues else ''))

print("=" * 60)
print("Section 1: Basic non-ring bonds")
print("=" * 60)
test('amide C(=O)-N → COOH + amine', 'CC(=O)Nc1ccccc1', (1,3), 'O', '[H]',
     expect_ring=False, expect_frag_count=2)
test('amide → acid chloride + amine', 'CC(=O)Nc1ccccc1', (1,3), 'Cl', '[H]',
     expect_ring=False, expect_frag_count=2)
test('ester C(=O)-O → COOH + ROH', 'CC(=O)OC', (1,3), 'O', '[H]',
     expect_ring=False, expect_frag_count=2)
test('ether Ar-O-Me → ArOH + MeBr', 'COc1ccccc1', (1,0), '[H]', 'Br',
     expect_ring=False, expect_frag_count=2)
test('N-alkyl → amine + RBr', 'CNCc1ccccc1', (1,0), '[H]', 'Br',
     expect_ring=False, expect_frag_count=2)

print()
print("=" * 60)
print("Section 2: Cross-coupling caps (multi-atom)")
print("=" * 60)
test('Suzuki Br + B(O)O', 'c1ccc(-c2ccccc2)cc1', (3,4), 'Br', 'B(O)O',
     expect_ring=False, expect_frag_count=2)
test('Suzuki reversed B(O)O + Br', 'c1ccc(-c2ccccc2)cc1', (3,4), 'B(O)O', 'Br',
     expect_ring=False, expect_frag_count=2)
test('Negishi Br + [Zn]Cl', 'c1ccc(-c2ccccc2)cc1', (3,4), 'Br', '[Zn]Cl',
     expect_ring=False, expect_frag_count=2)
test('Stille Br + [Sn](C)(C)C', 'c1ccc(-c2ccccc2)cc1', (3,4), 'Br', '[Sn](C)(C)C',
     expect_ring=False, expect_frag_count=2)
test('Kumada Br + [Mg]Br', 'c1ccc(-c2ccccc2)cc1', (3,4), 'Br', '[Mg]Br',
     expect_ring=False, expect_frag_count=2)
test('Hiyama Br + [Si](C)(C)C', 'c1ccc(-c2ccccc2)cc1', (3,4), 'Br', '[Si](C)(C)C',
     expect_ring=False, expect_frag_count=2)

print()
print("=" * 60)
print("Section 3: Special cap types ([H], =O)")
print("=" * 60)
test('[H] + [H] simple', 'CCCC', (1,2), '[H]', '[H]',
     expect_ring=False, expect_frag_count=2)
test('=O reductive amination', 'CNCC', (1,2), '[H]', '=O',
     expect_ring=False, expect_frag_count=2)
test('F cap (SNAr)', 'COc1ccccc1', (1,0), '[H]', 'F',
     expect_ring=False, expect_frag_count=2)
test('I cap', 'CCc1ccccc1', (1,2), 'I', '[H]',
     expect_ring=False, expect_frag_count=2)

print()
print("=" * 60)
print("Section 4: Ring-opening (in-ring bonds)")
print("=" * 60)
test('cyclohexane [H]+[H]', 'C1CCCCC1', (0,1), '[H]', '[H]',
     expect_ring=True, expect_frag_count=1)
test('cyclohexane Br+Cl', 'C1CCCCC1', (0,1), 'Br', 'Cl',
     expect_ring=True, expect_frag_count=1)
test('cyclohexanone ring =O+[H]', 'O=C1CCCCC1', (1,6), '[H]', 'O',
     expect_ring=True, expect_frag_count=1)
test('piperidine ring N-C', 'C1CCNCC1', (3,2), '[H]', '=O',
     expect_ring=True, expect_frag_count=1)
test('FC acylation retro (complex)', 'COc1ccc2c(c1)[C@@]1(C)CC(=O)CC(C)(C)[C@@H]1CC2=O',
     (5,19), '[H]', 'O', expect_ring=True, expect_frag_count=1)
test('naphthalene ring Ar-Ar', 'c1ccc2ccccc2c1', (4,5), '[H]', '[H]',
     expect_ring=True, expect_frag_count=1)
test('lactone ring C(=O)-O', 'O=C1CCCCO1', (1,6), 'O', '[H]',
     expect_ring=True, expect_frag_count=1)
test('lactam ring C(=O)-N', 'O=C1CCCCN1', (1,6), 'O', '[H]',
     expect_ring=True, expect_frag_count=1)

print()
print("=" * 60)
print("Section 5: Stereochemistry preservation")
print("=" * 60)
test('chiral center preserved', 'C[C@@H](O)c1ccccc1', (0,1), 'Br', '[H]',
     expect_ring=False, expect_frag_count=2)
test('E/Z alkene adjacent', 'C/C=C/Cc1ccccc1', (3,4), 'Br', '[H]',
     expect_ring=False, expect_frag_count=2)

print()
print("=" * 60)
print("Section 6: Heterocycles and complex scaffolds")
print("=" * 60)
test('pyridine Ar-C', 'Cc1ccncc1', (0,1), 'Br', '[H]',
     expect_ring=False, expect_frag_count=2)
test('indole ring bond', 'c1ccc2[nH]ccc2c1', (3,4), '[H]', '[H]',
     expect_ring=True, expect_frag_count=1)
test('furan ring bond', 'c1ccoc1', (0,4), '[H]', '[H]',
     expect_ring=True, expect_frag_count=1)
test('thiophene Ar-C', 'Cc1ccsc1', (0,1), 'Br', '[H]',
     expect_ring=False, expect_frag_count=2)
test('steroid-like ring bond', 'O=C1CC[C@@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C',
     (1,2), '[H]', 'O', expect_ring=True, expect_frag_count=1)
test('bridged bicyclic', 'C1CC2CCC1C2', (0,5), '[H]', '[H]',
     expect_ring=True, expect_frag_count=1)
test('spiro ring bond', 'C1CCC2(CC1)CCCCC2', (0,5), '[H]', '[H]',
     expect_ring=True, expect_frag_count=1)

print()
print("=" * 60)
print("Section 7: Charged species")
print("=" * 60)
test('quaternary N', 'CC[N+](C)(C)C', (1,2), '[H]', 'Br',
     expect_ring=False, expect_frag_count=2)
test('carboxylate', 'CC(=O)[O-]', (0,1), '[H]', 'Cl',
     expect_ring=False, expect_frag_count=2)

print()
print("=" * 60)
print("Section 8: Edge cases and error handling")
print("=" * 60)
test('invalid SMILES', 'NOT_A_SMILES', (0,1), 'Br', '[H]', expect_ok=False)
test('no bond between atoms', 'CCCC', (0,3), 'Br', '[H]', expect_ok=False)
test('invalid cap SMILES', 'CCCC', (1,2), 'XXXXX', '[H]', expect_ok=False)
test('atom out of range', 'CCCC', (0,99), 'Br', '[H]', expect_ok=False)
test('negative atom idx', 'CCCC', (-1,2), 'Br', '[H]', expect_ok=False)
test('same atom', 'CCCC', (1,1), 'Br', '[H]', expect_ok=False)

print()
print("=" * 60)
print("Section 9: suggest_capping regression")
print("=" * 60)
# Verify suggest_capping still works for common patterns
r = suggest_capping('c1ccc(-c2ccccc2)cc1', (3,4))
sc_ok = r.get('ok') and r.get('n_proposals', 0) >= 2
print(f'{"PASS" if sc_ok else "FAIL"} | suggest_capping biphenyl Ar-Ar (n={r.get("n_proposals",0)})')
if not sc_ok: fails.append(('suggest biphenyl', ['n_proposals < 2'], [], ''))
total += 1

r2 = suggest_capping('CC(=O)Nc1ccccc1', (1,3))
sc2_ok = r2.get('ok') and r2.get('n_proposals', 0) >= 1
print(f'{"PASS" if sc2_ok else "FAIL"} | suggest_capping amide C-N (n={r2.get("n_proposals",0)})')
if not sc2_ok: fails.append(('suggest amide', ['n_proposals < 1'], [], ''))
total += 1

r3 = suggest_capping('COc1ccc2c(c1)[C@@]1(C)CC(=O)CC(C)(C)[C@@H]1CC2=O', (5,19))
sc3_ok = r3.get('ok') and r3.get('n_proposals', 0) >= 1
for p in r3.get('proposals', []):
    if p.get('ring_opening') is not True:
        sc3_ok = False
print(f'{"PASS" if sc3_ok else "FAIL"} | suggest_capping ring bond (n={r3.get("n_proposals",0)}, ring_opening)')
if not sc3_ok: fails.append(('suggest ring', ['ring_opening not set'], [], ''))
total += 1

# ═══ Summary ═══
print()
print("=" * 60)
passed = total - len(fails)
print(f'RESULT: {passed}/{total} passed')
if fails:
    print(f'\nFailed ({len(fails)}):')
    for name, issues, frags, err in fails:
        print(f'  {name}: {issues}')
        if frags: print(f'    frags={frags}')
        if err: print(f'    err={err}')
else:
    print('All tests passed!')
