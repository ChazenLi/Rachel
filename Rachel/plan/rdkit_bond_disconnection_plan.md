# RDKit Bond Disconnection & Molecule Map Validation Plan

## Overview

This document provides a comprehensive guide for using RDKit to disconnect bonds freely and validate molecule maps correctly, suitable for LLM/agent molecular manipulation workflows.

---

## Part 1: Bond Disconnection Methods

### 1.1 RWMol (Read-Write Molecule) - Primary Method

```python
from rdkit import Chem

# Create RWMol from existing molecule
mol = Chem.MolFromSmiles('CCOCC')
rwmol = Chem.RWMol(mol)

# Remove a bond between two atoms by indices
rwmol.RemoveBond(1, 2)

# Get the modified molecule
result = rwmol.GetMol()
Chem.SanitizeMol(result)
```

**Key Methods:**
| Method | Description |
|--------|-------------|
| `RWMol.RemoveBond(beginIdx, endIdx)` | Remove bond between two atoms |
| `RWMol.AddBond(beginIdx, endIdx, BondType)` | Add new bond |
| `RWMol.RemoveAtom(idx)` | Remove atom and all its bonds |
| `RWMol.ReplaceAtom(idx, newAtom)` | Replace atom |
| `RWMol.ReplaceBond(idx, newBond)` | Replace bond |

---

### 1.2 EditableMol - Alternative Approach

```python
from rdkit import Chem

mol = Chem.MolFromSmiles('CCOCC')
edmol = Chem.EditableMol(mol)

# Remove bond
edmol.RemoveBond(1, 2)

# Add bond
edmol.AddBond(2, 4, order=Chem.rdchem.BondType.SINGLE)

# Convert back to Mol
new_mol = edmol.GetMol()
```

---

### 1.3 SMARTS Reactions for Bond Cleavage

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Define bond cleavage reaction using SMARTS
# Break C-C single bond
cc_bond_break = AllChem.ReactionFromSmarts('[C:1][C:2]>>[C:1].[C:2]')

# Run reaction
products = cc_bond_break.RunReactants((mol,))
for product in products[0]:
    print(Chem.MolToSmiles(product))

# Track broken bonds using properties:
# - `_ReactionDegreeChanged`: atoms where degree changed
# - `react_atom_idx`: index of atom in reacting molecule
```

---

### 1.4 FragmentOnBonds - Multiple Bond Breaking

```python
from rdkit import Chem

mol = Chem.MolFromSmiles("O-C-C-C-C-N")

# Fragment at bond indices, adds dummy atoms (*) at cleavage sites
mol_f = Chem.FragmentOnBonds(mol, (0, 2, 4), addDummies=True)

# Split into individual fragments
fragments = Chem.GetMolFrags(mol_f, asMols=True)
# Output: fragments as separate Mol objects

# Get atom indices for each fragment
atom_indices = Chem.GetMolFrags(mol_f, asMols=False)
```

---

### 1.5 FragmentOnSomeBonds - Permutations

```python
# Break exactly N bonds at a time - returns all combinations
result = Chem.FragmentOnSomeBonds(mol, bondIndices, numToBreak=2)

for fragments in result:
    print(Chem.MolToSmiles(fragments))
```

---

### 1.6 BRICS (Breaking of Retrosynthetically Interesting Chemical Substructures)

```python
from rdkit import Chem
from rdkit.Chem import BRICS, Recap

mol = Chem.MolFromSmiles("CCOCCC")

# Find chemically cleavable bonds
brics_bonds = list(BRICS.FindBRICSBonds(mol))
# Returns list of (atom_idx1, atom_idx2) tuples

# Fragment using BRICS
fragments = BRICS.BRICSBuild(mol)

# RECAP decomposition
recap = Recap.RecapDecompose(mol)
# Returns hierarchical fragment tree
```

---

### 1.7 Delete Substructures

```python
from rdkit import Chem

# Delete specific substructure
mol = Chem.MolFromSmiles('ClCCCl')
result = Chem.DeleteSubstructs(mol, Chem.MolFromSmiles('Cl'))

# Delete substructure and get fragments
mol_frags = Chem.GetMolFrags(result, asMols=True)
```

---

## Part 2: Atom Map Methods (Validation & Tracking)

### 2.1 Set/Get Atom Map Numbers

```python
from rdkit import Chem

mol = Chem.MolFromSmiles('CCO')

# Set atom map number (used in reaction mapping)
for atom in mol.GetAtoms():
    atom.SetAtomMapNum(atom.GetIdx() + 1)  # Start from 1

# Get atom map number
for atom in mol.GetAtoms():
    map_num = atom.GetAtomMapNum()
    print(f"Atom {atom.GetIdx()}: map {map_num}")

# Check if atom has map number
atom.HasProp('molAtomMapNumber')  # or: GetAtomMapNum() != 0
```

### 2.2 Clear Atom Map Numbers

```python
# Clear all atom map numbers
for atom in mol.GetAtoms():
    atom.SetAtomMapNum(0)

# Or use molecule-level clear
mol.ClearProp('molAtomMapNumber')
```

### 2.3 Substructure Matching for Mapping

```python
# Get single substructure match (returns atom indices in query order)
mol = Chem.MolFromSmiles('c1ccccc1O')
patt = Chem.MolFromSmarts('ccO')
match = mol.GetSubstructMatch(patt)
# Returns: (0, 5, 6)

# Get all matches
all_matches = mol.GetSubstructMatches(patt)
# Returns: ((0, 5, 6), (4, 5, 6))

# Quick validation
has_match = mol.HasSubstructMatch(patt)
```

### 2.4 Maximum Common Substructure (MCS) Mapping

```python
from rdkit import Chem
from rdkit.Chem import rdFMCS

mol1 = Chem.MolFromSmiles('CCO')
mol2 = Chem.MolFromSmiles('CCCO')

# Find MCS
params = rdFMCS.MCSParameters()
params.AtomTyper = rdFMCS.AtomCompare.CompareElements
params.BondTyper = rdFMCS.BondCompare.CompareOrder

res = rdFMCS.FindMCS([mol1, mol2], params)
# res.queryMol contains the MCS as a molecule

# Get atom mappings
map1 = mol1.GetSubstructMatch(res.queryMol)
map2 = mol2.GetSubstructMatch(res.queryMol)

# Create mapping as list of tuples
atom_mapping = list(zip(map1, map2))
# Output: [(0, 0), (1, 1), (2, 2), ...]
```

### 2.5 Utility Functions (MayaChemTools)

```python
# Check if atom map numbers present
has_maps = any(atom.GetAtomMapNum() != 0 for atom in mol.GetAtoms())

# Get atom indices sorted by atom map numbers
def get_atom_map_indices(mol):
    indices = []
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() > 0:
            indices.append((atom.GetAtomMapNum(), atom.GetIdx()))
    return [idx for _, idx in sorted(indices)]
```

---

## Part 3: Complete LLM/Agent Workflow

### 3.1 MoleculeEditor Class

```python
from rdkit import Chem
from rdkit.Chem import AllChem, BRICS, rdFMCS

class MoleculeEditor:
    def __init__(self, smiles):
        """Initialize with SMILES string"""
        self.original = Chem.MolFromSmiles(smiles)
        self.current = Chem.Mol(self.original)
        self.atom_mapping_history = []
    
    # ============= BOND OPERATIONS =============
    
    def disconnect_bond(self, atom1_idx, atom2_idx):
        """Disconnect bond between two atoms"""
        rwmol = Chem.RWMol(self.current)
        rwmol.RemoveBond(atom1_idx, atom2_idx)
        self.current = rwmol.GetMol()
        Chem.SanitizeMol(self.current)
        self.atom_mapping_history.append({
            'action': 'disconnect',
            'atoms': (atom1_idx, atom2_idx)
        })
        return self
    
    def disconnect_bonds(self, bond_list):
        """Disconnect multiple bonds"""
        for a1, a2 in bond_list:
            self.disconnect_bond(a1, a2)
        return self
    
    def reconnect_bond(self, atom1_idx, atom2_idx, bond_type='SINGLE'):
        """Reconnect two atoms with a bond"""
        bt = getattr(Chem.rdchem.BondType, bond_type)
        rwmol = Chem.RWMol(self.current)
        rwmol.AddBond(atom1_idx, atom2_idx, bt)
        self.current = rwmol.GetMol()
        Chem.SanitizeMol(self.current)
        return self
    
    # ============= FRAGMENT OPERATIONS =============
    
    def get_fragments(self):
        """Get all fragments after bond breaking"""
        return Chem.GetMolFrags(self.current, asMols=True)
    
    def get_fragment_smiles(self):
        """Get SMILES for all fragments"""
        frags = self.get_fragments()
        return [Chem.MolToSmiles(f) for f in frags]
    
    def find_brics_bonds(self):
        """Find BRICS cleavable bonds"""
        return list(BRICS.FindBRICSBonds(self.current))
    
    # ============= VALIDATION =============
    
    def validate_structure(self):
        """Validate molecule is chemically correct"""
        try:
            Chem.SanitizeMol(self.current)
            return True, "Valid molecule"
        except Exception as e:
            return False, str(e)
    
    def is_valid(self):
        """Quick validity check"""
        return self.validate_structure()[0]
    
    # ============= MAPPING =============
    
    def set_atom_maps(self):
        """Set atom map numbers (1-based)"""
        for atom in self.current.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        return self
    
    def clear_atom_maps(self):
        """Clear all atom map numbers"""
        for atom in self.current.GetAtoms():
            atom.SetAtomMapNum(0)
        return self
    
    def get_atom_mapping(self):
        """Get current atom map as dict"""
        return {atom.GetIdx(): atom.GetAtomMapNum() 
                for atom in self.current.GetAtoms()}
    
    def validate_against(self, target_mol):
        """Validate mapping against target molecule using MCS"""
        params = rdFMCS.MCSParameters()
        params.AtomTyper = rdFMCS.AtomCompare.CompareElements
        
        res = rdFMCS.FindMCS([self.current, target_mol], params)
        
        if res.queryMol is None:
            return False, "No MCS found", None
        
        map1 = self.current.GetSubstructMatch(res map2 = target_mol.GetSubstructMatch(res.queryMol)
        
       .queryMol)
        mapping = list(zip(map1, map2))
        return True, f"MCS atoms: {len(mapping)}", mapping
    
    # ============= OUTPUT =============
    
    def get_smiles(self):
        """Get canonical SMILES"""
        return Chem.MolToSmiles(self.current)
    
    def get_smarts(self):
        """Get SMARTS representation"""
        return Chem.MolToSmarts(self.current)
    
    def get_2d_coords(self):
        """Generate 2D coordinates for visualization"""
        AllChem.Compute2DCoords(self.current)
        return self
    
    def copy(self):
        """Create a copy of current editor"""
        new_editor = MoleculeEditor(self.get_smiles())
        new_editor.original = Chem.Mol(self.original)
        return new_editor
```

### 3.2 Usage Example

```python
# Initialize
editor = MoleculeEditor('CCOCC')

# Set atom maps for tracking
editor.set_atom_maps()

# Disconnect a bond
editor.disconnect_bond(1, 2)

# Validate
print(f"Valid: {editor.is_valid()}")
print(f"SMILES: {editor.get_smiles()}")

# Get fragments
frags = editor.get_fragment_smiles()
print(f"Fragments: {frags}")

# Validate against target
target = Chem.MolFromSmiles('CCO')
valid, msg, mapping = editor.validate_against(target)
print(f"Validation: {msg}")
```

---

## Part 4: Key Functions Summary Table

| Purpose | Function/Method | Module |
|---------|----------------|--------|
| Remove bond | `RWMol.RemoveBond(idx1, idx2)` | rdkit.Chem.rdchem |
| Add bond | `RWMol.AddBond(idx1, idx2, BondType)` | rdkit.Chem.rdchem |
| Remove atom | `RWMol.RemoveAtom(idx)` | rdkit.Chem.rdchem |
| Delete substructure | `Chem.DeleteSubstructs(mol, pattern)` | rdkit.Chem |
| Fragment on bonds | `Chem.FragmentOnBonds(mol, bondIndices)` | rdkit.Chem |
| Fragment permutations | `Chem.FragmentOnSomeBonds(mol, indices, n)` | rdkit.Chem |
| Get fragments | `Chem.GetMolFrags(mol, asMols=True)` | rdkit.Chem |
| BRICS bonds | `BRICS.FindBRICSBonds(mol)` | rdkit.Chem.BRICS |
| RECAP decompose | `Recap.RecapDecompose(mol)` | rdkit.Chem.Recap |
| SMARTS reaction | `AllChem.ReactionFromSmarts(smarts)` | rdkit.Chem.AllChem |
| Set atom map | `atom.SetAtomMapNum(num)` | rdkit.Chem.rdchem |
| Get atom map | `atom.GetAtomMapNum()` | rdkit.Chem.rdchem |
| Substruct match | `mol.GetSubstructMatch(query)` | rdkit.Chem |
| MCS find | `rdFMCS.FindMCS(mols, params)` | rdkit.Chem.rdFMCS |
| Validate/sanitize | `Chem.SanitizeMol(mol)` | rdkit.Chem |
| To SMILES | `Chem.MolToSmiles(mol)` | rdkit.Chem.rdmolfiles |
| To SMARTS | `Chem.MolToSmarts(mol)` | rdkit.Chem.rdmolfiles |

---

## Part 5: Related LLM/Agent Frameworks

For integrating with LLM/agent workflows:

| Project | Description | Use Case |
|---------|-------------|----------|
| **ChemGraph** | LangGraph + RDKit + ASE | Automated simulation |
| **DrugAgent** | Multi-agent collaboration | Drug discovery |
| **Molpert** | Graph-based perturbation | Molecule editing |
| **LocalRetro** | Retrosynthesis with atom mapping | Reaction planning |

---

## References

- RDKit Documentation: https://www.rdkit.org/docs/
- RDKit Cookbook: https://www.rdkit.org/docs/Cookbook.html
- RDKit GitHub: https://github.com/rdkit/rdkit
- MayaChemTools: https://www.mayachemtools.org/

---

*Last Updated: 2026-02-25*
