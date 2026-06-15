import logging
from typing import Dict, List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

logger = logging.getLogger("ChemicalReactor.Utils")

try:
    from rdkit.Chem.rdMolStandardize import SaltRemover
except ImportError:
    SaltRemover = None

# =====================================================================
# THE DICTIONARY DEFINITIONS (58 Library Entities)
# =====================================================================
DEFAULT_REACTIONS = {
    "Schotten-Baumann_amide": "[C:1](=[O:2])[O;H1].[N;!H0;X3;!$(N-[!#6;!#1]):3]>>[C:1](=[O:2])[N:3]",
    "1,2,4-triazole_acetohydrazide": "[C&H0&$(C-[#6]):1]#[N&H0:2].[N&H2:3]-[N&H1:4]-[C&H0&$(C-[#6])&R0:5]=[O&D1]>>[N:2]1-[C:1]=[N:3]-[N:4]-[C:5]=1",
    "1,2,4-triazole_carboxylic-acid/ester": "[C&H0&$(C-[#6]):1]#[N&H0:2].[C&H0&$(C-[#6])&R0:5](=[O&D1])-[#8;H1,$(O-[C&H3]),$(O-[C&H2]-[C&H3])]>>[N:2]1-[C:1]=N-N-[C:5]=1",
    "3-nitrile-pyridine": "[#6&!$([#6](-C=O)-C=O):4]-[C&H0:1](=[O&D1])-[C;H1&!$(C-[!#6])&!$(C-C(=O)O),H2:2]-[C&H0&R0:3](=[O&D1])-[#6&!$([#6](-C=O)-C=O):5]>>[c:1]1(-[#6:4]):[c:2]:[c:3](-[#6:5]):n:c(-O):c:1-C#N",
    "Buchwald-Hartwig": "[Cl,Br,I][c&$(c1:[c,n]:[c,n]:[c,n]:[c,n]:[c,n]:1):1].[N;$(NC)&!$(N=*)&!$([N&-])&!$(N#*)&!$([N&D3])&!$([N&D4])&!$(N[c,O])&!$(N[C,S]=[S,O,N]),H2&$(Nc1:[c,n]:[c,n]:[c,n]:[c,n]:[c,n]:1):2]>>[c:1][N:2]",
    "Fischer indole": "[N&H1&$(N-c1ccccc1):1](-[N&H2])-[c:5]:[c&H1:4].[C&$(C(-,:[#6])[#6]):2](=[O&D1])-[C&H2&$(C(-,:[#6])[#6])&!$(C(-,:C=O)C=O):3]>>[C:5]1-[N:1]-[C:2]=[C:3]-[C:4]:1",
    "Friedlaender chinoline": "[N&H2&$(N-c1ccccc1):1]-[c:2]:[c:3]-[C&H1:4]=[O&D1].[C&$(C(-,:[#6])[#6]):6](=[O&D1])-[C&H2&$(C(-,:[#6])[#6])&!$(C(-,:C=O)C=O):5]>>[N:1]1-[c:2]:[c:3]-[C:4]=[C:5]-[C:6]:1",
    "Grignard_alcohol": "[#6:1][C;H1,$(C(-,:[#6])[#6]):2]=[O&D1:3].[Cl,Br,I][#6&$([#6]~[#6])&!$([#6](-,:[Cl,Br,I])[Cl,Br,I])&!$([#6]=O):4]>>[C:1][#6:2](-,:[O&H1:3])[#6:4]",
    "Grignard_carbonyl": "[#6:1][C:2]#[#7&D1].[Cl,Br,I][#6&$([#6]~[#6])&!$([#6](-,:[Cl,Br,I])[Cl,Br,I])&!$([#6]=O):3]>>[#6:1][C:2](=O)[#6:3]",
    "Heck_non-terminal_vinyl": "[#6;c,$(C(=O)O),$(C#N):3][#6:2](-,:[#6:5])=[#6&H1&$([#6][#6]):1].[#6;$([#6]=[#6]),$(c:c):4][Cl,Br,I]>>[#6:4][#6&H0:1]=[#6:2](-,:[#6:5])[#6:3]",
    "Heck_terminal_vinyl": "[#6;c,$(C(=O)O),$(C#N):3][#6&H1:2]=[#6&H2:1].[#6;$([#6]=[#6]),$(c:c):4][Cl,Br,I]>>[#6:4]/[#6:1]=[#6:2]/[#6:3]",
    "Huisgen_Cu-catalyzed_1,4-subst": "[C&H0&$(C-[#6]):1]#[C&H1:2].[C;H1,H2;A;!$(C=O):3]-[#17,#35,#53,O&H1]>>[C:1]1=[C:2]-N(-[C:3])-N=N-1",
    "Huisgen_Ru-catalyzed_1,5_subst": "[C&H0&$(C-[#6]):1]#[C&H1:2].[C;H1,H2;A;!$(C=O):3]-[#17,#35,#53,O&H1]>>[C:1]1=[C:2]-N=NN-1-[C:3]",
    "Huisgen_disubst-alkyne": "[C&H0&$(C-[#6]):1]#[C&H0&$(C-[#6]):2].[C;H1,H2;A;!$(C=O):3]-[#17,#35,#53,O&H1]>>[C:1]1=[C:2]-N=NN-1-[C:3]",
    "Mitsunobu_imide": "[C;H1&$(C(-,:[#6])[#6]),H2&$(C[#6]):1][O&H1].[N&H1&$(N(-,:C=O)C=O):2]>>[C:1][N:2]",
    "Mitsunobu_phenole": "[C;H1&$(C(-,:[#6])[#6]),H2&$(C[#6]):1][O&H1].[O&H1&$(Oc1ccccc1):2]>>[C:1][O:2]",
    "Mitsunobu_sulfonamide": "[C;H1&$(C(-,:[#6])[#6]),H2&$(C[#6]):1][O&H1].[N&H1&$(N(-,:[#6])S(=O)=O):2]>>[C:1][N:2]",
    "Mitsunobu_tetrazole_1": "[C;H1&$(C(-,:[#6])[#6]),H2&$(C[#6]):1][O&H1].[#7&H1:2]1~[#7:3]~[#7:4]~[#7:5]~[#6:6]~1>>[C:1][#7:2]1:[#7:3]:[#7:4]:[#7:5]:[#6:6]:1",
    "Mitsunobu_tetrazole_2": "[C;H1&$(C(-,:[#6])[#6]),H2&$(C[#6]):1][O&H1].[#7&H1:2]1~[#7:3]~[#7:4]~[#7:5]~[#6:6]~1>>[#7&H0:2]1:[#7:3]:[#7&H0:4](-,:[C:1]):[#7:5]:[#6:6]:1",
    "Mitsunobu_tetrazole_3": "[C;H1&$(C(-,:[#6])[#6]),H2&$(C[#6]):1][O&H1].[#7:2]1~[#7:3]~[#7&H1:4]~[#7:5]~[#6:6]~1>>[C:1][#7&H0:2]1:[#7:3]:[#7&H0:4]:[#7:5]:[#6:6]:1",
    "Mitsunobu_tetrazole_4": "[C;H1&$(C(-,:[#6])[#6]),H2&$(C[#6]):1][O&H1].[#7:2]1~[#7:3]~[#7&H1:4]~[#7:5]~[#6:6]~1>>[#7:2]1:[#7:3]:[#7:4](-,:[C:1]):[#7:5]:[#6:6]:1",
    "N-alkylation_heterocycles": "[C:1][Cl,Br,I].[n&H1&+0&r5&!$(n[#6]=[O,S,N])&!$(n~n~n)&!$(n~n~c~n)&!$(n~c~n~n):2]>>[C:1][n:2]",
    "N-arylation_heterocycles": "[c:1]B(-,:O)O.[n&H1&+0&r5&!$(n[#6]=[O,S,N])&!$(n~n~n)&!$(n~n~c~n)&!$(n~c~n~n):2]>>[c:1][n:2]",
    "Negishi": "[#6&$([#6]~[#6])&!$([#6]~[S,N,O,P]):1][Cl,Br,I].[Cl,Br,I][#6&$([#6]~[#6])&!$([#6]~[S,N,O,P]):2]>>[#6:2][#6:1]",
    "Niementowski_quinazoline": "[c:1](-[C&$(C-c1ccccc1):2](=[O&D1:3])-[O&H1]):[c:4]-[N&H2:5].[N&!H0&!$(N-N)&!$(N-C=N)&!$(N(-C=O)-C=O):6]-[C;H1,$(C-[#6]):7]=[O&D1]>>[c:4]1:[c:1]-[C:2](=[O:3])-[N:6]-[C:7]=[N:5]-1",
    "Paal-Knorr pyrrole": "[#6:5]-[C&R0:1](=[O&D1])-[C;H1,H2:2]-[C;H1,H2:3]-[C:4](=[O&D1])-[#6:6].[N&H2&$(N-[C,N])&!$(NC=[O,S,N])&!$(N(-,:[#6])[#6])&!$(N~N~N):7]>>[C:1]1(-[#6:5])=[C:2]-[C:3]=[C:4](-[#6:6])-[N:7]-1",
    "Pictet-Spengler": "[c&H1:1]1:[c:2](-[C&H2:7]-[C&H2:8]-[N&H2:9]):[c:3]:[c:4]:[c:5]:[c:6]:1.[#6:11]-[C&H1&R0:10]=[O&D1]>>[c:1]12:[c:2](-[C&H2:7]-[C&H2:8]-[N&H1:9]-[C:10]-1-[#6:11]):[c:3]:[c:4]:[c:5]:[c:6]:2",
    "SNAr": "[N&X3&!H0&!$(N-[C,S]=O):1].[F,Cl]-[c&R1:2]>>[N&X3:1]-[c&R1:2]",
    "Sonogashira": "[c&R1:3][Br,Cl,I].[C&H1&X2:2]#[C&X2&$(C(#C)a):1]>>[c&R1:3][C&H0&X2:2]#[C&X2:1]",
    "Stille": "[#6;$(C=C-[#6]),$(c:c):1][Br,I].[Cl,Br,I][c:2]>>[c:2][#6:1]",
    "Suzuki": "[#6&H0&D3&$([#6](~[#6])~[#6]):1]B(-,:O)O.[#6&H0&D3&$([#6](~[#6])~[#6]):2][Cl,Br,I]>>[#6:2][#6:1]",
    "Williamson ether": "[#6&$([#6]~[#6])&!$([#6]=O):2][#8&H1:3].[Cl,Br,I][#6&H2&$([#6]~[#6]):4]>>[C&H2:4][O:3][#6:2]",
    "Wittig": "[#6:3]-[C;H1,$([C&H0](-[#6])[#6]);!$(CC=O):1]=[O&D1].[Cl,Br,I][C&H2&$(C-[#6])&!$(CC[I,Br])&!$(CCO[C&H3]):2]>>[C:3][C:1]=[C:2]",
    "alkyl_halide": "[C&X4&H2:1][Br,Cl,I].[O&H1&$(Oa):2]>>[C&X4&H2:1]-[O&H0:2]",
    "benzimidazole_derivatives_aldehyde": "[c&r6:1](-[N&H1&$(N-[#6]):2]):[c&r6:3]-[N&H2:4].[#6:6]-[C&H1&R0:5]=[O&D1]>>[c:3]1:[c:1]:[n:2]:[c:5](-[#6:6]):[n:4]@1",
    "benzimidazole_derivatives_carboxylic-acid/ester": "[c&r6:1](-[N&H1&$(N-[#6]):2]):[c&r6:3]-[N&H2:4].[#6:6]-[C&R0:5](=[O&D1])-[#8;H1,$(O-[C&H3])]>>[c:3]1:[c:1]:[n:2]:[c:5](-[#6:6]):[n:4]@1",
    "benzofuran": "[Br,I;$(*c1ccccc1)]-[c:1]:[c:2]-[O&H1:3].[C&H1:5]#[C&$(C-[#6]):4]>>[c:1]1:[c:2]-[O:3]-[C:4]=[C:5]-1",
    "benzothiazole": "[c&r6:1](-[S&H1:2]):[c&r6:3]-[N&H2:4].[#6:6]-[C&H1&R0:5]=[O&D1]>>[c:3]1:[c:1]:[s:2]:[c:5](-[#6:6]):[n:4]@1",
    "benzothiophene": "[Br,I;$(*c1ccccc1)]-[c:1]:[c:2]-[S&D2:3]-[C&H3].[C&H1:5]#[C&$(C-[#6]):4]>>[c:1]1:[c:2]-[S:3]-[C:4]=[C:5]-1",
    "benzoxazole_arom-aldehyde": "[c:1](-[O&H1&$(Oc1ccccc1):2]):[c&r6:3]-[N&H2:4].[c:6]-[C&H1&R0:5]=[O&D1]>>[c:3]1:[c:1]:[o:2]:[c:5](-[c:6]):[n:4]@1",
    "benzoxazole_carboxylic-acid": "[c&r6:1](-[O&H1:2]):[c&r6:3]-[N&H2:4].[#6:6]-[C&R0:5](=[O&D1])-[O&H1]>>[c:3]1:[c:1]:[o:2]:[c:5](-[#6:6]):[n:4]@1",
    "decarboxylative_coupling": "[c&$(c1[c&$(c[C,S,N](=[O&D1])[R0&!O&H1])]cccc1):1][C&$(C(=O)[O&H1])].[c&$(c1aaccc1):2][Cl,Br,I]>>[c:1][c:2]",
    "heteroaromatic_nuc_sub": "[c&!$(c1ccccc1)&$(c1[n,c]c[n,c]c[n,c]1):1][Cl,F].[N&$(NC)&!$(N=*)&!$([N&-])&!$(N#*)&!$([N&D3])&!$([N&D4])&!$(N[c,O])&!$(N[C,S]=[S,O,N]):2]>>[c:1][N:2]",
    "imidazole": "[C&$(C(-,:[#6])[#6&!$([#6]Br)]):4](=[O&D1])[C&H1&$(C(-,:[#6])[#6]):5]Br.[#7&H2:3][C&$(C(=N)(-,:N)[c,#7]):2]=[#7&H1&D1:1]>>[C:4]1=[C&H0:5][N&H1:3][C:2]=[N:1]1",
    "indole": "[Br,I;$(*c1ccccc1)]-[c:1]:[c:2]-[N&H2:3].[C&H1:5]#[C&$(C-[#6]):4]>>[c:1]1:[c:2]-[N:3]-[C:4]=[C:5]-1",
    "nucl_sub_aromatic_ortho_nitro": "[c&$(c1c(-,:N(~O)~O)cccc1):1][Cl,F].[N&$(NC)&!$(N=*)&!$([N&-])&!$(N#*)&!$([N&D3])&!$([N&D4])&!$(N[c,O])&!$(N[C,S]=[S,O,N]):2]>>[c:1][N:2]",
    "nucl_sub_aromatic_para_nitro": "[c&$(c1ccc(-,:N(~O)~O)cc1):1][Cl,F].[N&$(NC)&!$(N=*)&!$([N&-])&!$(N#*)&!$([N&D3])&!$([N&D4])&!$(N[c,O])&!$(N[C,S]=[S,O,N]):2]>>[c:1][N:2]",
    "oxadiazole": "[#6:6][C:5]#[#7&D1:4].[#6:1][C:2](=[O&D1:3])[O&H1]>>[#6:6][c:5]1[n:4][o:3][c:2](-,:[#6:1])n1",
    "phthalazinone": "[c&r6:1](-[C&$(C=O):6]-[O&H1]):[c&r6:2]-[C;H1,$(C-C):3]=[O&D1].[N&H2:4]-[N&H1&$(N-[#6])&!$(NC=[O,S,N]):5]>>[c:1]1:[c:2]-[C:3]=[N:4]-[N:5]-[C:6]-1",
    "piperidine_indole": "[c&H1:3]1:[c:4]:[c:5]:[c&H1:6]:[c:7]2:[n&H1:8]:[c:9]:[c&H1:1]:[c:2]:1:2.O=[C:10]1[#6&H2:11][#6&H2:12][N:13][#6&H2:14][#6&H2:15]1>>[#6&H2:12]1[#6&H1:11]=[C:10](-,:[c:1]2:[c:9]:[n:8]:[c:7]3:[c:6]:[c:5]:[c:4]:[c:3]:[c:2]:2:3)[#6&H2:15][#6&H2:14][N:13]1",
    "pyrazole": "[#6&!$([#6](-C=O)-C=O):4]-[C&H0:1](=[O&D1])-[C;H1&!$(C-[!#6])&!$(C-C(=O)O),H2:2]-[C&H0&R0:3](=[O&D1])-[#6&!$([#6](-C=O)-C=O):5].[N&H2:6]-[N&!H0;$(N-[#6]),H2:7]>>[C:1]1(-[#6:4])-[C:2]=[C:3](-[#6:5])-[N:7]-[N:6]=1",
    "reductive amination": "[#6:4]-[C;H1,$([C&H0](-[#6])[#6]):1]=[O&D1].[N;H2,$([N&H1&D2](-,:C)C);!$(N-[#6]=*):3]-[C:5]>>[#6:4][C:1]-[N:3]-[C:5]",
    "spiro-chromanone": "[c:1](-[C&$(C-c1ccccc1):2](=[O&D1:3])-[C&H3:4]):[c:5]-[O&H1:6].[C&$(C1-[C&H2]-[C&H2]-[N,C]-[C&H2]-[C&H2]-1):7]=[O&D1]>>[O:6]1-[c:5]:[c:1]-[C:2](=[O&D1:3])-[C:4]-[C:7]-1",
    "sulfon_amide": "[S&$(S(=O)(=O)[C,N]):1]Cl.[N&$(NC)&!$(N=*)&!$([N&-])&!$(N#*)&!$([N&D3])&!$([N&D4])&!$(N[c,O])&!$(N[C,S]=[S,O,N]):2]>>[S:1][N&+0:2]",
    "tetrazole_connect_regioisomere_1": "[C&H0&$(C-[#6]):1]#[N&H0:2].[C&A&!$(C=O):3]-[#17,#35,#53]>>[C:1]1=[N:2]-N(-[C:3])-N=N-1",
    "tetrazole_connect_regioisomere_2": "[C&H0&$(C-[#6]):1]#[N&H0:2].[C&A&!$(C=O):3]-[#17,#35,#53]>>[C:1]1=[N:2]-N=N-N-1-[C:3]",
    "tetrazole_terminal": "[C&H0&$(C-[#6]):1]#[N&H0:2]>>[C:1]1=[N:2]-N-N=N-1",
    "thiazole": "[#6:6]-[C&R0:1](=[O&D1])-[C&H1&R0:5](-[#6:7])-[#17,#35,#53].[N&H2:2]-[C:3]=[S&D1:4]>>[c:1]1(-[#6:6]):[n:2]:[c:3]:[s:4][c:5]:1-,:[#6:7]",
    "thiourea": "[N&$(N-[#6]):3]=[C&$(C=S):1].[N&$(N[#6])&!$(N=*)&!$([N&-])&!$(N#*)&!$([N&D3])&!$([N&D4])&!$(N[O,N])&!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N&+0:2]",
    "triaryl-imidazole": "[C&$(C-c1ccccc1):1](=[O&D1])-[C&D3&$(C-c1ccccc1):2]~[O;D1,H1].[C&H1&$(C-c):3]=[O&D1]>>[C:1]1-N=[C:3]-[N&H1]-[C:2]=1",
    "urea": "[N&$(N-[#6]):3]=[C&$(C=O):1].[N&$(N[#6])&!$(N=*)&!$([N&-])&!$(N#*)&!$([N&D3])&!$([N&D4])&!$(N[O,N])&!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N&+0:2]"
}

class StructureEngine:
    def __init__(self):
        if SaltRemover is not None:
            try:
                self.remover = SaltRemover()
            except Exception:
                self.remover = None
        else:
            self.remover = None

    def clean_and_sanitize(self, smiles: str) -> Optional[Chem.Mol]:
        try:
            smiles = str(smiles).strip()
            if not smiles or smiles == "nan":
                return None
            if "." in smiles and self.remover is None:
                smiles = max(smiles.split("."), key=len)
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            if self.remover is not None:
                mol = self.remover.StripMol(mol, dontRemoveEverything=True)
            mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(mol)
            return mol
        except Exception:
            return None

    @staticmethod
    def calculate_properties(mol: Chem.Mol) -> Dict[str, float]:
        return {
            "MW": round(float(Descriptors.MolWt(mol)),1),
            "LogP": round(float(Descriptors.MolLogP(mol)),1),
            "HBD": int(Descriptors.NumHDonors(mol)),
            "HBA": int(Descriptors.NumHAcceptors(mol)),
            "RotBonds": int(Descriptors.NumRotatableBonds(mol))
        }

class Reaction:
    def __init__(self, name: str, smarts: str):
        self.name = name
        self.smarts = smarts
        self._rxn = AllChem.ReactionFromSmarts(smarts)

    def run_bimolecular(self, mol_a: Chem.Mol, mol_b: Chem.Mol) -> Optional[Chem.Mol]:
        """Runs transformation cross-checking BOTH reactant alignment configurations dynamically."""
        if mol_a is None or mol_b is None:
            return None

        # Try Alignment 1: (mol_a, mol_b)
        try:
            ps = self._rxn.RunReactants((mol_a, mol_b))
            if ps and len(ps) > 0 and len(ps[0]) > 0:
                return ps[0][0]
        except Exception:
            pass
            
        # Try Alignment 2: (mol_b, mol_a) - handles flipped positional layouts
        try:
            ps = self._rxn.RunReactants((mol_b, mol_a))
            if ps and len(ps) > 0 and len(ps[0]) > 0:
                return ps[0][0]
        except Exception:
            pass

        return None
