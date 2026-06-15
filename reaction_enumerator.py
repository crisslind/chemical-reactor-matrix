#!/usr/bin/env python3
import os
import json
import logging
import argparse
from typing import Dict, List, Optional, Union, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from utils import DEFAULT_REACTIONS, StructureEngine, Reaction

logging.basicConfig(
    format="%(asctime)s - [%(levelname)s] - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("ChemicalReactor")
logger.setLevel(logging.INFO)

MAX_WORKERS = os.cpu_count() or 1

def _parallel_worker(
    scaffold_id: str, 
    scaffold_smiles: str, 
    bb_catalog: Dict[str, str], 
    reactions_subset: Dict[str, str],
    max_mw: float,
    max_logp: float
) -> List[Dict]:
    """Processes a single scaffold and collapses internal product duplicates."""
    engine = StructureEngine()
    # Key: Product_SMILES, Value: row data dictionary
    local_duplicates: Dict[str, Dict] = {}
    
    scaffold_mol = engine.clean_and_sanitize(scaffold_smiles)
    if not scaffold_mol:
        return []

    rxn_objs = [Reaction(name, smarts) for name, smarts in reactions_subset.items()]
    
    parsed_bbs = []
    for bb_id, bb_smi in bb_catalog.items():
        bb_mol = engine.clean_and_sanitize(bb_smi)
        if bb_mol:
            canonical_bb_smi = Chem.MolToSmiles(bb_mol, canonical=True)
            parsed_bbs.append((str(bb_id), canonical_bb_smi, bb_mol))

    scaffold_canonical = Chem.MolToSmiles(scaffold_mol, canonical=True)

    for rxn in rxn_objs:
        for bb_id, bb_smi_str, bb_mol in parsed_bbs:
            prod_mol = rxn.run_bimolecular(scaffold_mol, bb_mol)
            if not prod_mol:
                continue
                
            try:
                prod_mol.UpdatePropertyCache(strict=False)
                Chem.SanitizeMol(prod_mol)
                prod_smiles = Chem.MolToSmiles(prod_mol, canonical=True)
                
                props = engine.calculate_properties(prod_mol)
                if props["MW"] > max_mw or props["LogP"] > max_logp:
                    continue
                
                # Deduplication logic inside the worker
                if prod_smiles in local_duplicates:
                    existing = local_duplicates[prod_smiles]
                    # Append unique BB IDs and SMILES using a semicolon delimiter
                    if bb_id not in existing["Building_Block_ID"].split("; "):
                        existing["Building_Block_ID"] += f"; {bb_id}"
                    if bb_smi_str not in existing["Building_Block_SMILES"].split("; "):
                        existing["Building_Block_SMILES"] += f"; {bb_smi_str}"
                else:
                    row = {
                        "Scaffold_ID": scaffold_id,
                        "Input_Scaffold_SMILES": scaffold_canonical,
                        "Reaction_Executed": rxn.name,
                        "Building_Block_ID": bb_id,
                        "Building_Block_SMILES": bb_smi_str,
                        "Product_SMILES": prod_smiles
                    }
                    row.update(props)
                    local_duplicates[prod_smiles] = row
            except Exception:
                continue

    return list(local_duplicates.values())


class HighPerformanceReactor:
    def __init__(self, selected_reactions: Optional[List[str]] = None):
        if selected_reactions:
            self.rxn_subset = {k: v for k, v in DEFAULT_REACTIONS.items() if k in selected_reactions}
            if not self.rxn_subset:
                self.rxn_subset = DEFAULT_REACTIONS
        else:
            self.rxn_subset = DEFAULT_REACTIONS

    def _parse_scaffolds(self, source: str) -> List[Tuple[str, str]]:
        engine = StructureEngine()
        extracted = []
        
        if not os.path.exists(source):
            mol = engine.clean_and_sanitize(source)
            if mol:
                return [("Lead_Scaffold_0", Chem.MolToSmiles(mol, canonical=True))]
            raise FileNotFoundError(f"Input path or SMILES invalid: '{source}'")

        if source.lower().endswith(('.sdf', '.sd')):
            suppl = Chem.SDMolSupplier(source)
            for idx, mol in enumerate(suppl):
                if mol:
                    nm = mol.GetProp("_Name") if mol.HasProp("_Name") else f"SDF_{idx}"
                    extracted.append((nm, Chem.MolToSmiles(mol, canonical=True)))
        else:
            sep = '\t' if source.lower().endswith(('.tsv', '.txt')) else ','
            df = pd.read_csv(source, sep=sep)
            df.columns = [str(c).strip() for c in df.columns]
            
            smi_col = next((c for c in df.columns if 'smiles' in c.lower() or 'structure' in c.lower()), df.columns[0])
            id_col = next((c for c in df.columns if 'id' in c.lower() or 'name' in c.lower()), None)
            
            for idx, row in df.iterrows():
                smi = str(row[smi_col]).strip()
                lbl = str(row[id_col]).strip() if (id_col and pd.notna(row[id_col])) else f"Row_{idx}"
                extracted.append((lbl, smi))
                
        return extracted

    def _parse_catalog(self, source: Union[str, Dict[str, str]]) -> Dict[str, str]:
        if isinstance(source, dict):
            return source
        if not os.path.exists(source):
            raise FileNotFoundError(f"Reagent library catalog path invalid: '{source}'")
        if source.lower().endswith('.json'):
            with open(source, 'r') as f:
                return json.load(f)
                
        sep = '\t' if source.lower().endswith(('.tsv', '.txt')) else ','
        df = pd.read_csv(source, sep=sep)
        df.columns = [str(c).strip() for c in df.columns]
        
        smi_col = next((c for c in df.columns if 'smiles' in c.lower() or 'structure' in c.lower()), df.columns[0])
        id_col = next((c for c in df.columns if 'id' in c.lower() or 'name' in c.lower()), None)
        
        if id_col:
            return df.dropna(subset=[smi_col]).set_index(id_col)[smi_col].astype(str).to_dict()
        return {f"BB_{idx}": str(smi).strip() for idx, smi in enumerate(df[smi_col].dropna())}

    def execute_library(
        self, 
        scaffold_source: str, 
        bb_source: Union[str, Dict[str, str]], 
        max_mw: float = 2000.0,
        max_logp: float = 20.0,
        workers: int = MAX_WORKERS
    ) -> pd.DataFrame:
        scaffolds = self._parse_scaffolds(scaffold_source)
        catalog = self._parse_catalog(bb_source)
        
        logger.info(f"Loaded {len(scaffolds)} target scaffolds & {len(catalog)} building block reagents.")
        logger.info(f"Running library generation on {workers} parallel processing pool workers...")
        
        global_duplicates: Dict[str, Dict] = {}
        
        # Multiprocessing execution accompanied by tracking loop visualizer
        with ProcessPoolExecutor(max_workers=workers) as executor:
            future_to_scaffold = {
                executor.submit(
                    _parallel_worker, 
                    scaf_id, scaf_smi, catalog, self.rxn_subset, max_mw, max_logp
                ): scaf_id for scaf_id, scaf_smi in scaffolds
            }
            
            # Progress bar tracks completed scaffolds
            with tqdm(total=len(future_to_scaffold), desc="Synthesizing Library Matrices", unit="scaffold") as pbar:
                for future in as_completed(future_to_scaffold):
                    scaf_id = future_to_scaffold[future]
                    pbar.update(1)
                    try:
                        worker_rows = future.result()
                        if not worker_rows:
                            continue
                            
                        # Final cross-scaffold duplicate consolidation pass
                        for row in worker_rows:
                            prod_smi = row["Product_SMILES"]
                            if prod_smi in global_duplicates:
                                existing = global_duplicates[prod_smi]
                                
                                # Consolidate IDs across different scaffold lineages if needed
                                for chunk in row["Building_Block_ID"].split("; "):
                                    if chunk not in existing["Building_Block_ID"].split("; "):
                                        existing["Building_Block_ID"] += f"; {chunk}"
                                        
                                # Consolidate SMILES structures
                                for chunk in row["Building_Block_SMILES"].split("; "):
                                    if chunk not in existing["Building_Block_SMILES"].split("; "):
                                        existing["Building_Block_SMILES"] += f"; {chunk}"
                                        
                                # Append different Scaffold IDs if the exact same product is shared
                                if row["Scaffold_ID"] not in existing["Scaffold_ID"].split("; "):
                                    existing["Scaffold_ID"] += f"; {row['Scaffold_ID']}"
                            else:
                                global_duplicates[prod_smi] = row
                                
                    except Exception as exc:
                        logger.error(f"Worker tracking line item '{scaf_id}' raised unexpected panic: {exc}")

        master_records = list(global_duplicates.values())
        logger.info(f"Execution complete. Total unique structures synthesized: {len(master_records)}")
        return pd.DataFrame(master_records)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parallel Chemical Reaction Enumeration Matrix.")
    parser.add_argument("-s", "--scaffolds", required=True, type=str, help="Single SMILES or file path (.csv/.tsv/.sdf)")
    parser.add_argument("-b", "--building_blocks", required=True, type=str, help="Catalog file path (.csv/.tsv/.json)")
    parser.add_argument("-o", "--output", default="library_output.csv", type=str, help="Output target path destination for results CSV")
    parser.add_argument("-r", "--reactions", nargs="*", default=None, help="Optional functional reaction subset filters")
    parser.add_argument("--max_mw", default=2000.0, type=float, help="Maximum allowed molecular weight filtering cutoff")
    parser.add_argument("--max_logp", default=20.0, type=float, help="Maximum allowed LogP filtering cutoff")
    parser.add_argument("-w", "--workers", default=MAX_WORKERS, type=int, help="Override default multiprocessing allocation values")
    
    args = parser.parse_args()
    
    reactor = HighPerformanceReactor(selected_reactions=args.reactions)
    results_df = reactor.execute_library(
        scaffold_source=args.scaffolds,
        bb_source=args.building_blocks,
        max_mw=args.max_mw,
        max_logp=args.max_logp,
        workers=args.workers
    )
    
    results_df.to_csv(args.output, index=False)
    logger.info(f"Successfully exported data to location matrix destination: '{args.output}'")
