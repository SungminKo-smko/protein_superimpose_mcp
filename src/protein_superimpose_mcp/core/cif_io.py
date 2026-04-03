"""
CIF 파일 읽기/쓰기 및 좌표 변환 유틸리티.
gemmi를 사용하여 원본 CIF 메타데이터(pLDDT, _ma_qa_metric_local 등)를 완전 보존.
"""

import time
from pathlib import Path

import numpy as np
import gemmi
from Bio.PDB import MMCIFParser


def parse_structure_with_retry(parser, name, filepath, retries=5, delay=10):
    """타임아웃 시 재시도하며 CIF 구조를 파싱."""
    for attempt in range(retries):
        try:
            return parser.get_structure(name, str(filepath))
        except (TimeoutError, OSError) as e:
            if attempt < retries - 1:
                print(f"    [재시도 {attempt+1}/{retries}] {Path(filepath).name}: {e}")
                time.sleep(delay)
            else:
                raise


def apply_transform_to_cif(input_path, output_path, rot, tran):
    """
    원본 CIF 파일의 모든 메타데이터(pLDDT, _ma_qa_metric_local 등)를 유지하면서
    rotation/translation을 좌표(_atom_site.Cartn_x/y/z)에만 적용하여 저장.

    변환식: new_coord = old_coord @ rot + tran  (BioPython Superimposer 규약)
    """
    doc = gemmi.cif.read(str(input_path))
    block = doc.sole_block()

    table = block.find("_atom_site.", ["Cartn_x", "Cartn_y", "Cartn_z"])
    for row in table:
        coord = np.array([float(row[0]), float(row[1]), float(row[2])])
        new_coord = coord @ rot + tran
        row[0] = f"{new_coord[0]:.3f}"
        row[1] = f"{new_coord[1]:.3f}"
        row[2] = f"{new_coord[2]:.3f}"

    doc.write_file(str(output_path))


def inspect_cif(filepath):
    """
    CIF 파일의 구조 정보를 반환.

    Returns
    -------
    dict
        chains: 체인 ID 목록
        residue_counts: {chain_id: 잔기 수}
        ca_counts: {chain_id: Cα 원자 수}
        total_atoms: 전체 원자 수
        valid: 파싱 성공 여부
        error: 에러 메시지 (파싱 실패 시)
    """
    filepath = Path(filepath)
    result = {
        "file": str(filepath),
        "filename": filepath.name,
        "chains": [],
        "residue_counts": {},
        "ca_counts": {},
        "total_atoms": 0,
        "valid": False,
        "error": None,
    }

    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("query", str(filepath))
        model = structure[0]

        chains = []
        residue_counts = {}
        ca_counts = {}
        total_atoms = 0

        for chain in model.get_chains():
            cid = chain.id
            chains.append(cid)

            residues = [r for r in chain.get_residues() if r.id[0] == " "]
            residue_counts[cid] = len(residues)

            ca_count = sum(1 for r in residues if "CA" in r)
            ca_counts[cid] = ca_count

            total_atoms += sum(len(list(r.get_atoms())) for r in chain.get_residues())

        result["chains"] = chains
        result["residue_counts"] = residue_counts
        result["ca_counts"] = ca_counts
        result["total_atoms"] = total_atoms
        result["valid"] = True

    except Exception as e:
        result["error"] = str(e)

    return result
