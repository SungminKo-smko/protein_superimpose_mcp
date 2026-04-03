"""
단백질 구조 Superimposition 핵심 로직.

두 가지 모드 지원:
1. superimpose_group: 디자인 ID별 그룹화 후 그룹 내 정렬
2. superimpose_all: 모든 CIF 파일을 단일 reference에 정렬
"""

import re
import shutil
import warnings
from pathlib import Path
from collections import defaultdict

import numpy as np
from Bio.PDB import MMCIFParser, Superimposer

from .cif_io import apply_transform_to_cif, parse_structure_with_retry

warnings.filterwarnings("ignore")


def get_ca_atoms(chain):
    """Chain에서 Cα 원자 리스트를 반환 (standard residue only)."""
    ca_atoms = []
    for residue in chain.get_residues():
        if residue.id[0] == " " and "CA" in residue:
            ca_atoms.append(residue["CA"])
    return ca_atoms


def get_ca_dict(chain):
    """
    Chain에서 {residue_seq_num: CA_atom} dict 반환.
    삽입 코드가 있는 잔기 및 비표준 잔기(HETATM)는 제외.
    """
    return {
        res.id[1]: res["CA"]
        for res in chain.get_residues()
        if res.id[0] == " " and "CA" in res
    }


def get_matched_ca_pairs(ref_ca_dict, mob_chain):
    """
    Reference CA dict와 mobile chain에서
    공통 잔기 번호(residue seq num) 기준 Cα 원자쌍을 반환.
    CDR 길이가 달라도 공통 프레임워크 잔기로 alignment 가능.
    """
    mob_ca_dict = get_ca_dict(mob_chain)
    common = sorted(set(ref_ca_dict.keys()) & set(mob_ca_dict.keys()))
    return [ref_ca_dict[i] for i in common], [mob_ca_dict[i] for i in common]


def superimpose_group(
    input_dir: str,
    output_dir: str,
    chain_id: str = "A",
    reference_model_idx: int = 0,
) -> dict:
    """
    input_dir 내 CIF 파일을 디자인 ID별로 그룹화하고,
    각 그룹에서 reference_model을 기준으로 나머지 모델을 지정 Chain의 Cα로 superimpose.
    원본 CIF의 모든 메타데이터(pLDDT 등)를 유지하며 결과를 output_dir에 저장.

    Returns
    -------
    dict
        success_count, error_count, group_count, messages
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    chain_id = chain_id.upper()

    # 파일 그룹화: design_id -> {model_idx: filepath}
    pattern = re.compile(r"^(.+)_model_(\d+)\.cif$")
    groups = defaultdict(dict)

    for fname in sorted(input_dir.glob("*.cif")):
        m = pattern.match(fname.name)
        if m:
            design_id = m.group(1)
            model_idx = int(m.group(2))
            groups[design_id][model_idx] = fname

    messages = []
    messages.append(f"총 디자인 그룹 수: {len(groups)}")
    messages.append(f"Superimpose 기준 Chain: {chain_id}")
    messages.append(f"Reference 모델 인덱스: {reference_model_idx}")

    parser = MMCIFParser(QUIET=True)
    sup = Superimposer()

    success_count = 0
    error_count = 0

    for design_id, model_files in sorted(groups.items()):
        if reference_model_idx not in model_files:
            messages.append(f"[WARN] {design_id}: reference model_{reference_model_idx} 없음, skip")
            error_count += 1
            continue

        # Reference 구조 로드
        ref_path = model_files[reference_model_idx]
        try:
            ref_struct = parse_structure_with_retry(parser, "ref", ref_path)
        except (TimeoutError, OSError) as e:
            messages.append(f"[ERROR] {design_id} reference: {e}, skip")
            error_count += 1
            continue
        ref_model = ref_struct[0]

        if chain_id not in ref_model:
            messages.append(f"[WARN] {design_id}: reference에 Chain {chain_id} 없음, skip")
            error_count += 1
            continue

        ref_ca = get_ca_atoms(ref_model[chain_id])
        if not ref_ca:
            messages.append(f"[WARN] {design_id}: reference Chain {chain_id}에 Cα 없음, skip")
            error_count += 1
            continue

        # Reference 모델은 변환 없이 원본 그대로 복사
        shutil.copy2(str(ref_path), str(output_dir / ref_path.name))

        # 나머지 모델 superimpose
        for idx, mob_path in sorted(model_files.items()):
            if idx == reference_model_idx:
                continue

            try:
                mob_struct = parse_structure_with_retry(parser, "mob", mob_path)
            except (TimeoutError, OSError) as e:
                messages.append(f"[ERROR] {design_id} model_{idx}: {e}, skip")
                error_count += 1
                continue
            mob_model = mob_struct[0]

            if chain_id not in mob_model:
                messages.append(f"[WARN] {design_id} model_{idx}: Chain {chain_id} 없음, skip")
                error_count += 1
                continue

            mob_ca = get_ca_atoms(mob_model[chain_id])

            if len(mob_ca) != len(ref_ca):
                messages.append(
                    f"[WARN] {design_id} model_{idx}: "
                    f"Cα 수 불일치 (ref={len(ref_ca)}, mob={len(mob_ca)}), skip"
                )
                error_count += 1
                continue

            # 지정 Chain Cα 기준으로 superimpose 계산
            sup.set_atoms(ref_ca, mob_ca)
            rot, tran = sup.rotran

            # 원본 CIF 메타데이터를 유지하면서 좌표만 변환하여 저장
            out_path = output_dir / mob_path.name
            try:
                apply_transform_to_cif(mob_path, out_path, rot, tran)
            except Exception as e:
                messages.append(f"[ERROR] {design_id} model_{idx} 저장 실패: {e}, skip")
                error_count += 1
                continue

            success_count += 1

    messages.append(f"완료: 성공 {success_count}개, 오류 {error_count}개")
    messages.append(f"결과 저장 위치: {output_dir}")

    return {
        "success_count": success_count,
        "error_count": error_count,
        "group_count": len(groups),
        "output_dir": str(output_dir),
        "messages": messages,
    }


def superimpose_all(
    input_root: str,
    output_root: str,
    chain_id: str = "A",
    reference_path: str = None,
) -> dict:
    """
    input_root 하위의 모든 CIF 파일을 단일 reference 기준으로 superimpose.

    Parameters
    ----------
    input_root     : CIF 파일을 포함하는 최상위 디렉토리 (재귀 탐색)
    output_root    : 결과 저장 디렉토리 (입력 구조 미러링)
    chain_id       : Superimpose 기준 chain ID (대문자로 자동 변환)
    reference_path : 기준 CIF 파일 경로. 미지정 시 input_root 내 첫 번째 파일 사용

    Returns
    -------
    dict
        success, skip, error, total, reference, messages
    """
    input_root = Path(input_root).resolve()
    output_root = Path(output_root).resolve()
    chain_id = chain_id.upper()

    messages = []

    # 파일 수집
    all_files = sorted(input_root.rglob("*.cif"))
    if not all_files:
        return {
            "success": 0,
            "skip": 0,
            "error": 0,
            "total": 0,
            "reference": None,
            "output_root": str(output_root),
            "messages": [f"[ERROR] CIF 파일을 찾을 수 없습니다: {input_root}"],
        }

    total = len(all_files)
    messages.append(f"총 CIF 파일 수: {total}")
    messages.append(f"Superimpose Chain: {chain_id}")

    # Reference 결정
    if reference_path:
        ref_path = Path(reference_path).resolve()
        if not ref_path.exists():
            raise FileNotFoundError(f"Reference 파일 없음: {ref_path}")
    else:
        ref_path = all_files[0]

    try:
        messages.append(f"Reference 구조: {ref_path.relative_to(input_root)}")
    except ValueError:
        messages.append(f"Reference 구조: {ref_path}")

    # Reference 파싱
    parser = MMCIFParser(QUIET=True)
    try:
        ref_struct = parser.get_structure("ref", str(ref_path))
    except Exception as e:
        raise RuntimeError(f"Reference 파싱 실패: {e}")

    ref_model = ref_struct[0]
    if chain_id not in ref_model:
        available = [c.id for c in ref_model.get_chains()]
        raise ValueError(
            f"Reference에 Chain '{chain_id}' 없음. "
            f"사용 가능한 chain: {available}"
        )

    ref_ca_dict = get_ca_dict(ref_model[chain_id])
    if not ref_ca_dict:
        raise ValueError(f"Reference Chain '{chain_id}'에 Cα 원자 없음")
    messages.append(f"Reference Cα 수: {len(ref_ca_dict)}")

    # Reference 복사
    ref_out = output_root / ref_path.relative_to(input_root)
    ref_out.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(str(ref_path), str(ref_out))

    # Superimpose
    sup = Superimposer()
    success = skip = error = 0

    for i, cif_path in enumerate(all_files, 1):
        if cif_path == ref_path:
            continue

        out_path = output_root / cif_path.relative_to(input_root)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            mob_struct = parser.get_structure("mob", str(cif_path))
            mob_model = mob_struct[0]

            if chain_id not in mob_model:
                messages.append(f"[WARN] Chain '{chain_id}' 없음: {cif_path.relative_to(input_root)}")
                skip += 1
                continue

            ref_ca, mob_ca = get_matched_ca_pairs(ref_ca_dict, mob_model[chain_id])

            if len(ref_ca) < 3:
                messages.append(
                    f"[WARN] 공통 Cα {len(ref_ca)}개 (최소 3 필요): "
                    f"{cif_path.relative_to(input_root)}"
                )
                skip += 1
                continue

            sup.set_atoms(ref_ca, mob_ca)
            rot, tran = sup.rotran

            apply_transform_to_cif(cif_path, out_path, rot, tran)
            success += 1

        except Exception as e:
            messages.append(f"[ERROR] {cif_path.relative_to(input_root)}: {e}")
            error += 1

    messages.append(f"완료: 성공 {success}개, 스킵 {skip}개, 오류 {error}개")
    messages.append(f"결과 저장 위치: {output_root}")

    return {
        "success": success,
        "skip": skip,
        "error": error,
        "total": total,
        "reference": str(ref_path),
        "output_root": str(output_root),
        "messages": messages,
    }
