"""
Protein Superimpose MCP Server.

단백질 구조 superimposition을 위한 MCP 서버.
FastMCP 패턴을 사용하여 CIF 구조 검사, 그룹별 정렬, 전체 정렬 도구를 제공.
"""

from pathlib import Path

from mcp.server.fastmcp import FastMCP

from .core import (
    inspect_cif,
    superimpose_group as _superimpose_group,
    superimpose_all as _superimpose_all,
)

mcp = FastMCP(
    "protein-superimpose",
    description="단백질 구조 superimposition MCP 서버 - CIF 파일의 Cα 원자 기반 구조 정렬",
)


@mcp.tool()
def inspect_structure(path: str) -> dict:
    """CIF 파일의 구조 정보를 검사합니다.

    체인 목록, 잔기 수, Cα 원자 수, 전체 원자 수, 유효성 등을 반환합니다.

    Args:
        path: CIF 파일 경로
    """
    return inspect_cif(path)


@mcp.tool()
def superimpose_group(
    input_dir: str,
    output_dir: str,
    chain: str = "A",
    reference_model: int = 0,
) -> dict:
    """디자인 ID별로 CIF 파일을 그룹화하여 superimpose합니다.

    파일명 패턴 '{design_id}_model_{N}.cif'를 기준으로 그룹화하고,
    각 그룹 내에서 reference 모델을 기준으로 지정 Chain의 Cα 원자로 정렬합니다.
    원본 CIF 메타데이터(pLDDT 등)는 완전히 보존됩니다.

    Args:
        input_dir: 입력 CIF 파일 디렉토리 경로
        output_dir: 출력 CIF 파일 디렉토리 경로 (없으면 자동 생성)
        chain: Superimpose 기준 chain ID (예: A, B). 기본값 "A"
        reference_model: 기준 모델 인덱스. 기본값 0
    """
    return _superimpose_group(
        input_dir=input_dir,
        output_dir=output_dir,
        chain_id=chain,
        reference_model_idx=reference_model,
    )


@mcp.tool()
def superimpose_all(
    input_root: str,
    output_root: str,
    chain: str = "A",
    reference: str | None = None,
) -> dict:
    """디렉토리 트리 내 모든 CIF 파일을 단일 reference 기준으로 superimpose합니다.

    하위 디렉토리를 재귀적으로 탐색하여 모든 CIF 파일을 수집하고,
    지정한 chain의 Cα 원자 기준으로 정렬합니다.
    잔기 번호 기준 공통 Cα만 사용하므로 길이가 다른 설계도 허용됩니다.
    입력 디렉토리 구조를 유지하며 결과를 저장합니다.

    Args:
        input_root: CIF 파일이 있는 최상위 디렉토리 (하위 폴더 포함 재귀 탐색)
        output_root: 결과를 저장할 최상위 디렉토리 (입력 디렉토리 구조 미러링)
        chain: Superimpose 기준 chain ID (예: A, B). 기본값 "A"
        reference: 기준 CIF 파일 경로. 미지정 시 알파벳 순 첫 번째 파일 자동 선택
    """
    return _superimpose_all(
        input_root=input_root,
        output_root=output_root,
        chain_id=chain,
        reference_path=reference,
    )


@mcp.tool()
def list_cif_files(directory: str) -> dict:
    """디렉토리 내 CIF 파일 목록을 반환합니다.

    지정된 디렉토리를 재귀적으로 탐색하여 모든 .cif 파일의 경로와 파일 크기를 반환합니다.

    Args:
        directory: 탐색할 디렉토리 경로
    """
    dir_path = Path(directory).resolve()
    if not dir_path.exists():
        return {"error": f"디렉토리가 존재하지 않습니다: {directory}", "files": []}
    if not dir_path.is_dir():
        return {"error": f"디렉토리가 아닙니다: {directory}", "files": []}

    files = []
    for cif_path in sorted(dir_path.rglob("*.cif")):
        files.append({
            "path": str(cif_path),
            "name": cif_path.name,
            "size_bytes": cif_path.stat().st_size,
            "relative_path": str(cif_path.relative_to(dir_path)),
        })

    return {
        "directory": str(dir_path),
        "count": len(files),
        "files": files,
    }


def main():
    """MCP 서버 실행."""
    import sys

    transport = "sse" if "--transport" in sys.argv and "sse" in sys.argv else "stdio"
    mcp.run(transport=transport)


if __name__ == "__main__":
    main()
