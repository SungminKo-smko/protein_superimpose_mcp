"""
Protein Superimpose MCP Server.

단백질 구조 superimposition을 위한 MCP 서버.
FastMCP 패턴을 사용하여 CIF 구조 검사, 그룹별 정렬, 전체 정렬 도구를 제공.
Azure Files를 통한 파일 업로드/다운로드를 지원.
"""

import os
import base64
import shutil
from pathlib import Path

from mcp.server.fastmcp import FastMCP

from .core import (
    inspect_cif,
    superimpose_group as _superimpose_group,
    superimpose_all as _superimpose_all,
)

_is_container = os.environ.get("CONTAINER_APP_NAME") or os.environ.get("MCP_HOST")

# Azure Files 마운트 경로 (ACA에서는 /data, 로컬에서는 ./data)
DATA_DIR = Path(os.environ.get("DATA_DIR", "/data" if _is_container else "./data"))
UPLOAD_DIR = DATA_DIR / "upload"
OUTPUT_DIR = DATA_DIR / "output"

mcp = FastMCP(
    "protein-superimpose",
    host="0.0.0.0" if _is_container else "127.0.0.1",
    port=8000,
)


@mcp.tool()
def upload_file(filename: str, content_base64: str, subfolder: str = "") -> dict:
    """CIF 파일을 서버에 업로드합니다.

    파일 내용을 base64로 인코딩하여 전송합니다.
    업로드된 파일은 서버의 upload 디렉토리에 저장됩니다.

    Args:
        filename: 저장할 파일명 (예: structure.cif)
        content_base64: 파일 내용의 base64 인코딩 문자열
        subfolder: upload 디렉토리 내 하위 폴더 (선택)
    """
    target_dir = UPLOAD_DIR / subfolder if subfolder else UPLOAD_DIR
    target_dir.mkdir(parents=True, exist_ok=True)
    target_path = target_dir / filename

    try:
        content = base64.b64decode(content_base64)
        target_path.write_bytes(content)
        return {
            "status": "success",
            "path": str(target_path),
            "size_bytes": len(content),
            "filename": filename,
        }
    except Exception as e:
        return {"status": "error", "error": str(e)}


@mcp.tool()
def download_file(path: str) -> dict:
    """서버의 파일을 base64로 인코딩하여 반환합니다.

    주로 정렬 결과 파일을 다운로드할 때 사용합니다.

    Args:
        path: 다운로드할 파일의 서버 내 경로
    """
    file_path = Path(path)
    if not file_path.exists():
        return {"status": "error", "error": f"파일이 존재하지 않습니다: {path}"}

    try:
        content = file_path.read_bytes()
        return {
            "status": "success",
            "filename": file_path.name,
            "size_bytes": len(content),
            "content_base64": base64.b64encode(content).decode("utf-8"),
        }
    except Exception as e:
        return {"status": "error", "error": str(e)}


@mcp.tool()
def list_server_files(directory: str = "") -> dict:
    """서버의 데이터 디렉토리 내 파일 목록을 반환합니다.

    upload, output 등 서버 데이터 디렉토리의 파일을 조회합니다.

    Args:
        directory: 조회할 하위 경로. 빈 문자열이면 데이터 루트 전체를 조회합니다.
    """
    target = DATA_DIR / directory if directory else DATA_DIR
    if not target.exists():
        return {"directory": str(target), "exists": False, "files": []}

    files = []
    for p in sorted(target.rglob("*")):
        if p.is_file():
            files.append({
                "path": str(p),
                "name": p.name,
                "size_bytes": p.stat().st_size,
                "relative_path": str(p.relative_to(DATA_DIR)),
            })

    return {
        "directory": str(target),
        "data_root": str(DATA_DIR),
        "upload_dir": str(UPLOAD_DIR),
        "output_dir": str(OUTPUT_DIR),
        "count": len(files),
        "files": files,
    }


@mcp.tool()
def cleanup(directory: str = "output") -> dict:
    """서버의 데이터 디렉토리를 정리합니다.

    Args:
        directory: 정리할 디렉토리 ("upload", "output", 또는 "all")
    """
    removed = 0
    errors = []

    if directory == "all":
        targets = [UPLOAD_DIR, OUTPUT_DIR]
    elif directory == "upload":
        targets = [UPLOAD_DIR]
    elif directory == "output":
        targets = [OUTPUT_DIR]
    else:
        return {"status": "error", "error": f"지원하지 않는 디렉토리: {directory}. 'upload', 'output', 'all' 중 선택"}

    for t in targets:
        if t.exists():
            for p in t.iterdir():
                try:
                    if p.is_file():
                        p.unlink()
                    elif p.is_dir():
                        shutil.rmtree(p)
                    removed += 1
                except Exception as e:
                    errors.append(str(e))

    return {"status": "success", "removed": removed, "errors": errors}


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

    if "--transport" in sys.argv and "sse" in sys.argv:
        mcp.run(transport="sse")
    else:
        mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
