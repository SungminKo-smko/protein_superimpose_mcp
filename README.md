# Protein Superimpose MCP Server

단백질 구조 superimposition을 위한 MCP(Model Context Protocol) 서버.
CIF 파일의 Cα 원자를 기반으로 단백질 구조를 정렬합니다.

## 기능

- **inspect_structure**: CIF 파일의 구조 정보 검사 (체인, 잔기, Cα 원자 수)
- **superimpose_group**: 디자인 ID별 그룹화 후 그룹 내 구조 정렬
- **superimpose_all**: 디렉토리 트리 내 모든 CIF 파일을 단일 reference 기준으로 정렬
- **list_cif_files**: 디렉토리 내 CIF 파일 목록 조회

## 설치

```bash
# uv 사용 (권장)
cd /path/to/protein_superimpose_mcp
uv sync

# pip 사용
pip install -e .
```

## Claude Desktop 설정

`claude_desktop_config.json`에 다음을 추가:

```json
{
  "mcpServers": {
    "protein-superimpose": {
      "command": "uv",
      "args": ["run", "--directory", "/path/to/protein_superimpose_mcp", "protein-superimpose-mcp"]
    }
  }
}
```

## 사용 예시

Claude Desktop에서 다음과 같이 사용할 수 있습니다:

- "이 CIF 파일의 구조 정보를 확인해줘" → `inspect_structure` 호출
- "이 폴더의 모든 모델을 Chain A 기준으로 정렬해줘" → `superimpose_group` 호출
- "하위 폴더 전체의 CIF 파일을 하나의 reference에 맞춰 정렬해줘" → `superimpose_all` 호출

## 의존성

- [mcp](https://github.com/modelcontextprotocol/python-sdk) - Model Context Protocol SDK
- [BioPython](https://biopython.org/) - 구조 파싱 및 Superimposer
- [gemmi](https://gemmi.readthedocs.io/) - CIF 읽기/쓰기 (메타데이터 보존)
- [NumPy](https://numpy.org/) - 좌표 변환 연산
