FROM python:3.12-slim AS base

WORKDIR /app

COPY pyproject.toml .
COPY src/ src/

RUN pip install --no-cache-dir .

ENV MCP_HOST=0.0.0.0

EXPOSE 8000

CMD ["protein-superimpose-mcp", "--transport", "sse"]
