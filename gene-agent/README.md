# Gene Agent · BioReAct Demo

A bioinformatics gene query agent built as a learning project to explore
core AI agent engineering concepts: **ReAct**, **Skills**, **MCP**, and
**multi-provider LLM support**.

---

## What This Project Demonstrates

| Concept | Implementation |
|---|---|
| ReAct Loop | Reason → Act (tool call) → Observe → repeat, up to 5 rounds per skill |
| Skills | Behaviour cards as `.md` files with YAML frontmatter; body = system prompt |
| LLM Router | A small LLM call selects which skill(s) to activate per query |
| MCP | Tools exposed as a standalone stdio server, discovered dynamically by the agent |
| Multi-provider | One `AsyncOpenAI` client handles Local Ollama, VectorEngine, or any OpenAI-compatible API |

---

## Architecture

```
Browser
  │  POST /run  {query, provider, api_key, model_name, base_url}
  ▼
gene_agent_server.py        ← Flask HTTP server (port 5001)
  │
  ▼
agent.py  ·  GeneAgent
  ├─ _route()               ← LLM selects skills for this query
  └─ _run_skill()           ← ReAct loop per skill
       │  tool calls
       ▼
mcp_server.py               ← MCP stdio server
  ├─ search_gene_ncbi        → NCBI Entrez API
  └─ get_uniprot_function    → UniProt REST API
```

**Three layers, each with a single responsibility:**

- `gene_agent_server.py` — HTTP interface only; no business logic
- `agent.py` — Orchestration only; no API calls, no HTTP
- `mcp_server.py` — Tool execution only; no LLM, no agent logic

---

## Skills

Skills are Markdown files (`skills/`) with a YAML frontmatter header.
The frontmatter declares which tools the skill is allowed to use;
the Markdown body becomes the system prompt for the ReAct loop.

```
skills/
├── gene_genomics.md     # Identify gene location, chromosome, aliases
└── gene_proteomics.md   # Protein function, UniProt accession, disease associations
```

### `gene_genomics`
- **Tool**: `search_gene_ncbi`
- **Output**: Genomic Identity Card — Gene ID, chromosome location, aliases, NCBI summary

### `gene_proteomics`
- **Tool**: `get_uniprot_function`
- **Output**: Protein Function & Disease Report — function summary, disease list, UniProt link

The LLM router inspects both skill descriptions and decides which to run:

| Query | Skills activated |
|---|---|
| "Which chromosome is BRCA1 on?" | `gene_genomics` only |
| "What diseases are linked to TP53?" | `gene_proteomics` only |
| "Tell me everything about EGFR" | `gene_genomics` + `gene_proteomics` |

---

## Providers

The agent uses the **OpenAI-compatible chat completions API** for all providers.
`providers.py` is the single source of truth; `agent.py` never touches API keys or URLs.

| Provider | Base URL | Model | Key required |
|---|---|---|---|
| Local | `http://localhost:11434/v1` | `qwen3:30b-a3b` | No (Ollama) |
| Online | Configurable via UI | Configurable via UI | Yes |

### Adding a new provider

Edit `PROVIDER_CONFIGS` in `providers.py`:

```python
"my_provider": {
    "label": "My Provider",
    "base_url": "https://api.example.com/v1",
    "api_key": None,       # None = user supplies at runtime
    "model": "my-model",
},
```

---

## Quick Start

### 1. Prerequisites

- Python ≥ 3.13, [`uv`](https://docs.astral.sh/uv/)
- [Ollama](https://ollama.com) running locally with `qwen3:30b-a3b`
  ```bash
  ollama pull qwen3:30b-a3b
  ```

### 2. Install dependencies

```bash
cd gene-agent
uv sync
```

### 3. Start the server

```bash
uv run python gene_agent_server.py
```

Open **http://localhost:5001** in your browser.

### 4. CLI test (optional)

```bash
uv run python main.py
```

---

## Web UI Walkthrough

```
┌─────────────────────────────────┬──────────────────────────────────────┐
│  Provider                       │                                      │
│  [ Local ]  [ Online ]          │  ROUTE · Router selected 2 skill(s) │
│                                 │  ┌──────────────┐ ┌───────────────┐ │
│  Query Gene                     │  │gene_genomics │ │gene_proteomics│ │
│  ┌───────────────────────────┐  │  └──────────────┘ └───────────────┘ │
│  │ Query detailed info       │  │                                      │
│  │ about TP53 gene           │  │  SKILL · Running: gene_genomics      │
│  └───────────────────────────┘  │                                      │
│  [ ▶ Run Agent ]                │  ACT · Call tool → search_gene_ncbi  │
│                                 │    gene_name  "TP53"                 │
│  Quick Select                   │    organism   "human"                │
│  TP53  BRCA1  EGFR  PTEN ...   │                                      │
│                                 │  OBSERVE · search_gene_ncbi returned │
│  Available Tools                │  { gene_id: "7157", chromosome: "17" │
│  search_gene_ncbi               │    location: "17p13.1", ... }        │
│  get_uniprot_function           │                                      │
│                                 │  FINAL ANSWER · Report               │
└─────────────────────────────────┴──────────────────────────────────────┘
```

**Selecting Online provider** reveals three additional inputs:

| Field | Description |
|---|---|
| URL | API base URL (dropdown; default: VectorEngine) |
| API Key | Your API key for the selected service |
| Model | Dropdown populated automatically per URL |

**Supported online models** (via VectorEngine proxy):

```
deepseek-v3.2          deepseek-v3.2-thinking   deepseek-v3.2-fast
MiniMax-M2.5           kimi-k2.5
qwen3-coder-plus       qwen3.5-plus             qwen3.5-397b-a17b
qwen3-coder            qwen3-coder-480b-a35b-instruct
```

---

## Project Structure

```
gene-agent/
├── agent.py                  # GeneAgent: router + ReAct loop
├── providers.py              # Provider configs + get_client() factory
├── skill_loader.py           # Parses YAML frontmatter from skill files
├── mcp_server.py             # MCP stdio server (tools)
├── gene_agent_server.py      # Flask HTTP API
├── gene_agent_ui.html        # Web UI
├── main.py                   # CLI entry point
├── skills/
│   ├── gene_genomics.md
│   └── gene_proteomics.md
├── _archive/                 # Historical reference (V1 code)
│   ├── easy_function_calling_version/
│   └── without_mcp_tools_calling/
├── record_and_thinking/      # Observations on model behaviour
└── pyproject.toml
```

---

## Evolution History

This project was built incrementally as a learning exercise:

| Version | What was added |
|---|---|
| V1 | Simple function calling — monolithic agent, hard-coded tools |
| V2 | Skills system — behaviour in `.md` files, LLM router selects skills |
| V3 | MCP — tools extracted to a standalone stdio server |
| V4 | Multi-provider — `providers.py` + OpenAI SDK for all providers |

The V1 code is preserved in `_archive/easy_function_calling_version/`
for comparison.

---

## Environment Management

```bash
uv sync              # install / update all dependencies
uv add <package>     # add a dependency
uv remove <package>  # remove a dependency
uv run python <file> # run inside the managed virtualenv
```
