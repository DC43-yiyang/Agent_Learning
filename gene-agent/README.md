# Gene Agent · BioReAct Demo

A bioinformatics gene query agent implementing **ReAct** (Reasoning + Acting) structure with **simple function calling**.

## Features

- **ReAct Pattern**: Reason → Act (call tools) → Observe → Reason → ... until final answer
- **Function Calling**: Two tools exposed via JSON schema to the LLM:
  - `search_gene_ncbi`: Query gene info from NCBI Entrez API
  - `get_uniprot_function`: Query protein function from UniProt API
- **Web UI**: Flask backend + HTML frontend to visualize the reasoning steps

## Requirements

- Python >= 3.13
- [Ollama](https://ollama.com) with `qwen3:30b-a3b` model
- Dependencies: flask, flask-cors, ollama, requests

## Quick Start

```bash
# Install dependencies
uv sync   # or: pip install -e .

# Start backend server
python gene_agent_server.py
```

Then open http://localhost:5001 in your browser.

## Environment Management with uv

Use `uv` only for creating and maintaining the local environment (you can manage versions yourself):

```bash
# (Re)create / update the virtualenv from pyproject.toml / uv.lock
uv sync

# Add or update a dependency (example)
uv add flask-cors

# Remove a dependency (example)
uv remove flask-cors
```

To run code inside this environment you can either:

- `uv run python gene_agent_skills.py`  
  or  
- `.venv/bin/python gene_agent_skills.py`

## Project Structure

```
gene-agent/
├── gene_agent.py         # ReAct agent + function calling
├── gene_agent_server.py  # Flask HTTP API
├── gene_agent_ui.html    # Web UI
├── main.py
└── pyproject.toml
```
