"""
MCP Server — Gene Tools
========================
This file turns the gene query functions into an MCP-compliant tool server.

Key MCP concepts demonstrated here:
  - Server:      the process that exposes tools
  - @list_tools: tells any client "here is what I can do"
  - @call_tool:  actually executes a tool when a client asks
  - stdio:       the transport layer (client and server talk via stdin/stdout)

Once running, ANY MCP-compatible client can use these tools:
  - Your own agent (agent.py, Phase 2)
  - Claude Desktop
  - Cursor
  - Any future agent you build
"""

import asyncio
import json

import requests
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp import types


# ─────────────────────────────────────────────
# 1. Tool implementations
#    Identical logic to tools.py — same API calls.
#    The difference: these are now exposed via a standard protocol,
#    not imported as Python functions.
# ─────────────────────────────────────────────

def search_gene_ncbi(gene_name: str, organism: str = "human") -> dict:
    """Query basic gene information via NCBI Entrez API."""
    try:
        search_resp = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params={
                "db": "gene",
                "term": f"{gene_name}[Gene Name] AND {organism}[Organism]",
                "retmode": "json",
                "retmax": 1,
            },
            timeout=10,
        )
        id_list = search_resp.json().get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return {"error": f"Gene not found: {gene_name}"}

        gene_id = id_list[0]
        summary_resp = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params={"db": "gene", "id": gene_id, "retmode": "json"},
            timeout=10,
        )
        result = summary_resp.json().get("result", {}).get(gene_id, {})
        summary = result.get("summary", "N/A")

        return {
            "gene_id": gene_id,
            "name": result.get("name", "N/A"),
            "full_name": result.get("description", "N/A"),
            "organism": result.get("organism", {}).get("scientificname", "N/A"),
            "chromosome": result.get("chromosome", "N/A"),
            "location": result.get("maplocation", "N/A"),
            "summary": summary[:500] + "..." if len(summary) > 500 else summary,
            "aliases": result.get("otheraliases", "N/A"),
        }
    except Exception as e:
        return {"error": f"NCBI API call failed: {str(e)}"}


def get_uniprot_function(gene_name: str, organism: str = "human") -> dict:
    """Get protein function description via UniProt API."""
    try:
        organism_map = {"human": "Homo sapiens", "mouse": "Mus musculus"}
        resp = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={
                "query": f"gene:{gene_name} AND organism_name:{organism_map.get(organism, organism)} AND reviewed:true",
                "fields": "gene_names,protein_name,cc_function,cc_disease,go_p",
                "format": "json",
                "size": 1,
            },
            timeout=10,
        )
        results = resp.json().get("results", [])
        if not results:
            return {"error": f"UniProt found no reviewed entry: {gene_name}"}

        entry = results[0]
        function_text, disease_text = "", ""
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    function_text = texts[0].get("value", "")[:600]
            if comment.get("commentType") == "DISEASE":
                disease_name = comment.get("disease", {}).get("diseaseName", "")
                if disease_name:
                    disease_text += disease_name + "; "

        protein_name = (
            entry.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", "N/A")
        )
        return {
            "uniprot_id": entry.get("primaryAccession", "N/A"),
            "protein_name": protein_name,
            "function": function_text or "N/A",
            "associated_diseases": disease_text.rstrip("; ") or "N/A",
        }
    except Exception as e:
        return {"error": f"UniProt API call failed: {str(e)}"}


# ─────────────────────────────────────────────
# 2. MCP Server setup
#    app = the server instance, identified by a name.
#    Any client that connects will see this name.
# ─────────────────────────────────────────────

app = Server("gene-tools")


# ─────────────────────────────────────────────
# 3. list_tools handler
#    When a client connects and asks "what can you do?",
#    MCP calls this function and sends the result back.
#    This is how clients discover tools dynamically —
#    no hardcoding on the client side.
# ─────────────────────────────────────────────

@app.list_tools()
async def list_tools() -> list[types.Tool]:
    return [
        types.Tool(
            name="search_gene_ncbi",
            description=(
                "Query basic gene information from the NCBI database, including "
                "gene ID, chromosome location, function summary, and aliases."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "gene_name": {
                        "type": "string",
                        "description": "Gene symbol, e.g. TP53, BRCA1",
                    },
                    "organism": {
                        "type": "string",
                        "description": "Organism name, default: human",
                        "default": "human",
                    },
                },
                "required": ["gene_name"],
            },
        ),
        types.Tool(
            name="get_uniprot_function",
            description=(
                "Query the UniProt database for detailed protein function "
                "description and associated diseases for a given gene."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "gene_name": {
                        "type": "string",
                        "description": "Gene symbol",
                    },
                    "organism": {
                        "type": "string",
                        "description": "Organism name, default: human",
                        "default": "human",
                    },
                },
                "required": ["gene_name"],
            },
        ),
    ]


# ─────────────────────────────────────────────
# 4. call_tool handler
#    When a client says "call search_gene_ncbi with these args",
#    MCP routes it here. We dispatch to the real function,
#    then return the result as a TextContent (MCP's standard result type).
#
#    Note: MCP results are always a list of Content objects.
#    TextContent is the most common — it just wraps a string.
# ─────────────────────────────────────────────

@app.call_tool()
async def call_tool(name: str, arguments: dict) -> list[types.TextContent]:
    dispatch = {
        "search_gene_ncbi": search_gene_ncbi,
        "get_uniprot_function": get_uniprot_function,
    }

    if name not in dispatch:
        return [types.TextContent(
            type="text",
            text=json.dumps({"error": f"Unknown tool: {name}"}),
        )]

    result = dispatch[name](**arguments)
    return [types.TextContent(type="text", text=json.dumps(result, ensure_ascii=False))]


# ─────────────────────────────────────────────
# 5. Entry point
#    stdio_server() sets up the stdin/stdout pipes.
#    app.run() starts the event loop and waits for client connections.
#    The server stays alive until the client disconnects.
# ─────────────────────────────────────────────

async def main():
    async with stdio_server() as (read_stream, write_stream):
        await app.run(
            read_stream,
            write_stream,
            app.create_initialization_options(),
        )


if __name__ == "__main__":
    asyncio.run(main())
