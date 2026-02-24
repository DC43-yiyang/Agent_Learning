"""
Bioinformatics Gene Query Agent
================================
First Agent Demo: ReAct structure + simple function calling
Dependencies: pip install ollama requests
"""

import json
import requests
import ollama

# ─────────────────────────────────────────────
# 1. Define Tools
#    Each tool = a Python function + a JSON schema description
# ─────────────────────────────────────────────

def search_gene_ncbi(gene_name: str, organism: str = "human") -> dict:
    """Query gene basic information via NCBI Entrez API"""
    try:
        # Step 1: Search to get Gene ID
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        search_params = {
            "db": "gene",
            "term": f"{gene_name}[Gene Name] AND {organism}[Organism]",
            "retmode": "json",
            "retmax": 1
        }
        search_resp = requests.get(search_url, params=search_params, timeout=10)
        search_data = search_resp.json()
        
        id_list = search_data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return {"error": f"Gene not found: {gene_name}"}
        
        gene_id = id_list[0]
        
        # Step 2: Get detailed info using Gene ID
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        summary_params = {
            "db": "gene",
            "id": gene_id,
            "retmode": "json"
        }
        summary_resp = requests.get(summary_url, params=summary_params, timeout=10)
        summary_data = summary_resp.json()
        
        result = summary_data.get("result", {}).get(gene_id, {})
        
        return {
            "gene_id": gene_id,
            "name": result.get("name", "N/A"),
            "full_name": result.get("description", "N/A"),
            "organism": result.get("organism", {}).get("scientificname", "N/A"),
            "chromosome": result.get("chromosome", "N/A"),
            "location": result.get("maplocation", "N/A"),
            "summary": result.get("summary", "N/A")[:500] + "..." if len(result.get("summary", "")) > 500 else result.get("summary", "N/A"),
            "aliases": result.get("otheraliases", "N/A")
        }
    except Exception as e:
        return {"error": f"NCBI API call failed: {str(e)}"}


def get_uniprot_function(gene_name: str, organism: str = "human") -> dict:
    """Get protein function description via UniProt API"""
    try:
        organism_map = {"human": "Homo sapiens", "mouse": "Mus musculus"}
        org_name = organism_map.get(organism, organism)
        
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            "query": f"gene:{gene_name} AND organism_name:{org_name} AND reviewed:true",
            "fields": "gene_names,protein_name,cc_function,cc_disease,go_p",
            "format": "json",
            "size": 1
        }
        resp = requests.get(url, params=params, timeout=10)
        data = resp.json()
        
        results = data.get("results", [])
        if not results:
            return {"error": f"UniProt found no reviewed entry: {gene_name}"}
        
        entry = results[0]
        
        # Extract function description
        comments = entry.get("comments", [])
        function_text = ""
        disease_text = ""
        for comment in comments:
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    function_text = texts[0].get("value", "")[:600]
            if comment.get("commentType") == "DISEASE":
                disease_name = comment.get("disease", {}).get("diseaseName", "")
                if disease_name:
                    disease_text += disease_name + "; "
        
        # Extract protein name
        protein_name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "N/A")
        
        # UniProt accession
        accession = entry.get("primaryAccession", "N/A")
        
        return {
            "uniprot_id": accession,
            "protein_name": protein_name,
            "function": function_text if function_text else "N/A",
            "associated_diseases": disease_text.rstrip("; ") if disease_text else "N/A",
        }
    except Exception as e:
        return {"error": f"UniProt API call failed: {str(e)}"}


# ─────────────────────────────────────────────
# 2. JSON Schema for Tools
#    Metadata for the model to know "what tools are available"
# ─────────────────────────────────────────────

TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "search_gene_ncbi",
            "description": "Query gene basic information from NCBI database, including gene ID, chromosome location, function summary, aliases, etc.",
            "parameters": {
                "type": "object",
                "properties": {
                    "gene_name": {
                        "type": "string",
                        "description": "Gene name, e.g. TP53, BRCA1, EGFR"
                    },
                    "organism": {
                        "type": "string",
                        "description": "Organism, default human, can also be mouse, etc.",
                        "default": "human"
                    }
                },
                "required": ["gene_name"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "get_uniprot_function",
            "description": "Get detailed protein function description and associated diseases from UniProt database.",
            "parameters": {
                "type": "object",
                "properties": {
                    "gene_name": {
                        "type": "string",
                        "description": "Gene name"
                    },
                    "organism": {
                        "type": "string",
                        "description": "Organism, default human",
                        "default": "human"
                    }
                },
                "required": ["gene_name"]
            }
        }
    }
]

# Tool name to function mapping
TOOL_REGISTRY = {
    "search_gene_ncbi": search_gene_ncbi,
    "get_uniprot_function": get_uniprot_function
}


# ─────────────────────────────────────────────
# 3. Agent Core Loop (ReAct Pattern)
#    Reason -> Act -> Observe -> Reason -> ...
# ─────────────────────────────────────────────

SYSTEM_PROMPT = """You are a professional bioinformatics assistant skilled at querying and analyzing gene information.

When users ask questions about genes, you should:
1. First use search_gene_ncbi to get basic gene information
2. Then use get_uniprot_function to get protein function description
3. Synthesize both sources and provide a comprehensive, professional analysis report

The report should include: basic info, chromosome location, protein function, related diseases (if any), and research significance.
Respond in English with a professional yet accessible style."""


def run_agent(user_query: str, verbose: bool = True) -> str:
    """
    Run Agent main loop

    This is the core of ReAct:
    - Model reasons (Reason)
    - Decides to call tools (Act)
    - Gets tool results (Observe)
    - Continues reasoning until final answer
    """
    
    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": user_query}
    ]
    
    steps = []  # Record each step for visualization
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"🧬 User query: {user_query}")
        print(f"{'='*60}")
    
    # ReAct loop, max 5 iterations to prevent infinite loop
    for iteration in range(5):
        if verbose:
            print(f"\n[Round {iteration + 1}] Model reasoning...")
        
        # Call model
        response = ollama.chat(
            model="qwen3:30b-a3b",
            messages=messages,
            tools=TOOLS,
            options={
                "temperature": 0.6,
                "num_ctx": 8192
            }
        )
        
        message = response.message
        
        # Check for tool call
        if message.tool_calls:
            for tool_call in message.tool_calls:
                tool_name = tool_call.function.name
                tool_args = tool_call.function.arguments
                if isinstance(tool_args, str):
                    tool_args = json.loads(tool_args) if tool_args else {}
                
                if verbose:
                    print(f"\n🔧 Calling tool: {tool_name}")
                    print(f"   Args: {json.dumps(tool_args, ensure_ascii=False)}")
                
                steps.append({
                    "type": "tool_call",
                    "tool": tool_name,
                    "args": tool_args
                })
                
                # Execute tool
                if tool_name in TOOL_REGISTRY:
                    tool_result = TOOL_REGISTRY[tool_name](**tool_args)
                else:
                    tool_result = {"error": f"Unknown tool: {tool_name}"}
                
                if verbose:
                    print(f"   Result: {json.dumps(tool_result, ensure_ascii=False, indent=2)[:300]}...")
                
                steps.append({
                    "type": "tool_result",
                    "tool": tool_name,
                    "result": tool_result
                })
                
                # Add tool result back to messages
                messages.append({"role": "assistant", "content": "", "tool_calls": message.tool_calls})
                messages.append({
                    "role": "tool",
                    "content": json.dumps(tool_result, ensure_ascii=False)
                })
        
        else:
            # No tool call, model gave final answer
            final_answer = message.content
            
            if verbose:
                print(f"\n{'='*60}")
                print("📋 Final report:")
                print(f"{'='*60}")
                print(final_answer)
            
            steps.append({
                "type": "final_answer",
                "content": final_answer
            })
            
            return final_answer, steps
    
    return "Max iterations reached, query incomplete.", steps


# ─────────────────────────────────────────────
# 4. Main entry point
# ─────────────────────────────────────────────

if __name__ == "__main__":
    # Test with classic cancer-related genes
    test_queries = [
        "Query detailed information about TP53 gene",
        # "Look up BRCA1 gene",
        # "What is the role of EGFR gene in cancer?"
    ]
    
    for query in test_queries:
        answer, steps = run_agent(query, verbose=True)
        print(f"\n✅ Done, {len([s for s in steps if s['type'] == 'tool_call'])} tool calls executed")
