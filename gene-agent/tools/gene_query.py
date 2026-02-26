"""
tools/gene_query.py
-------------------
Gene and protein query functions via NCBI Entrez and UniProt APIs.
"""

import requests


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
