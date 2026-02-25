"""
gene_agent_server.py
Wraps the gene agent as an HTTP service for the frontend interface.
Dependencies: pip install flask flask-cors
"""

import asyncio
import os

from flask import Flask, request, jsonify, send_file
from flask_cors import CORS

from agent import GeneAgent
from providers import PROVIDER_CONFIGS

app = Flask(__name__)
CORS(app)

_script_dir = os.path.dirname(os.path.abspath(__file__))
_agent = GeneAgent(
    mcp_server_script=os.path.join(_script_dir, "mcp_server.py"),
    skills_dir=os.path.join(_script_dir, "skills"),
)


@app.route("/run", methods=["POST"])
def run():
    data = request.get_json()
    query = data.get("query", "").strip()
    if not query:
        return jsonify({"error": "query is required"}), 400

    provider_name = data.get("provider", "local")
    api_key = data.get("api_key", None) or None
    model_name = data.get("model_name", None) or None
    base_url = data.get("base_url", None) or None

    if provider_name not in PROVIDER_CONFIGS:
        return jsonify({"error": f"Unknown provider: {provider_name}"}), 400

    try:
        answer, steps = asyncio.run(
            _agent.run(query, provider_name=provider_name, api_key=api_key, model_name=model_name, base_url=base_url, verbose=True)
        )
        serializable_steps = [{k: v for k, v in step.items()} for step in steps]
        used_model = model_name or PROVIDER_CONFIGS[provider_name]["model"] or "—"
        return jsonify({
            "answer": answer,
            "steps": serializable_steps,
            "provider": PROVIDER_CONFIGS[provider_name]["label"],
            "model": used_model,
        })
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


@app.route("/health", methods=["GET"])
def health():
    return jsonify({
        "status": "ok",
        "providers": {k: v["label"] for k, v in PROVIDER_CONFIGS.items()},
    })


@app.route("/")
def index():
    return send_file(os.path.join(_script_dir, "gene_agent_ui.html"))


if __name__ == "__main__":
    print("🧬 Gene Agent Server starting...")
    print(f"   Providers: {[v['label'] for v in PROVIDER_CONFIGS.values()]}")
    print("   API:    http://localhost:5001/run")
    print("   Health: http://localhost:5001/health")
    app.run(host="0.0.0.0", port=5001, debug=False)
