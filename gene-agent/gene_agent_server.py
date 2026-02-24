"""
gene_agent_server.py
Wraps gene_agent.py as HTTP service for frontend interface
Dependencies: pip install flask flask-cors
"""

from flask import Flask, request, jsonify
from flask_cors import CORS
from gene_agent import run_agent

app = Flask(__name__)
CORS(app)  # Allow local frontend cross-origin requests

@app.route("/run", methods=["POST"])
def run():
    data = request.get_json()
    query = data.get("query", "").strip()
    if not query:
        return jsonify({"error": "query is required"}), 400
    
    try:
        answer, steps = run_agent(query, verbose=True)
        # Ensure all steps are JSON serializable
        serializable_steps = []
        for step in steps:
            serializable_steps.append({k: v for k, v in step.items()})
        return jsonify({"answer": answer, "steps": serializable_steps})
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route("/health", methods=["GET"])
def health():
    return jsonify({"status": "ok", "model": "qwen3:30b-a3b"})

from flask import send_file
import os

@app.route("/")
def index():
    return send_file(os.path.join(os.path.dirname(__file__), "gene_agent_ui.html"))

if __name__ == "__main__":
    print("🧬 Gene Agent Server starting...")
    print("   API:  http://localhost:5001/run")
    print("   Health: http://localhost:5001/health")
    app.run(host="0.0.0.0", port=5001, debug=False)
