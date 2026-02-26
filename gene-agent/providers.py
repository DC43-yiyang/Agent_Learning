"""
providers.py — LLM Provider Configurations
==========================================
All providers use the OpenAI-compatible chat completions API format.
Local Ollama is accessed via its built-in OpenAI-compatible endpoint
(http://localhost:11434/v1), so the agent code is identical for all providers.

To add a new provider: add an entry to PROVIDER_CONFIGS.
"""

import os

from dotenv import load_dotenv
from openai import AsyncOpenAI

load_dotenv()

PROVIDER_CONFIGS = {
    "local": {
        "label": "Local (Ollama)",
        "base_url": "http://localhost:11434/v1",
        "api_key": "ollama",        # Ollama ignores the key but the SDK requires one
        "model": "qwen3:30b-a3b",
    },
    "online": {
        "label": "Online",
        "base_url": os.getenv("ONLINE_BASE_URL", "https://api.vectorengine.ai/v1"),
        "api_key": os.getenv("ONLINE_API_KEY"),
        "model": os.getenv("ONLINE_MODEL", "MiniMax-M2.5"),
    },
}


def get_client(
    provider_name: str,
    api_key: str | None = None,
    model_name: str | None = None,
    base_url: str | None = None,
) -> tuple[AsyncOpenAI, str]:
    """Return (AsyncOpenAI client, model_name) for the given provider.

    api_key, model_name, and base_url can be supplied at runtime (e.g. from
    the web UI); they fall back to the values in PROVIDER_CONFIGS if not provided.
    """
    if provider_name not in PROVIDER_CONFIGS:
        raise ValueError(
            f"Unknown provider '{provider_name}'. "
            f"Choose from: {list(PROVIDER_CONFIGS.keys())}"
        )

    cfg = PROVIDER_CONFIGS[provider_name]
    key = api_key or cfg["api_key"]
    model = model_name or cfg["model"]
    url = base_url or cfg["base_url"]

    if key is None:
        raise ValueError(f"Provider '{provider_name}' requires an API key.")
    if model is None:
        raise ValueError(f"Provider '{provider_name}' requires a model name.")

    client = AsyncOpenAI(base_url=url, api_key=key)
    return client, model
