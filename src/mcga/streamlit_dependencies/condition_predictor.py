# condition_predictor.py

import requests

ASKCOS_ENDPOINT = "http://askcos.mit.edu/context"

class AskcosClient:
    """Client minimal pour l’API ASKCOS (gratuit, pas de clé)."""

    def __init__(self, endpoint: str = ASKCOS_ENDPOINT):
        self.endpoint = endpoint.rstrip("/")

    def predict_conditions(self,
                           reactants: list[str],
                           products: list[str],
                           agents: list[str] = None,
                           top_k: int = 3) -> list[dict]:
        """
        Envoie reactants/products/agents (SMILES) et renvoie les top_k
        prédictions sous forme de liste de dicts:
        {"solvent": "...", "catalyst": "...", ...}
        """
        payload = {
            "reactants": reactants,
            "products": products,
        }
        if agents:
            payload["agents"] = agents

        resp = requests.post(self.endpoint, json=payload, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        return data.get("predictions", [])[:top_k]
