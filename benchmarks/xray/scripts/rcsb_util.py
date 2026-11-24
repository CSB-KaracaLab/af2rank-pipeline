"""Small helpers for the public RCSB PDB APIs.

Uses the macOS system curl (/usr/bin/curl) because the Anaconda Python/curl CA
bundle on this machine cannot verify the RCSB certificate chain.
"""
import json, subprocess, time

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
CURL = "/usr/bin/curl"


def _curl(args, data=None, retries=5):
    for attempt in range(retries):
        cmd = [CURL, "-sS", "-m", "300", "-w", "\n%{http_code}"] + args
        p = subprocess.run(cmd, input=data, capture_output=True, text=True)
        if p.returncode == 0:
            body, _, code = p.stdout.rpartition("\n")
            code = int(code)
            if code == 204:          # search API: no hits
                return None
            if code == 200:
                return json.loads(body)
            if code in (400, 404):
                raise RuntimeError(f"HTTP {code}: {body[:500]}")
        time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"request failed after {retries} tries: {args[-1]} {p.stderr[:300]}")


def post_json(url, payload):
    return _curl(["-H", "Content-Type: application/json", "--data-binary", "@-", url],
                 data=json.dumps(payload))


def get_json(url):
    return _curl([url])


def search_all(query):
    """Run a search query with return_all_hits; return list of identifiers (or [])."""
    q = json.loads(json.dumps(query))
    q.setdefault("request_options", {})["return_all_hits"] = True
    res = post_json(SEARCH_URL, q)
    if not res:
        return []
    return [r["identifier"] for r in res.get("result_set", [])]


def search_count(query):
    q = json.loads(json.dumps(query))
    q["request_options"] = {"return_counts": True}
    res = post_json(SEARCH_URL, q)
    return 0 if not res else res["total_count"]


def graphql(query, variables=None):
    res = post_json(GRAPHQL_URL, {"query": query, "variables": variables or {}})
    if res.get("errors"):
        raise RuntimeError(json.dumps(res["errors"])[:1000])
    return res["data"]
