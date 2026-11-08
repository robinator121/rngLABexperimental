#!/usr/bin/env python3
# RNGLab v0.18 (written by Claude Code, 2026-11-07). Output for the Stage B follow-up: docs/v018/evidence/agy_literature_check.json
"""Check AI-proposed references against Crossref: does the DOI exist, and do title, year and first author match?

Usage: python3 verify_dois.py <file containing a JSON array> > verified.json
"""
import difflib
import json
import re
import sys
import time
import urllib.parse
import urllib.request


def crossref(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    req = urllib.request.Request(url, headers={"User-Agent": "rnglab-citation-check/1.0 (mailto:none)"})
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            return json.load(r)["message"]
    except Exception as e:  # noqa: BLE001
        return {"error": str(e)}


def main():
    text = open(sys.argv[1]).read()
    items = json.loads(re.search(r"\[.*\]", text, re.S).group(0))
    out = []
    for it in items:
        m = crossref(it.get("doi", ""))
        if "error" in m:
            out.append({**it, "verdict": "DOI_NOT_FOUND", "detail": m["error"]})
            continue
        title = (m.get("title") or [""])[0]
        year = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
        first = ((m.get("author") or [{}])[0].get("family") or "")
        sim = difflib.SequenceMatcher(None, title.lower(), it.get("title", "").lower()).ratio()
        ok = sim >= 0.85 and (year == it.get("year")) and first.lower() in it.get("authors", "").lower()
        out.append({**it, "verdict": "MATCH" if ok else "MISMATCH", "crossref_title": title, "crossref_year": year,
                    "crossref_first_author": first, "crossref_venue": (m.get("container-title") or [""])[0],
                    "title_similarity": round(sim, 2)})
        time.sleep(0.3)
    json.dump(out, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
