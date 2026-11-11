"""Crossref lookups that complement scripts/v018/verify_dois.py: metadata of a DOI, or a bibliographic query.

    python docs/v019/evidence/crossref_lookup.py doi:<DOI> ... query:"<free text>" ... > crossref_lookup.json

Records what Crossref returns (title, authors, venue, volume, pages, year). It reads no full text.
"""
import json
import sys
import urllib.parse
import urllib.request


def get(url):
    request = urllib.request.Request(url, headers={"User-Agent": "rnglab-citation-check/1.0 (mailto:none)"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)["message"]


def brief(item):
    return {"doi": item.get("DOI"), "title": (item.get("title") or [""])[0],
            "authors": [a.get("family") for a in item.get("author", [])],
            "venue": (item.get("container-title") or [""])[0], "volume": item.get("volume"), "issue": item.get("issue"),
            "page": item.get("page"), "year": (item.get("issued", {}).get("date-parts") or [[None]])[0][0],
            "published_print": (item.get("published-print", {}).get("date-parts") or [[None]])[0][0]}


def main():
    out = []
    arguments = []
    for argument in sys.argv[1:]:  # "@file" reads one lookup per line
        arguments += [line.strip() for line in open(argument[1:]) if line.strip()] if argument.startswith("@") else [argument]
    for argument in arguments:
        kind, _, value = argument.partition(":")
        try:
            if kind == "doi":
                out.append({"lookup": argument, "result": brief(get("https://api.crossref.org/works/" + urllib.parse.quote(value)))})
            else:
                items = get("https://api.crossref.org/works?rows=4&query.bibliographic=" + urllib.parse.quote(value))["items"]
                out.append({"lookup": argument, "result": [brief(i) for i in items]})
        except Exception as error:  # noqa: BLE001
            out.append({"lookup": argument, "error": str(error)})
    json.dump(out, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
