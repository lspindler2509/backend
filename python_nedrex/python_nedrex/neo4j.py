import json as _json
from typing import Any as _Any
from typing import Dict as _Dict
from typing import Generator as _Generator

from requests.exceptions import ChunkedEncodingError  # type: ignore

from python_nedrex import config as _config
from python_nedrex.common import http as _http
from python_nedrex.exceptions import NeDRexError


def neo4j_query(query: str) -> _Generator[_Dict[str, _Any], None, None]:
    url = f"{_config.url_base}/neo4j/query"
    params = {"query": query}

    resp = _http.get(url, params=params, headers={"x-api-key": _config.api_key}, stream=True)
    if resp.status_code != 200:
        raise NeDRexError("Querying Neo4j returned a non-200 status code.")

    try:
        for line in resp.iter_lines():
            data = _json.loads(line.decode())
            for item in data:
                yield item
    except ChunkedEncodingError as exc:
        raise NeDRexError("cypher query could not be executed") from exc
