from typing import List as _List

from python_nedrex import config as _config
from python_nedrex.common import check_response as _check_response
from python_nedrex.common import check_status_factory as _check_status_factory
from python_nedrex.common import http as _http


def domino_submit(seeds: _List[str], network: str = "DEFAULT") -> str:
    url = f"{_config.url_base}/domino/submit"
    body = {"seeds": seeds, "network": network}

    resp = _http.post(url, json=body, headers={"x-api-key": _config.api_key})
    result: str = _check_response(resp)
    return result


check_domino_status = _check_status_factory("/domino/status")
