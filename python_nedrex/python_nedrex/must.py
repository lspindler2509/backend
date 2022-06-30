from typing import List as _List

from python_nedrex import config as _config
from python_nedrex.common import check_response as _check_response
from python_nedrex.common import check_status_factory as _check_status_factory
from python_nedrex.common import http as _http


# pylint: disable=R0913
def must_request(
    seeds: _List[str],
    hubpenalty: float,
    multiple: bool,
    trees: int,
    maxit: int,
    network: str = "DEFAULT",
) -> str:
    body = {
        "seeds": seeds,
        "network": network,
        "hubpenalty": hubpenalty,
        "multiple": multiple,
        "trees": trees,
        "maxit": maxit,
    }

    url = f"{_config.url_base}/must/submit"
    resp = _http.post(url, json=body, headers={"x-api-key": _config.api_key})
    result: str = _check_response(resp)
    return result


# pylint: enable=R0913


check_must_status = _check_status_factory("/must/status")
