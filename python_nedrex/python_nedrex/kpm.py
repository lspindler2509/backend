from typing import Any as _Any
from typing import Dict as _Dict
from typing import List as _List

from python_nedrex import config as _config
from python_nedrex.common import check_response as _check_response
from python_nedrex.common import http as _http


def kpm_submit(seeds: _List[str], k: int, network: str = "DEFAULT") -> str:

    url = f"{_config.url_base}/kpm/submit"
    body = {"seeds": seeds, "k": k, "network": network}

    resp = _http.post(url, json=body, headers={"x-api-key": _config.api_key})
    result: str = _check_response(resp)
    return result


def check_kpm_status(uid: str) -> _Dict[str, _Any]:

    url = f"{_config.url_base}/kpm/status"
    params = {"uid": uid}

    resp = _http.get(url, params=params, headers={"x-api-key": _config.api_key})
    result: _Dict[str, _Any] = _check_response(resp)
    return result
