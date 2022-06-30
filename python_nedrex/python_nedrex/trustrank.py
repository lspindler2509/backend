from typing import List as _List
from typing import Optional as _Optional

from python_nedrex import config as _config
from python_nedrex.common import check_response as _check_response
from python_nedrex.common import check_status_factory as _check_status_factory
from python_nedrex.common import http as _http


def trustrank_submit(
    seeds: _List[str],
    damping_factor: float = 0.85,
    only_direct_drugs: bool = True,
    only_approved_drugs: bool = True,
    n: _Optional[int] = None,  # pylint: disable=C0103
) -> str:
    url = f"{_config.url_base}/trustrank/submit"

    body = {
        "seeds": seeds,
        "damping_factor": damping_factor,
        "only_direct_drugs": only_direct_drugs,
        "only_approved_drugs": only_approved_drugs,
        "N": n,
    }

    resp = _http.post(url, json=body, headers={"x-api-key": _config.api_key})
    result: str = _check_response(resp)
    return result


check_trustrank_status = _check_status_factory("/trustrank/status")


def download_trustrank_results(uid: str) -> str:
    url = f"{_config.url_base}/trustrank/download"
    params = {"uid": uid}

    resp = _http.get(url, params=params, headers={"x-api-key": _config.api_key})
    result: str = _check_response(resp, return_type="text")
    return result
