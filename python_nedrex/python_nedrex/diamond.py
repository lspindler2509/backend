from typing import List as _List

from python_nedrex import config as _config
from python_nedrex.common import check_response as _check_response
from python_nedrex.common import check_status_factory as _check_status_factory
from python_nedrex.common import http as _http


def diamond_submit(
    seeds: _List[str],
    n: int,  # pylint: disable=C0103
    alpha: int = 1,
    network: str = "DEFAULT",
    edges: str = "all",
) -> str:
    if edges not in {"limited", "all"}:
        raise ValueError(f"invalid value for argument edges ({edges!r}), should be all|limited")

    url = f"{_config.url_base}/diamond/submit"
    body = {
        "seeds": seeds,
        "n": n,
        "alpha": alpha,
        "network": network,
        "edges": edges,
    }

    resp = _http.post(url, json=body, headers={"x-api-key": _config.api_key})
    result: str = _check_response(resp)
    return result


check_diamond_status = _check_status_factory("/diamond/status")


def download_diamond_results(uid: str) -> str:
    url = f"{_config.url_base}/diamond/download"
    params = {"uid": uid}
    resp = _http.get(url, params=params, headers={"x-api-key": _config.api_key})
    result: str = _check_response(resp, return_type="text")
    return result
