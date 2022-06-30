from typing import List as _List

from python_nedrex import config as _config
from python_nedrex.common import check_response as _check_response
from python_nedrex.common import check_status_factory as _check_status_factory
from python_nedrex.common import http as _http


# pylint: disable=R0913
def robust_submit(
    seeds: _List[str],
    network: str = "DEFAULT",
    initial_fraction: float = 0.25,
    reduction_factor: float = 0.9,
    num_trees: int = 30,
    threshold: float = 0.1,
) -> str:

    body = {
        "seeds": seeds,
        "network": network,
        "initial_fraction": initial_fraction,
        "reduction_factor": reduction_factor,
        "num_trees": num_trees,
        "threshold": threshold,
    }
    url = f"{_config.url_base}/robust/submit"

    resp = _http.post(url, json=body, headers={"x-api-key": _config.api_key})
    result: str = _check_response(resp)
    return result


# pylint: enable=R0913

check_robust_status = _check_status_factory("/robust/status")


def download_robust_results(uid: str) -> str:
    url = f"{_config.url_base}/robust/results"
    params = {"uid": uid}

    resp = _http.get(url, params=params, headers={"x-api-key": _config.api_key})
    result: str = _check_response(resp, return_type="text")
    return result
