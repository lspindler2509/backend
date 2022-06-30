from typing import List as _List

from python_nedrex import config as _config
from python_nedrex.common import check_response as _check_response
from python_nedrex.common import check_status_factory as _check_status_factory
from python_nedrex.common import http as _http

check_validation_status = _check_status_factory("/validation/status")


def check_module_member_type(mmt: str) -> None:
    if mmt not in {"gene", "protein"}:
        raise ValueError(f"module_member_type {mmt!r} is invalid (should be 'gene' or 'protein'")


# pylint: disable=R0913
def joint_validation_submit(
    module_members: _List[str],
    module_member_type: str,
    test_drugs: _List[str],
    true_drugs: _List[str],
    permutations: int,
    only_approved_drugs: bool = True,
) -> str:
    check_module_member_type(module_member_type)

    url = f"{_config.url_base}/validation/joint"
    body = {
        "module_members": module_members,
        "module_member_type": module_member_type,
        "test_drugs": test_drugs,
        "true_drugs": true_drugs,
        "permutations": permutations,
        "only_approved_drugs": only_approved_drugs,
    }
    resp = _http.post(url, json=body, headers={"x-api-key": _config.api_key})
    result: str = _check_response(resp)
    return result


# pylint: enable=R0913


def module_validation_submit(
    module_members: _List[str],
    module_member_type: str,
    true_drugs: _List[str],
    permutations: int,
    only_approved_drugs: bool = True,
) -> str:
    check_module_member_type(module_member_type)

    url = f"{_config.url_base}/validation/module"
    body = {
        "module_members": module_members,
        "module_member_type": module_member_type,
        "true_drugs": true_drugs,
        "permutations": permutations,
        "only_approved_drugs": only_approved_drugs,
    }
    resp = _http.post(url, json=body, headers={"x-api-key": _config.api_key})
    result: str = _check_response(resp)
    return result


def drug_validation_submit(
    test_drugs: _List[str], true_drugs: _List[str], permutations: int, only_approved_drugs: bool = True
) -> str:

    url = f"{_config.url_base}/validation/drug"
    body = {
        "test_drugs": test_drugs,
        "true_drugs": true_drugs,
        "permutations": permutations,
        "only_approved_drugs": only_approved_drugs,
    }

    resp = _http.post(url, json=body, headers={"x-api-key": _config.api_key})
    result: str = _check_response(resp)
    return result
