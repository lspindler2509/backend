from typing import Any as _Any
from typing import Dict as _Dict
from typing import Generator as _Generator
from typing import List as _List
from typing import Optional as _Optional

from python_nedrex import config as _config
from python_nedrex.common import check_pagination_limit as _check_pagination_limit
from python_nedrex.common import check_response as _check_response
from python_nedrex.common import get_pagination_limit as _get_pagination_limit
from python_nedrex.common import http as _http


def get_effect_choices() -> _List[str]:
    url = f"{_config.url_base}/variants/get_effect_choices"
    resp = _http.get(url, headers={"x-api-key": _config.api_key})
    result: _List[str] = _check_response(resp)
    return result


def get_review_status_choices() -> _List[str]:
    url = f"{_config.url_base}/variants/get_review_choices"
    resp = _http.get(url, headers={"x-api-key": _config.api_key})
    result: _List[str] = _check_response(resp)
    return result


# pylint: disable=R0913
def get_variant_disorder_associations(
    variant_ids: _Optional[_List[str]] = None,
    disorder_ids: _Optional[_List[str]] = None,
    review_status: _Optional[_List[str]] = None,
    effect: _Optional[_List[str]] = None,
    limit: _Optional[int] = None,
    offset: int = 0,
) -> _List[_Dict[str, _Any]]:
    max_limit = _get_pagination_limit()
    if isinstance(limit, int):
        _check_pagination_limit(limit, max_limit)
    else:
        limit = max_limit

    params = {
        "variant_id": variant_ids,
        "disorder_id": disorder_ids,
        "review_status": review_status,
        "effect": effect,
        "limit": limit,
        "offset": offset,
    }

    url = f"{_config.url_base}/variants/get_variant_disorder_associations"
    resp = _http.get(url, params=params, headers={"x-api-key": _config.api_key})
    result: _List[_Dict[str, _Any]] = _check_response(resp)
    return result


# pylint: enable=R0913


def iter_variant_disorder_associations(
    variant_ids: _Optional[_List[str]] = None,
    disorder_ids: _Optional[_List[str]] = None,
    review_status: _Optional[_List[str]] = None,
    effect: _Optional[_List[str]] = None,
) -> _Generator[_Dict[str, _Any], None, None]:
    max_limit = _get_pagination_limit()
    offset = 0

    kwargs = {
        "variant_ids": variant_ids,
        "disorder_ids": disorder_ids,
        "review_status": review_status,
        "effect": effect,
        "limit": max_limit,
        "offset": offset,
    }

    while True:
        results = get_variant_disorder_associations(**kwargs)
        if len(results) == 0:
            return
        yield from results

        offset += max_limit
        kwargs["offset"] = offset


def get_variant_gene_associations(
    variant_ids: _Optional[_List[str]] = None,
    gene_ids: _Optional[_List[str]] = None,
    limit: _Optional[int] = None,
    offset: int = 0,
) -> _List[_Dict[str, _Any]]:

    max_limit = _get_pagination_limit()
    if isinstance(limit, int):
        _check_pagination_limit(limit, max_limit)
    else:
        limit = max_limit

    params = {
        "variant_id": variant_ids,
        "gene_id": gene_ids,
        "offset": offset,
        "limit": limit,
    }

    url = f"{_config.url_base}/variants/get_variant_gene_associations"
    resp = _http.get(url, params=params, headers={"x-api-key": _config.api_key})
    result: _List[_Dict[str, _Any]] = _check_response(resp)
    return result


def iter_variant_gene_associations(
    variant_ids: _Optional[_List[str]] = None,
    gene_ids: _Optional[_List[str]] = None,
) -> _Generator[_Dict[str, _Any], None, None]:
    max_limit = _get_pagination_limit()
    offset = 0

    kwargs = {
        "variant_ids": variant_ids,
        "gene_ids": gene_ids,
        "limit": max_limit,
        "offset": offset,
    }

    while True:
        results = get_variant_gene_associations(**kwargs)
        if len(results) == 0:
            return
        yield from results

        offset += max_limit
        kwargs["offset"] = offset


def get_variant_based_disorder_associated_genes(
    disorder_id: str, review_status: _Optional[_List[str]] = None, effect: _Optional[_List[str]] = None
) -> _List[_Dict[str, _Any]]:
    params = {"disorder_id": disorder_id, "review_status": review_status, "effect": effect}

    url = f"{_config.url_base}/variants/variant_based_disorder_associated_genes"

    resp = _http.get(url, params=params, headers={"x-api-key": _config.api_key})
    result: _List[_Dict[str, _Any]] = _check_response(resp)
    return result


def get_variant_based_gene_associated_disorders(
    gene_id: str, review_status: _Optional[_List[str]] = None, effect: _Optional[_List[str]] = None
) -> _List[_Dict[str, _Any]]:
    params = {"gene_id": gene_id, "review_status": review_status, "effect": effect}

    url = f"{_config.url_base}/variants/variant_based_gene_associated_disorders"

    resp = _http.get(url, params=params, headers={"x-api-key": _config.api_key})
    result: _List[_Dict[str, _Any]] = _check_response(resp)
    return result
