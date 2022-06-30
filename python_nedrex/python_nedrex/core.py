"""Module containing functions relating to the general routes in the NeDRex API

Additionally, this module also contains routes for obtaining API keys.
"""

from typing import Any as _Any
from typing import Dict as _Dict
from typing import Generator as _Generator
from typing import List as _List
from typing import Optional as _Optional
from typing import cast as _cast

from python_nedrex import config as _config
from python_nedrex.common import check_pagination_limit as _check_pagination_limit
from python_nedrex.common import check_response as _check_response
from python_nedrex.common import get_pagination_limit as _get_pagination_limit
from python_nedrex.common import http as _http
from python_nedrex.decorators import check_url_base as _check_url_base
from python_nedrex.exceptions import NeDRexError as _NeDRexError


def _check_type(coll_name: str, coll_type: str) -> bool:
    if coll_type == "edge":
        if coll_name in get_edge_types():
            return True
        raise _NeDRexError(f"type={coll_name!r} not in NeDRex edge types")

    if coll_type == "node":
        if coll_name in get_node_types():
            return True
        raise _NeDRexError(f"type={coll_name!r} not in NeDRex node types")

    raise _NeDRexError(f"_check_type received invalid coll_type={coll_type!r}")


@_check_url_base
def api_keys_active() -> bool:
    """Checks whether API keys are active for the instance of NeDRex set in the config

    Returns True if the keys are active, False otherwise
    """
    url = f"{_config.url_base}/api_key_setting"
    response = _http.get(url)
    if response.status_code != 200:
        raise Exception("Unexpected non-200 status code")
    return _cast(bool, response.json())


@_check_url_base
def get_api_key(*, accept_eula: bool = False) -> _Any:
    """Obtains a new API key from the NeDRex API.

    This function will only return if accept_eula is explicitly set to True
    """
    if accept_eula is not True:
        raise _NeDRexError("an API key cannot be obtained unless accept_eula is set to True")

    url = f"{_config.url_base}/admin/api_key/generate"
    response = _http.post(url, json={"accept_eula": accept_eula})
    return response.json()


@_check_url_base
def get_node_types() -> _Any:
    """
    Returns the list of node types stored in NeDRexDB

        Returns:
            node_list (list[str]): List of node types in NeDRex
    """
    url: str = f"{_config.url_base}/list_node_collections"
    response = _http.get(url, headers={"x-api-key": _config.api_key})
    node_list = _check_response(response)
    return node_list


@_check_url_base
def get_edge_types() -> _Any:
    """
    Returns a list of edge types stored in NeDRexDB

        Returns:
            edge_list (list[str]): List of edge types in NeDRex
    """
    url: str = f"{_config.url_base}/list_edge_collections"
    response = _http.get(url, headers={"x-api-key": _config.api_key})
    edge_list = _check_response(response)
    return edge_list


@_check_url_base
def get_collection_attributes(coll_type: str, include_counts: bool = False) -> _Any:
    """
    Retrurns a list of the available attributes stored in NeDRex for the given type

        Parameters:
            type (str): The node or edge type to get available attributes for

        Returns:
            attr_list (list[str]): The list of available attributes for the specified node/edge type
    """
    url: str = f"{_config.url_base}/{coll_type}/attributes"
    response = _http.get(url, headers={"x-api-key": _config.api_key}, params={"include_counts": include_counts})
    attr_list = _check_response(response)
    return attr_list


@_check_url_base
def get_node_ids(coll_type: str) -> _Any:
    """
    Returns a list of node identifiers in NeDRex for the given type

        Parameters:
            type(str): The node type to get IDs for
        Returns:
            node_ids (list[str]): The list of available node_ids for the specified node type
    """
    _check_type(coll_type, "node")

    url: str = f"{_config.url_base}/{coll_type}/attributes/primaryDomainId/json"

    resp = _http.get(url, headers={"x-api-key": _config.api_key})
    data = _check_response(resp)
    node_ids = [i["primaryDomainId"] for i in data]
    return node_ids


@_check_url_base
def get_nodes(
    node_type: str,
    attributes: _Optional[_List[str]] = None,
    node_ids: _Optional[_List[str]] = None,
    limit: _Optional[int] = None,
    offset: int = 0,
) -> _Any:
    """
    Returns nodes in NeDRex for the given type

        Parameters:
            node_type (str): The node type to collect
            attributes (Optional[list[str]]): A list of attributes to return for the collected nodes. The default
              (None) returns all attributes.
            node_ids (Optional[list[str]]): A list of the specific node IDs to return. The default (None) returns all
              nodes.
            limit (Optional[int]): A limit for the number of records to be returned. The maximum value for this is set
              by the API.
            offset (int): The number of records to skip before returning records. Default is 0 (no records skipped).
        Returns:
            node_ids (list[str]): The list of available node_ids for the specified node type
    """
    _check_type(node_type, "node")

    upper_limit = _get_pagination_limit()
    _check_pagination_limit(limit, upper_limit)

    params = {"node_id": node_ids, "attribute": attributes, "offset": offset, "limit": limit}

    resp = _http.get(
        f"{_config.url_base}/{node_type}/attributes/json", params=params, headers={"x-api-key": _config.api_key}
    )

    items = _check_response(resp)
    return items


@_check_url_base
def iter_nodes(
    node_type: str, attributes: _Optional[_List[str]] = None, node_ids: _Optional[_List[str]] = None
) -> _Generator[_Dict[str, _Any], None, None]:

    _check_type(node_type, "node")
    upper_limit = _get_pagination_limit()

    params: _Dict[str, _Any] = {"node_id": node_ids, "attribute": attributes, "limit": upper_limit}

    offset = 0
    while True:
        params["offset"] = offset
        resp = _http.get(
            f"{_config.url_base}/{node_type}/attributes/json", params=params, headers={"x-api-key": _config.api_key}
        )

        data = _check_response(resp)

        for doc in data:
            yield doc

        if len(data) < upper_limit:
            break
        offset += upper_limit


@_check_url_base
def get_edges(edge_type: str, limit: _Optional[int] = None, offset: _Optional[int] = None) -> _Any:
    """
    Returns edges in NeDRex of the given type

        Parameters:
            edge_type (str): The node type to collect
            limit (Optional[int]): A limit for the number of records to be returned. The maximum value for this is set
              by the API.
            offset (int): The number of records to skip before returning records. Default is 0 (no records skipped).
    """
    _check_type(edge_type, "edge")

    params = {"limit": limit, "offset": offset, "api_key": _config.api_key}

    resp = _http.get(f"{_config.url_base}/{edge_type}/all", params=params, headers={"x-api-key": _config.api_key})
    items = _check_response(resp)
    return items


@_check_url_base
def iter_edges(edge_type: str) -> _Generator[_Dict[str, _Any], None, None]:
    _check_type(edge_type, "edge")
    upper_limit = _get_pagination_limit()

    offset = 0
    while True:
        params = {"offset": offset, "limit": upper_limit}
        resp = _http.get(f"{_config.url_base}/{edge_type}/all", params=params, headers={"x-api-key": _config.api_key})
        data = _check_response(resp)

        for doc in data:
            yield doc

        if len(data) < upper_limit:
            break
        offset += upper_limit
