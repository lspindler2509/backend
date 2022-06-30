"""Module containing python functions to access the disorder routes in the NeDRex API"""

from typing import Any as _Any
from typing import Callable as _Callable
from typing import List as _List
from typing import Union as _Union

from python_nedrex import config as _config
from python_nedrex.common import check_response as _check_response
from python_nedrex.common import http as _http
from python_nedrex.decorators import check_url_base as _check_url_base


def _generate_route(path: str) -> _Callable[[_Union[str, _List[str]]], _Any]:
    @_check_url_base
    def new_func(codes: _Union[str, _List[str]]) -> _Any:
        if isinstance(codes, str):
            codes = [codes]

        url = f"{_config.url_base}/disorder/{path}"
        resp = _http.get(url, params={"q": codes}, headers={"x-api-key": _config.api_key})
        return _check_response(resp)

    return new_func


search_by_icd10 = _generate_route("get_by_icd10")
get_disorder_descendants = _generate_route("descendants")
get_disorder_ancestors = _generate_route("ancestors")
get_disorder_parents = _generate_route("parents")
get_disorder_children = _generate_route("children")
