from typing import Any as _Any
from typing import Dict as _Dict
from typing import Iterable as _Iterable
from typing import List as _List
from typing import Optional as _Optional

from python_nedrex import config as _config
from python_nedrex.common import check_pagination_limit as _check_pagination_limit
from python_nedrex.common import check_response as _check_response
from python_nedrex.common import get_pagination_limit as _get_pagination_limit
from python_nedrex.common import http as _http
from python_nedrex.exceptions import NeDRexError


def ppis(evidence: _Iterable[str], skip: int = 0, limit: _Optional[int] = None) -> _List[_Dict[str, _Any]]:
    evidence_set = set(evidence)
    extra_evidence = evidence_set - {"exp", "pred", "ortho"}
    if extra_evidence:
        raise NeDRexError(f"unexpected evidence types: {extra_evidence}")

    maximum_limit = _get_pagination_limit()
    _check_pagination_limit(limit, maximum_limit)

    params = {"iid_evidence": list(evidence_set), "skip": skip, "limit": limit}

    resp = _http.get(f"{_config.url_base}/ppi", params=params, headers={"x-api-key": _config.api_key})
    result: _List[_Dict[str, _Any]] = _check_response(resp)
    return result
