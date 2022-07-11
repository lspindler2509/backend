from typing import Any as _Any
from typing import Dict as _Dict
from typing import Optional as _Optional

from python_nedrex import config as _config
from python_nedrex.common import check_response as _check_response
from python_nedrex.common import download_file as _download_file
from python_nedrex.common import http as _http


def get_metadata() -> _Dict[str, _Any]:
    url = f"{_config.url_base}/static/metadata"
    resp = _http.get(url, headers={"x-api-key": _config.api_key})
    result: _Dict[str, _Any] = _check_response(resp)
    return result


def get_license() -> str:
    url = f"{_config.url_base}/static/licence"
    resp = _http.get(url)
    result: str = _check_response(resp=resp, return_type='text')
    return result


def download_lengths_map(target: _Optional[str] = None) -> None:
    if target is None:
        target = "lengths.map"

    url = f"{_config.url_base}/static/lengths.map"

    _download_file(url, target)
